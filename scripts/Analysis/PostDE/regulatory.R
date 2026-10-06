postde_read_network <- function(path, source_col = "source", target_col = "target", mor_col = "mor") {
    if (is.null(source_col)) source_col <- "source"
    if (is.null(target_col)) target_col <- "target"
    if (is.null(mor_col)) mor_col <- "mor"
    df <- postde_read_tsv(path, header = TRUE)
    if (!all(c(source_col, target_col, mor_col) %in% colnames(df))) {
        stop(paste0("network file must have columns ", source_col, ", ", target_col, ", ", mor_col))
    }
    df <- df[, c(source_col, target_col, mor_col)]
    colnames(df) <- c("source", "target", "mor")
    df
}

postde_validate_network <- function(net) {
    net <- as.data.frame(net)
    if (!all(c("source", "target", "mor") %in% colnames(net))) {
        stop("network must have source, target and mor columns")
    }
    net$source <- as.character(net$source)
    net$target <- as.character(net$target)
    net$mor <- as.numeric(net$mor)
    if (any(is.na(net$source)) || any(!nzchar(net$source)) || any(is.na(net$target)) || any(!nzchar(net$target))) {
        stop("network contains empty or NA source/target ids")
    }
    if (any(!is.finite(net$mor))) {
        stop("network contains non-finite mor values")
    }
    if (any(net$mor == 0)) {
        stop("network contains zero mor values")
    }
    key <- paste(net$source, net$target, sep = "\r")
    if (anyDuplicated(key)) {
        stop(paste0("network contains duplicate source-target pairs: ", paste(unique(net$source[duplicated(key)]), collapse = ", ")))
    }
    net
}

postde_check_mlm_rank <- function(net, expr) {
    net <- as.data.frame(net)
    expr <- as.matrix(expr)
    sources <- unique(net$source)
    if (length(sources) == 0) {
        stop("no sources in network")
    }
    mor_mat <- matrix(0, nrow = nrow(expr), ncol = length(sources), dimnames = list(rownames(expr), sources))
    for (i in seq_len(nrow(net))) {
        if (net$target[i] %in% rownames(expr)) {
            mor_mat[net$target[i], net$source[i]] <- net$mor[i]
        }
    }
    design <- cbind(1, mor_mat)
    r <- qr(design)$rank
    if (r < ncol(design)) {
        stop(paste0("MLM design is rank deficient (rank ", r, " < ", ncol(design), "): sources are collinear with the intercept or each other; remove collinear sources"))
    }
    if (nrow(design) - r <= 0) {
        stop("MLM design has no residual degrees of freedom")
    }
    invisible(TRUE)
}

postde_activity_matrix <- function(acts, sample_names) {
    acts <- as.data.frame(acts)
    if (!all(c("source", "condition", "score") %in% colnames(acts))) {
        stop("decoupler result must contain source, condition and score columns")
    }
    if (any(!is.finite(acts$score))) {
        stop("decoupler activity scores contain non-finite values")
    }
    sources <- unique(acts$source)
    mat <- matrix(NA_real_, nrow = length(sources), ncol = length(sample_names), dimnames = list(sources, sample_names))
    for (i in seq_len(nrow(acts))) {
        if (acts$condition[i] %in% sample_names) {
            mat[acts$source[i], acts$condition[i]] <- acts$score[i]
        }
    }
    if (any(is.na(mat))) {
        stop("decoupler activity matrix has missing values after pivoting")
    }
    mat
}

postde_run_decoupler <- function(entry, config, outdir) {
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        stop("decoupler enabled but package decoupleR is not installed")
    }
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("decoupler enabled but package limma is not installed")
    }
    dc <- config$decoupler
    net <- NULL
    if (!is.null(dc$network)) {
        net <- postde_read_network(dc$network, dc$source_col, dc$target_col, dc$mor_col)
    } else {
        if (!isTRUE(dc$allow_network)) {
            stop("decoupler network retrieval requires allow_network=true")
        }
        resource <- if (is.null(dc$resource)) "collectri" else dc$resource
        organism <- if (is.null(dc$organism)) "human" else dc$organism
        if (resource == "collectri") {
            net <- decoupleR::get_collectri(organism = organism)
        } else if (resource == "progeny") {
            top <- if (is.null(dc$top)) 500 else dc$top
            net <- decoupleR::get_progeny(organism = organism, top = top)
            if ("weight" %in% colnames(net) && !"mor" %in% colnames(net)) {
                net$mor <- net$weight
            }
        } else {
            stop(paste0("unknown decoupler resource '", resource, "' (use collectri or progeny)"))
        }
        net <- as.data.frame(net)
        if (!all(c("source", "target", "mor") %in% colnames(net))) {
            stop(paste0("retrieved resource '", resource, "' lacks source/target/mor columns"))
        }
        net <- net[, c("source", "target", "mor")]
    }
    net <- postde_validate_network(net)
    saveRDS(net, file.path(outdir, "decoupler_network.rds"))
    targets <- intersect(net$target, rownames(entry$expression))
    if (length(targets) == 0) {
        stop("no network targets present in expression matrix")
    }
    net <- net[net$target %in% targets, , drop = FALSE]
    net <- decoupleR::rename_net(net, source, target, mor)
    min_size <- if (is.null(dc$min_size)) config$min_size else dc$min_size
    net <- decoupleR::filt_minsize(rownames(entry$expression), net, minsize = min_size)
    if (nrow(net) == 0) {
        stop("no sources with at least min_size targets after filtering")
    }
    write.table(net, file.path(outdir, "decoupler_network_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    net_meta <- list(
        coverage = length(intersect(net$target, rownames(entry$expression))) / length(unique(net$target)),
        checksum = unname(tools::md5sum(file.path(outdir, "decoupler_network_filtered.tsv"))),
        package = as.character(utils::packageVersion("decoupleR"))
    )
    saveRDS(net_meta, file.path(outdir, "decoupler_network_meta.rds"))
    methods <- unlist(dc$methods)
    if (is.null(methods)) {
        methods <- c("ulm", "mlm")
    }
    out <- list()
    for (method in methods) {
        if (method == "ulm") {
            acts <- decoupleR::run_ulm(entry$expression, net, minsize = min_size)
        } else if (method == "mlm") {
            postde_check_mlm_rank(net, entry$expression)
            acts <- decoupleR::run_mlm(entry$expression, net, minsize = min_size)
        } else {
            stop(paste0("unknown decoupler method '", method, "' (use ulm or mlm)"))
        }
        act_mat <- postde_activity_matrix(acts, rownames(entry$design))
        write.table(act_mat, file.path(outdir, paste0("decoupler_", method, "_activities.tsv")), sep = "\t", col.names = NA, quote = FALSE)
        diff <- postde_differential_scores(act_mat, entry$design, entry$contrast)
        write.table(diff$results, file.path(outdir, paste0("decoupler_", method, "_differential.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(diff$design, file.path(outdir, "decoupler_design.tsv"), sep = "\t", col.names = NA, quote = FALSE)
        write.table(data.frame(term = names(diff$contrast), coefficient = diff$contrast), file.path(outdir, "decoupler_contrast.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
        out[[method]] <- list(status = "ok", activities = file.path(outdir, paste0("decoupler_", method, "_activities.tsv")), differential = file.path(outdir, paste0("decoupler_", method, "_differential.tsv")))
        if (isTRUE(dc$contrast_activity)) {
            ca <- data.frame(
                source = diff$results$feature_id,
                contrast_activity = diff$results$stat,
                pvalue = diff$results$pvalue,
                padj = diff$results$padj,
                stringsAsFactors = FALSE
            )
            write.table(ca, file.path(outdir, paste0("decoupler_", method, "_contrast_activity.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
            out[[method]]$contrast_activity <- file.path(outdir, paste0("decoupler_", method, "_contrast_activity.tsv"))
        }
    }
    out
}

postde_safe_formula <- function(formula) {
    if (!inherits(formula, "formula")) {
        stop("dream formula must be a formula")
    }
    allowed_ops <- c("~", "+", "-", "*", ":", "/", "|", "||", "(")
    walk <- function(e) {
        if (is.name(e) || is.numeric(e)) {
            return(invisible(TRUE))
        }
        if (is.call(e)) {
            op <- as.character(e[[1]])
            if (!op %in% allowed_ops) {
                stop(paste0("dream formula contains disallowed operator/call '", op, "'"))
            }
            for (i in 2:length(e)) {
                walk(e[[i]])
            }
            return(invisible(TRUE))
        }
        stop(paste0("dream formula contains disallowed element: ", deparse(e)))
    }
    walk(formula)
    invisible(TRUE)
}

postde_safe_contrast_expr <- function(expr) {
    if (!is.character(expr) || length(expr) != 1 || !nzchar(expr)) {
        stop("dream contrast must be a single non-empty string expression")
    }
    parsed <- tryCatch(parse(text = expr), error = function(err) {
        stop(paste0("dream contrast is not parseable: ", conditionMessage(err)))
    })
    allowed_ops <- c("+", "-", "*", "/", ":")
    walk <- function(e) {
        if (is.name(e) || is.numeric(e)) {
            return(invisible(TRUE))
        }
        if (is.call(e)) {
            op <- as.character(e[[1]])
            if (!op %in% allowed_ops) {
                stop(paste0("dream contrast contains disallowed call '", op, "'"))
            }
            for (i in 2:length(e)) {
                walk(e[[i]])
            }
            return(invisible(TRUE))
        }
        stop(paste0("dream contrast contains disallowed element: ", deparse(e)))
    }
    for (i in seq_along(parsed)) {
        walk(parsed[[i]])
    }
    invisible(TRUE)
}

postde_fixed_formula <- function(formula) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("lme4 is required for dream but is not installed")
    }
    fixed <- suppressWarnings(lme4::nobars(formula))
    tl <- attr(terms(fixed), "term.labels")
    if (length(tl) == 0 && attr(terms(fixed), "intercept") == 0) {
        stop("dream formula has no fixed effects")
    }
    fixed
}

postde_random_group <- function(formula) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("lme4 is required for dream but is not installed")
    }
    bars <- suppressWarnings(lme4::findbars(formula))
    if (length(bars) == 0) {
        stop("dream formula must include a random effect term like (1|subject)")
    }
    groups <- vapply(bars, function(b) {
        gsub(".*\\|\\s*", "", deparse(b))
    }, character(1))
    unique(groups)
}

postde_run_dream <- function(entry, config, outdir) {
    if (isTRUE(entry$normalized)) {
        return(list(status = "skipped: normalized entry; use the ordinary entry"))
    }
    if (!requireNamespace("variancePartition", quietly = TRUE)) {
        stop("dream enabled but package variancePartition is not installed")
    }
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("dream enabled but package limma is not installed")
    }
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("dream enabled but package edgeR is not installed")
    }
    dm <- config$dream
    if (is.null(dm$metadata)) {
        stop("dream enabled but metadata path missing")
    }
    if (is.null(dm$formula)) {
        stop("dream enabled but formula missing")
    }
    meta <- postde_read_tsv(dm$metadata, header = TRUE)
    if (!"sample" %in% colnames(meta)) {
        stop("dream metadata must contain a 'sample' column")
    }
    if (anyDuplicated(meta$sample)) {
        stop("dream metadata contains duplicate sample ids")
    }
    if ("condition" %in% colnames(meta)) {
        stop("dream metadata must not contain a 'condition' column (condition comes from the DE analysis)")
    }
    if (!all(rownames(entry$metadata) %in% meta$sample)) {
        stop("dream metadata is missing samples present in the entry")
    }
    entry_meta <- as.data.frame(entry$metadata)
    entry_meta$sample <- rownames(entry_meta)
    common_cols <- setdiff(intersect(colnames(meta), colnames(entry_meta)), "sample")
    for (cc in common_cols) {
        mm <- merge(entry_meta[, c("sample", cc)], meta[, c("sample", cc)], by = "sample", suffixes = c(".entry", ".dream"))
        inconsistent <- !is.na(mm[[paste0(cc, ".entry")]]) & !is.na(mm[[paste0(cc, ".dream")]]) & as.character(mm[[paste0(cc, ".entry")]]) != as.character(mm[[paste0(cc, ".dream")]])
        if (any(inconsistent)) {
            stop(paste0("dream metadata column '", cc, "' has values inconsistent with the entry metadata"))
        }
    }
    if (length(common_cols) > 0) {
        meta <- meta[, setdiff(colnames(meta), common_cols), drop = FALSE]
    }
    data <- merge(entry_meta, meta, by = "sample", all.x = TRUE)
    rownames(data) <- data$sample
    data <- data[colnames(entry$counts), , drop = FALSE]
    formula <- as.formula(dm$formula)
    postde_safe_formula(formula)
    groups <- postde_random_group(formula)
    for (group in groups) {
        if (!group %in% colnames(data)) {
            stop(paste0("dream random grouping variable '", group, "' not found in metadata"))
        }
        tab <- table(data[[group]])
        if (sum(tab >= 2) < 2) {
            stop(paste0("dream random grouping variable '", group, "' must have at least 2 groups with more than 1 sample each"))
        }
    }
    fixed <- postde_fixed_formula(formula)
    fixed_design <- model.matrix(fixed, data = data)
    r <- qr(fixed_design)$rank
    if (r < ncol(fixed_design)) {
        stop("dream fixed design is not full rank")
    }
    if (nrow(fixed_design) - r <= 0) {
        stop("dream fixed design has no residual degrees of freedom")
    }
    dge <- edgeR::DGEList(counts = entry$counts)
    keep <- edgeR::filterByExpr(dge, design = fixed_design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    vobj <- variancePartition::voomWithDreamWeights(dge, formula, data)
    if (!is.null(dm$contrasts)) {
        contrasts <- unlist(dm$contrasts)
        for (nm in names(contrasts)) {
            postde_safe_contrast_expr(contrasts[[nm]])
        }
        L <- variancePartition::makeContrastsDream(formula, data, contrasts = contrasts)
    } else {
        if (!"condition" %in% colnames(data)) {
            stop("dream default contrast A-B requires 'condition' in the formula/metadata")
        }
        cvec <- postde_contrast(fixed, data, entry$A, entry$B)
        if (!identical(names(cvec), colnames(fixed_design))) {
            stop("dream default contrast names do not match fixed design columns")
        }
        L <- matrix(cvec, ncol = 1, dimnames = list(names(cvec), "A-B"))
    }
    fit <- variancePartition::dream(vobj, formula, data, L = L)
    fit <- variancePartition::eBayes(fit)
    out <- list()
    if (!is.null(dm$contrasts)) {
        for (nm in colnames(L)) {
            tt <- limma::topTable(fit, number = Inf, sort.by = "none", coef = nm)
            res <- data.frame(
                gene_id = rownames(tt),
                logFC = tt$logFC,
                pvalue = tt$P.Value,
                padj = tt$adj.P.Val,
                stat = tt$t,
                mean = tt$AveExpr,
                stringsAsFactors = FALSE
            )
            res_file <- file.path(outdir, paste0("dream_result_", postde_sanitize_id(nm), ".tsv"))
            write.table(res, res_file, sep = "\t", row.names = FALSE, quote = FALSE)
            out[[nm]] <- list(status = "ok", result = res_file)
        }
    } else {
        tt <- limma::topTable(fit, number = Inf, sort.by = "none", coef = "A-B")
        res <- data.frame(
            gene_id = rownames(tt),
            logFC = tt$logFC,
            pvalue = tt$P.Value,
            padj = tt$adj.P.Val,
            stat = tt$t,
            mean = tt$AveExpr,
            stringsAsFactors = FALSE
        )
        res_file <- file.path(outdir, "dream_result.tsv")
        write.table(res, res_file, sep = "\t", row.names = FALSE, quote = FALSE)
        out[["A-B"]] <- list(status = "ok", result = res_file)
    }
    write.table(fixed_design, file.path(outdir, "dream_design.tsv"), sep = "\t", col.names = NA, quote = FALSE)
    writeLines(dm$formula, file.path(outdir, "dream_formula.txt"))
    writeLines(rownames(dge), file.path(outdir, "dream_filtered_genes.txt"))
    saveRDS(L, file.path(outdir, "dream_contrast.rds"))
    varpart <- tryCatch(
        variancePartition::fitExtractVarPartModel(vobj, formula, data),
        error = function(e) structure(conditionMessage(e), class = "postde_varpart_error")
    )
    if (inherits(varpart, "postde_varpart_error")) {
        writeLines(as.character(varpart), file.path(outdir, "dream_varpart_status.txt"))
        out$varpart <- list(status = "error", message = as.character(varpart))
    } else {
        saveRDS(varpart, file.path(outdir, "dream_varpart.rds"))
        out$varpart <- list(status = "ok", result = file.path(outdir, "dream_varpart.rds"))
    }
    out
}

postde_run_regulatory <- function(entry, config, outdir) {
    out <- list()
    if (is.list(config$decoupler) && isTRUE(config$decoupler$enabled)) {
        out$decoupler <- postde_run_decoupler(entry, config, outdir)
    }
    if (is.list(config$dream) && isTRUE(config$dream$enabled)) {
        out$dream <- postde_run_dream(entry, config, outdir)
    }
    out
}
