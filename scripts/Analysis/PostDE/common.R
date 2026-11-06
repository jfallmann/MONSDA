postde_validate_design <- function(design, contrast) {
    design <- as.matrix(design)
    if (nrow(design) < 2) {
        stop("design has fewer than 2 rows")
    }
    if (any(!is.finite(design))) {
        stop("design contains non-finite values")
    }
    if (is.null(colnames(design))) {
        stop("design has no column names")
    }
    if (is.null(names(contrast))) {
        contrast <- as.numeric(contrast)
    } else {
        if (!all(names(contrast) %in% colnames(design))) {
            stop(paste0("contrast names do not match design columns: ", paste(setdiff(names(contrast), colnames(design)), collapse = ", ")))
        }
        contrast <- as.numeric(contrast[colnames(design)])
    }
    if (length(contrast) != ncol(design)) {
        stop(paste0("contrast length (", length(contrast), ") does not match design columns (", ncol(design), ")"))
    }
    if (any(!is.finite(contrast))) {
        stop("contrast contains non-finite values")
    }
    if (all(contrast == 0)) {
        stop("contrast is all zero")
    }
    const_cols <- apply(design, 2, function(x) length(unique(x)) == 1)
    if (any(const_cols & contrast != 0)) {
        stop(paste0("contrast has nonzero coefficient on constant design column(s): ", paste(colnames(design)[const_cols & contrast != 0], collapse = ", ")))
    }
    r <- qr(design)$rank
    if (r < ncol(design)) {
        stop(paste0("design is not full rank (rank ", r, " < ", ncol(design), " columns)"))
    }
    if (nrow(design) - r <= 0) {
        stop("design has no residual degrees of freedom")
    }
    a <- design %*% solve(crossprod(design)) %*% contrast
    fitted <- crossprod(design, a)
    if (max(abs(fitted - contrast)) > 1e-8) {
        stop("contrast is not estimable from the design")
    }
    invisible(TRUE)
}

postde_contrast <- function(formula, metadata, A, B) {
    md <- as.data.frame(metadata)
    if (!"condition" %in% colnames(md)) {
        stop("postde_contrast: metadata has no condition column")
    }
    md$condition <- relevel(md$condition, ref = B)
    mdA <- md
    mdB <- md
    mdA$condition <- factor(rep(A, nrow(md)), levels = levels(md$condition))
    mdB$condition <- factor(rep(B, nrow(md)), levels = levels(md$condition))
    mA <- model.matrix(formula, data = mdA)
    mB <- model.matrix(formula, data = mdB)
    setNames(colMeans(mA) - colMeans(mB), colnames(mA))
}

postde_rank_genes <- function(results, stat_col = "stat") {
    if (!is.data.frame(results)) {
        results <- as.data.frame(results)
    }
    if (!"gene_id" %in% colnames(results)) {
        stop("results must contain a gene_id column")
    }
    if (!stat_col %in% colnames(results)) {
        stop(paste0("results must contain a ", stat_col, " column"))
    }
    if (any(is.na(results$gene_id)) || any(!nzchar(as.character(results$gene_id)))) {
        stop("results contain empty or NA gene ids")
    }
    if (anyDuplicated(results$gene_id)) {
        stop("results contain duplicate gene ids")
    }
    if (!is.numeric(results[[stat_col]])) {
        stop(paste0(stat_col, " column must be numeric"))
    }
    if (any(!is.finite(results[[stat_col]]))) {
        results <- results[is.finite(results[[stat_col]]), , drop = FALSE]
    }
    if (nrow(results) == 0) {
        stop("results contain no genes with finite ", stat_col)
    }
    results <- results[order(-results[[stat_col]], results$gene_id), , drop = FALSE]
    results$rank <- seq_len(nrow(results))
    results
}

postde_ora_universe <- function(results, custom = NULL) {
    if (!is.data.frame(results)) {
        results <- as.data.frame(results)
    }
    if (!"gene_id" %in% colnames(results)) {
        stop("results must contain a gene_id column")
    }
    if (!"pvalue" %in% colnames(results)) {
        stop("results must contain a pvalue column")
    }
    if (any(is.na(results$gene_id)) || any(!nzchar(as.character(results$gene_id)))) {
        stop("results contain empty or NA gene ids")
    }
    tested <- results$gene_id[is.finite(results$pvalue)]
    if (is.null(custom)) {
        unique(tested)
    } else {
        unique(intersect(custom, tested))
    }
}

postde_significant <- function(results, direction = "all", padj = 0.05, lfc = 1) {
    if (!is.data.frame(results)) {
        results <- as.data.frame(results)
    }
    if (!all(c("gene_id", "logFC", "pvalue", "padj") %in% colnames(results))) {
        stop("results must contain gene_id, logFC, pvalue and padj columns")
    }
    if (any(is.na(results$gene_id)) || any(!nzchar(as.character(results$gene_id)))) {
        stop("results contain empty or NA gene ids")
    }
    sig <- is.finite(results$padj) & results$padj < padj
    if (direction == "up") {
        sig <- sig & results$logFC >= lfc
    } else if (direction == "down") {
        sig <- sig & results$logFC <= -lfc
    } else if (direction == "all") {
        sig <- sig & abs(results$logFC) >= lfc
    } else {
        stop(paste0("unknown direction '", direction, "' (use all, up or down)"))
    }
    sig[is.na(sig)] <- FALSE
    results$gene_id[sig]
}

postde_write_empty <- function(path, status, header = c("gene_id", "logFC", "pvalue", "padj", "stat", "mean")) {
    con <- file(path, open = "wt")
    writeLines(paste(header, collapse = "\t"), con)
    close(con)
    status_path <- paste0(path, ".status.json")
    esc <- gsub("([\"\\\\])", "\\\\\\1", status)
    writeLines(paste0('{"status": "', esc, '"}'), status_path)
    invisible(status)
}

postde_read_tsv <- function(path, header = TRUE, ...) {
    if (!file.exists(path)) {
        stop(paste0("file not found: ", path))
    }
    read.delim(path, header = header, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

postde_sanitize_id <- function(id) {
    if (is.null(id) || length(id) != 1 || is.na(id) || !nzchar(id)) {
        stop("postde_sanitize_id: id must be a single non-empty string")
    }
    if (id == "." || id == "..") {
        stop(paste0("postde_sanitize_id: unsafe id '", id, "'"))
    }
    gsub("[^A-Za-z0-9._-]", "_", id)
}

postde_check_sanitized_ids <- function(ids) {
    san <- vapply(ids, postde_sanitize_id, character(1))
    if (anyDuplicated(san)) {
        stop(paste0("postde: duplicate sanitized entry ids: ", paste(unique(san[duplicated(san)]), collapse = ", ")))
    }
    invisible(san)
}

postde_differential_scores <- function(scores, design, contrast) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma is required for differential scores but is not installed")
    }
    scores <- as.matrix(scores)
    design <- as.matrix(design)
    postde_validate_design(design, contrast)
    if (ncol(scores) != nrow(design)) {
        stop(paste0("scores columns (", ncol(scores), ") do not match design rows (", nrow(design), ")"))
    }
    if (!is.null(colnames(scores)) && !is.null(rownames(design)) && !identical(colnames(scores), rownames(design))) {
        stop("scores column names do not match design row names")
    }
    if (nrow(scores) == 0) {
        stop("scores has no features")
    }
    if (any(!is.finite(scores))) {
        stop("scores contains non-finite values")
    }
    fit <- limma::lmFit(scores, design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    tt <- limma::topTable(fit2, number = nrow(scores), sort.by = "none")
    se <- fit2$stdev.unscaled[, 1] * fit2$sigma
    list(
        results = data.frame(
            feature_id = rownames(scores),
            score_difference = tt$logFC,
            pvalue = tt$P.Value,
            padj = tt$adj.P.Val,
            stat = tt$t,
            mean_score = tt$AveExpr,
            SE = se,
            stringsAsFactors = FALSE
        ),
        design = design,
        contrast = contrast
    )
}

postde_gsva_scores <- function(expr, gene_sets, method = "gsva", min_size = 10, max_size = 500) {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("GSVA is required for gsva scores but is not installed")
    }
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop("BiocParallel is required for gsva scores but is not installed")
    }
    expr <- as.matrix(expr)
    if (method == "gsva") {
        param <- GSVA::gsvaParam(exprData = expr, geneSets = gene_sets, minSize = min_size, maxSize = max_size, kcdf = "Gaussian")
    } else if (method == "ssgsea") {
        param <- GSVA::ssgseaParam(exprData = expr, geneSets = gene_sets, minSize = min_size, maxSize = max_size)
    } else {
        stop(paste0("unknown GSVA method '", method, "' (use gsva or ssgsea)"))
    }
    GSVA::gsva(param, BPPARAM = BiocParallel::SerialParam())
}

postde_or <- function(x, default) {
    if (is.null(x)) default else x
}

postde_background_active <- function(config) {
    bg <- config$background
    is.list(bg) && (!is.null(bg$file) || !is.null(bg$expressed))
}

postde_expressed_genes <- function(entry, ex) {
    metric <- ex$metric
    cutoff <- ex$min
    if (is.null(metric) || !metric %in% c("tpm", "cpm", "expression")) {
        stop("background.expressed.metric must be tpm, cpm or expression")
    }
    if (!is.numeric(cutoff) || length(cutoff) != 1 || !is.finite(cutoff)) {
        stop("background.expressed.min must be a single number")
    }
    if (metric == "expression") {
        if (!is.null(ex$min_samples) || !is.null(ex$per_group)) {
            stop("background.expressed: min_samples/per_group are not available for metric expression")
        }
        if (is.null(entry$mean_scale) || !entry$mean_scale %in% c("baseMean", "logCPM")) {
            stop(paste0("background.expressed.metric=expression needs DE results with an average expression column (entry '", entry$id, "' has mean scale '", entry$mean_scale, "')"))
        }
        r <- entry$results
        return(unique(r$gene_id[is.finite(r$mean) & r$mean > cutoff]))
    }
    mat <- entry[[metric]]
    if (is.null(mat)) {
        stop(paste0("background.expressed.metric=", metric, " needs the ", metric, " matrix in the DE bundle (entry '", entry$id, "'), rerun the DE step"))
    }
    mat <- as.matrix(mat)
    pass <- mat > cutoff
    pass[is.na(pass)] <- FALSE
    min_samples <- postde_or(ex$min_samples, 1)
    per_group <- isTRUE(postde_or(ex$per_group, TRUE))
    if (per_group) {
        cond <- as.character(entry$metadata[colnames(mat), "condition"])
        ok <- rep(FALSE, nrow(mat))
        for (cols in split(colnames(mat), cond)) {
            need <- min(min_samples, length(cols))
            ok <- ok | rowSums(pass[, cols, drop = FALSE]) >= need
        }
    } else {
        ok <- rowSums(pass) >= min(min_samples, ncol(mat))
    }
    unique(rownames(mat)[ok])
}

postde_build_background <- function(entry, config, outdir) {
    bg <- config$background
    tested <- postde_ora_universe(entry$results)
    file_genes <- NULL
    expr_genes <- NULL
    if (!is.null(bg$file)) {
        file_genes <- unique(as.character(postde_read_tsv(bg$file, header = FALSE)[[1]]))
    }
    if (!is.null(bg$expressed)) {
        expr_genes <- postde_expressed_genes(entry, bg$expressed)
    }
    sets <- Filter(Negate(is.null), list(file_genes, expr_genes))
    combine <- postde_or(bg$combine, "intersect")
    base <- if (combine == "union") Reduce(union, sets) else Reduce(intersect, sets)
    final <- intersect(base, tested)
    sig <- unique(postde_significant(entry$results, direction = "all", padj = config$padj, lfc = config$lfc))
    below <- setdiff(intersect(sig, tested), final)
    include_sig <- isTRUE(postde_or(bg$include_significant, TRUE))
    if (length(below) > 0) {
        warning(
            paste0(
                "postde[", entry$id, "]: ", length(below), " of ", length(sig),
                " significant genes are not in the expression-derived background (low expression): ",
                paste(head(below, 20), collapse = ", "), if (length(below) > 20) ", ..." else "",
                if (include_sig) ". They are added to the background (background.include_significant=true)." else ". They are dropped from the restricted analyses."
            ),
            call. = FALSE, immediate. = TRUE
        )
        if (include_sig) {
            final <- union(final, below)
        }
    }
    if (length(final) == 0) {
        stop(paste0("postde[", entry$id, "]: background is empty after intersecting with the tested genes"))
    }
    if (length(final) < 100) {
        warning(paste0("postde[", entry$id, "]: restricted background has only ", length(final), " genes"), call. = FALSE, immediate. = TRUE)
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(expr_genes)) {
        writeLines(sort(expr_genes), file.path(outdir, "background_expressed.tsv"))
    }
    writeLines(sort(final), file.path(outdir, "background_final.tsv"))
    summary <- list(
        n_tested = length(tested),
        n_file = if (is.null(file_genes)) NULL else length(file_genes),
        n_expressed = if (is.null(expr_genes)) NULL else length(expr_genes),
        n_final = length(final),
        n_significant = length(sig),
        n_significant_below_cutoff = length(below),
        significant_below_cutoff = below,
        include_significant = include_sig,
        combine = combine,
        criteria = bg
    )
    jsonlite::write_json(summary, file.path(outdir, "background_summary.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
    list(genes = final, summary = summary, dir = outdir)
}

postde_restrict_entry <- function(entry, genes) {
    keep <- rownames(entry$results)[rownames(entry$results) %in% genes]
    if (length(keep) == 0) {
        stop(paste0("postde[", entry$id, "]: no tested genes left in the restricted background"))
    }
    entry$results <- entry$results[keep, , drop = FALSE]
    entry$expression <- entry$expression[keep, , drop = FALSE]
    entry$counts <- entry$counts[rownames(entry$counts) %in% genes, , drop = FALSE]
    entry$cpm <- NULL
    entry$tpm <- NULL
    entry
}
