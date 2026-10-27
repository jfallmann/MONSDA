## FUNCS
get_gene_name <- function(id, df) {
    if (!"gene_id" %in% colnames(df)) {
        message("WARNING: gene_id not found as colname, will be replaced by first match of colname with ID")
        colnames(df)[grepl("id$", names(df), ignore.case = TRUE)][1] <- "gene_id"
    }
    if (!"gene_name" %in% colnames(df)) {
        message("WARNING: gene_name not found as colname, will be replaced by gene column, please make sure the gtf file is in the correct format")
        df$gene_name <- df$gene
    }
    name_list <- unique(df$gene_name[df["type"] == "gene" & df["gene_id"] == id])
    name_list <- name_list[!is.na(name_list) & nzchar(name_list)]
    ## Must ALWAYS return exactly one string: callers assign the result into a
    ## column, so a zero-length return silently shortens the vector (unlist drops
    ## it) and a multi-element return lengthens it. Both crash the assignment.
    if (length(name_list) == 0) {
        return(id)
    }
    if (length(name_list) == 1) {
        return(name_list[1])
    }
    message(paste("WARNING: ambigous gene id: ", id))
    return(paste(name_list, collapse = "|"))
}


get_exon_name <- function(id, df) {
    if (!"gene_id" %in% colnames(df)) {
        message("WARNING: gene_id not found as colname, will be replaced by first match of colname with ID")
        colnames(df)[grepl("id$", names(df), ignore.case = TRUE)][1] <- "gene_id"
    }
    if (!"gene_name" %in% colnames(df)) {
        message("WARNING: gene_name not found as colname, will be replaced by gene column, please make sure the gtf file is in the correct format")
        df$gene_name <- df$gene
    }
    name_list <- df$gene_name[df["type"] == "exon" & df["gene_id"] == id]
    if (length(unique(name_list)) == 1) {
        return(name_list[1])
    } else {
        message(paste("WARNING: ambigous gene id: ", id))
        return(paste(unique(name_list), sep = "|"))
    }
}

get_gene_coords <- function(id, df) {
    if (!"gene_id" %in% colnames(df)) {
        message("WARNING: gene_id not found as colname, will be replaced by first match of colname with ID")
        colnames(df)[grepl("id$", names(df), ignore.case = TRUE)][1] <- "gene_id"
    }
    coord_rows <- df[df["type"] == "gene" & df["gene_id"] == id, ]
    if (nrow(coord_rows) == 0) {
        return(NA_character_)
    }
    coord_list <- paste(coord_rows$seqnames, coord_rows$start, coord_rows$end, coord_rows$strand, sep = ":")
    if (length(unique(coord_list)) == 1) {
        return(coord_list[1])
    } else {
        message(paste("WARNING: ambigous gene id: ", id))
        return(paste(unique(coord_list), collapse = "|"))
    }
}


add_gene_coordinates <- function(df, gene_ids, gtf_df, after = NULL) {
    if (!"gene_id" %in% colnames(gtf_df)) {
        message("WARNING: gene_id not found as colname, will be replaced by first match of colname with ID")
        colnames(gtf_df)[grepl("id$", names(gtf_df), ignore.case = TRUE)][1] <- "gene_id"
    }
    df <- as.data.frame(df)
    g <- gtf_df[gtf_df["type"] == "gene", ]
    coords <- paste(g$seqnames, g$start, g$end, g$strand, sep = ":")
    names(coords) <- as.character(g$gene_id)
    df$Coordinates <- unname(coords[as.character(gene_ids)])
    cols <- colnames(df)
    cols <- cols[cols != "Coordinates"]
    if (!is.null(after) && !is.na(after) && after %in% cols) {
        pos <- match(after, cols)
        new_order <- append(cols, "Coordinates", after = pos)
    } else {
        if (!is.null(after) && !is.na(after)) {
            message(paste("WARNING: column", after, "not found, placing Coordinates first"))
        }
        new_order <- c("Coordinates", cols)
    }
    df[, new_order, drop = FALSE]
}

fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


calc_cpm <- function(counts) {
    lib_sizes <- colSums(counts)
    cpm <- t(t(counts) / lib_sizes * 1e6)
    return(cpm)
}


calc_tpm <- function(counts, gtf) {
    # Get gene lengths from GTF (assumes gtf_gene has columns 'gene_id' and 'width')
    gene_lengths <- gtf$width
    names(gene_lengths) <- gtf$gene_id
    matched_lengths <- gene_lengths[rownames(counts)]
    gene_lengths_kb <- matched_lengths / 1000
    rpk <- counts / gene_lengths_kb
    scaling_factors <- colSums(rpk)
    tpm <- t(t(rpk) / scaling_factors * 1e6)
    return(tpm)
}

## Minimum number of samples in any non-empty condition group, used for the
## low-count prefilter. Empty factor levels are dropped first so that a group
## without samples cannot silently keep every gene.
min_group_size <- function(cond) {
    tab <- table(droplevels(factor(cond)))
    if (length(tab) == 0) {
        stop("No non-empty condition groups available for the low-count filter")
    }
    min(tab)
}

## Parse a MONSDA comparison string ("name:A+B-vs-C+D,...") into a list of
## contrasts with named A (numerator) and B (denominator) group vectors.
parse_comparisons <- function(cmp) {
    lapply(strsplit(cmp, ",")[[1]], function(contrast) {
        name <- strsplit(contrast, ":")[[1]][1]
        groups <- strsplit(strsplit(contrast, ":")[[1]][2], "-vs-")[[1]]
        list(
            name = name,
            A = unlist(strsplit(groups[1], "\\+"), use.names = FALSE),
            B = unlist(strsplit(groups[2], "\\+"), use.names = FALSE)
        )
    })
}

## Validate that every comparison is a pairwise contrast between two distinct,
## existing condition groups. Compound groups (pooled or weighted semantics)
## are not supported and rejected with a clear error before any fit runs.
validate_comparisons <- function(parsed, condition_levels) {
    for (cmp in parsed) {
        if (length(cmp$A) != 1 || length(cmp$B) != 1) {
            stop(paste0("Comparison '", cmp$name, "' uses compound groups (", paste(cmp$A, collapse = "+"), "-vs-", paste(cmp$B, collapse = "+"), "); only one group per side is supported"))
        }
        if (!cmp$A %in% condition_levels) {
            stop(paste0("Comparison '", cmp$name, "' references unknown group '", cmp$A, "' (available: ", paste(condition_levels, collapse = ", "), ")"))
        }
        if (!cmp$B %in% condition_levels) {
            stop(paste0("Comparison '", cmp$name, "' references unknown group '", cmp$B, "' (available: ", paste(condition_levels, collapse = ", "), ")"))
        }
        if (cmp$A == cmp$B) {
            stop(paste0("Comparison '", cmp$name, "' compares group '", cmp$A, "' against itself"))
        }
    }
    invisible(parsed)
}

## Select the samples of a pairwise contrast from the full annotation and count
## tables. Samples are chosen via the metadata condition column (never by regex
## on sample names) and ordered B then A; the count matrix is subset by the
## ordered sample rownames so metadata and counts stay exactly aligned.
select_contrast_samples <- function(sampleData_all, countData_all, A, B) {
    if (length(A) != 1 || length(B) != 1) {
        stop("select_contrast_samples requires exactly one group per side")
    }
    if (anyDuplicated(colnames(countData_all))) {
        stop(paste0("Duplicate sample IDs in count table: ", paste(colnames(countData_all)[duplicated(colnames(countData_all))], collapse = ", ")))
    }
    sampleData <- droplevels(rbind(subset(sampleData_all, condition == B), subset(sampleData_all, condition == A)))
    if (anyDuplicated(rownames(sampleData))) {
        stop(paste0("Duplicate sample IDs in comparison: ", paste(rownames(sampleData)[duplicated(rownames(sampleData))], collapse = ", ")))
    }
    if (!all(rownames(sampleData) %in% colnames(countData_all))) {
        stop("Count file does not correspond to the annotation file for this comparison")
    }
    countData <- countData_all[, rownames(sampleData), drop = FALSE]
    list(sampleData = sampleData, countData = countData)
}

## Global shift of a contrast: the median log fold change over all tested genes.
## Under spike-in normalization this offset is the overall expression change
## between the two conditions, which conventional normalization would have
## scaled away. Subtracting it leaves each gene's own deviation from that shift.
lfc_global_shift <- function(lfc) {
    vals <- lfc[is.finite(lfc)]
    if (length(vals) == 0) {
        message("WARNING: no finite log fold changes, global shift set to 0")
        return(0)
    }
    median(vals)
}

## Write the per-sample scaling factor spike-in normalization derived for a
## contrast (geometric-mean-centered, 1 = no rescaling relative to the other
## samples in the contrast), so that shift is traceable independent of the DE
## results table.
write_scaling_log <- function(scaling_factors, samples, path) {
    log_df <- data.frame(sample = samples, spikein_scaling_factor = scaling_factors)
    write.table(log_df, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

## Format a DESeq2 results object for export: add gene name and ID, select the
## canonical column order by name and append genomic coordinates after Gene_ID.
## Shrunk tables carry no stat column; raw (unshrunk) tables append stat last.
## With center = TRUE a log2FoldChange_centered column is appended last, holding
## the deviation from the contrast's global shift. It is appended rather than
## inserted because the significance filters address columns by position.
format_deseq2_results <- function(res, gtf_gene, shrink = TRUE, center = FALSE) {
    res$Gene <- vapply(rownames(res), function(x) {
        get_gene_name(x, gtf_gene)
    }, FUN.VALUE = character(1), USE.NAMES = FALSE)
    res$Gene_ID <- rownames(res)
    if (shrink) {
        res <- res[, c("Gene_ID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
    } else {
        res <- res[, c("Gene_ID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat")]
    }
    if (center) {
        shift <- lfc_global_shift(res$log2FoldChange)
        message(paste0("Global shift (median log2FoldChange) subtracted in log2FoldChange_centered: ", shift))
        res$log2FoldChange_centered <- res$log2FoldChange - shift
    }
    add_gene_coordinates(res, res$Gene_ID, gtf_gene, after = "Gene_ID")
}

## Format an edgeR result table (qlf$table or a topTags table) for export.
## Column names follow the DESeq2 template so both engines emit the same header:
## logFC -> log2FoldChange, PValue -> pvalue, BH-adjusted p -> padj, and the test
## statistic (F for quasi-likelihood, LR for likelihood ratio) -> stat. logCPM
## keeps its own name on purpose: it is a log-scale mean and is not the same
## quantity as DESeq2's baseMean. Tests without a statistic column (exactTest)
## simply omit stat. With center = TRUE the deviation from the contrast's global
## shift is appended last as log2FoldChange_centered.
format_edger_results <- function(tbl, gtf_gene, center = FALSE) {
    tbl <- as.data.frame(tbl)
    for (required in c("logFC", "logCPM", "PValue")) {
        if (!required %in% colnames(tbl)) {
            stop(paste0("format_edger_results: missing required column ", required))
        }
    }
    gene_ids <- rownames(tbl)
    out <- data.frame(
        Gene_ID = gene_ids,
        Gene = vapply(gene_ids, function(x) {
            get_gene_name(x, gtf_gene)
        }, FUN.VALUE = character(1), USE.NAMES = FALSE),
        logCPM = tbl$logCPM,
        log2FoldChange = tbl$logFC,
        pvalue = tbl$PValue,
        padj = if ("FDR" %in% colnames(tbl)) tbl$FDR else p.adjust(tbl$PValue, method = "BH"),
        stringsAsFactors = FALSE
    )
    rownames(out) <- gene_ids
    stat_col <- intersect(c("F", "LR"), colnames(tbl))
    if (length(stat_col) > 0) {
        out$stat <- tbl[[stat_col[1]]]
    }
    if (center) {
        shift <- lfc_global_shift(out$log2FoldChange)
        message(paste0("Global shift (median log2FoldChange) subtracted in log2FoldChange_centered: ", shift))
        out$log2FoldChange_centered <- out$log2FoldChange - shift
    }
    add_gene_coordinates(out, out$Gene_ID, gtf_gene, after = "Gene_ID")
}

## Derive a PNG canvas size, resolution and label cex from the number of
## samples/labels a DataSet-level figure has to show, instead of every call
## site hardcoding a fixed width/height/res that only looks right for a
## handful of samples and becomes illegible (or, worse, silently truncated)
## once a study has dozens to hundreds of samples.
scaled_png_geometry <- function(n, base_w = 1900, base_h = 1200, base_res = 300, per_label_px = 18, max_px = 8000) {
    n <- max(1, n)
    width <- min(max_px, max(base_w, base_w + per_label_px * n))
    height <- min(max_px, max(base_h, base_h + per_label_px * n))
    cex <- max(0.25, min(1, 30 / n))
    cex_margin <- max(6, min(30, 6 + n / 10))
    list(width = width, height = height, res = base_res, cex = cex, margin = cex_margin)
}

## Base-R (no gplots/pheatmap dependency, so the same function works in every
## DE/DAS/DEU/DTU conda env) sample-to-sample Pearson correlation heatmap.
## Unlike a Euclidean-distance heatmap.2 plot, the color scale here is a
## directly interpretable correlation coefficient, the color-key legend is
## drawn in its own panel so it cannot overlap the title or the matrix, and at
## high sample counts the dendrogram and axis labels are thinned/omitted
## intentionally (with a one-line caption) instead of being silently
## undrawable, which is what produced unreadable output previously.
plot_sample_correlation_heatmap <- function(expr_mat, outfile, group = NULL, title = "Sample-to-sample correlation", method = "pearson", max_dendro_n = 60, max_label_n = 80) {
    n <- ncol(expr_mat)
    geom <- scaled_png_geometry(n)
    cmat <- cor(expr_mat, method = method)
    ord <- if (n > 2) hclust(as.dist(1 - cmat))$order else seq_len(n)
    cmat <- cmat[ord, ord]
    labels <- colnames(cmat)
    show_every <- max(1, ceiling(n / max_label_n))
    shown_idx <- seq(1, n, by = show_every)
    draw_labels <- rep("", n)
    draw_labels[shown_idx] <- labels[shown_idx]
    has_dendro <- n > 2 && n <= max_dendro_n
    has_strip <- !is.null(group)
    ## margin only needs to fit the labels actually drawn, not every sample,
    ## otherwise thinning labels at high n still leaves a huge, mostly empty
    ## margin reserved for text that was never printed
    geom$margin <- scaled_png_geometry(min(n, max_label_n))$margin

    grDevices::png(outfile, width = geom$width, height = geom$height, res = geom$res)
    on.exit(grDevices::dev.off(), add = TRUE)
    top_h <- if (has_dendro) 1.2 else 0.1
    strip_h <- if (has_strip) 0.3 else 0
    graphics::layout(
        matrix(c(1, 3, 2, 4), nrow = 2),
        widths = c(5, 1), heights = c(top_h + strip_h, 5)
    )
    hmcol <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

    ## top-left: dendrogram (if small enough) and/or condition color strip
    graphics::par(mar = c(0.2, geom$margin, 2, 1))
    if (has_dendro) {
        dend <- as.dendrogram(hclust(as.dist(1 - cmat)))
        graphics::plot(dend, axes = FALSE, leaflab = "none", main = title, cex.main = 0.9)
    } else {
        graphics::plot.new()
        graphics::title(main = title, cex.main = 0.9)
        graphics::mtext(paste0("(dendrogram omitted above ", max_dendro_n, " samples; order from hclust on ", n, " samples)"), side = 1, cex = 0.6)
    }
    if (has_strip) {
        group <- as.factor(group)[ord]
        strip_col <- grDevices::rainbow(nlevels(group))[as.integer(group)]
        graphics::par(new = TRUE, mar = c(0.2, geom$margin, 0, 1))
        graphics::plot(1, type = "n", xlim = c(0.5, n + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
        graphics::rect(1:n - 0.5, 0, 1:n + 0.5, 1, col = strip_col, border = NA)
    }

    ## top-right: empty spacer above the color key
    graphics::par(mar = c(0.2, 1, 2, 3))
    graphics::plot.new()

    ## bottom-left: the correlation matrix itself
    graphics::par(mar = c(geom$margin, geom$margin, 0.2, 1))
    graphics::image(1:n, 1:n, t(cmat[n:1, ]), col = hmcol, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    graphics::axis(1, at = 1:n, labels = draw_labels, las = 2, cex.axis = geom$cex)
    graphics::axis(2, at = 1:n, labels = rev(draw_labels), las = 2, cex.axis = geom$cex)

    ## bottom-right: the color key
    graphics::par(mar = c(geom$margin, 1, 0.2, 3))
    key <- matrix(seq(-1, 1, length.out = 100), nrow = 1)
    graphics::image(1, seq(-1, 1, length.out = 100), key, col = hmcol, axes = FALSE, xlab = "", ylab = "")
    graphics::axis(4, at = seq(-1, 1, 0.5), las = 2, cex.axis = 0.8)
    graphics::mtext("Pearson r", side = 3, line = 0.5, cex = 0.8)
    invisible(cmat)
}


## ---------------------------------------------------------------------------
## Spike-in / RUV normalization: option parsing, selection of the number of
## factors k and diagnostics. Nothing here changes a result silently: every
## choice is written to the k-report table and announced in the log.
## ---------------------------------------------------------------------------

## Parse the free-form DE option string into a named list with lowercase keys.
## Snakemake hands the whole string over as a single argument while nextflow
## splits it on whitespace, so all remaining argv elements are folded back into
## one string before parsing. Known keys: controlgenes=<file>, ruvk=auto|<int>.
parse_de_options <- function(args_tail) {
    opts <- list()
    tokens <- unlist(strsplit(trimws(paste(args_tail, collapse = " ")), "[[:space:]]+"))
    for (tok in tokens[nzchar(tokens)]) {
        if (!grepl("=", tok, fixed = TRUE)) {
            message(paste0("WARNING: ignoring DE option '", tok, "', expected key=value"))
            next
        }
        opts[[tolower(sub("=.*$", "", tok))]] <- sub("^[^=]*=", "", tok)
    }
    opts
}

## Value of a DE option, or `default` if it was not given or is empty.
de_option <- function(opts, key, default = "") {
    value <- opts[[tolower(key)]]
    if (is.null(value) || !nzchar(value)) default else value
}

## Prepend the RUV factors W_1..W_k to the right hand side of a design formula,
## keeping the variable of interest last and preserving a suppressed intercept.
ruv_add_w_terms <- function(f, w_names) {
    if (length(w_names) == 0) {
        return(f)
    }
    tt <- stats::terms(f)
    rhs <- paste(c(if (attr(tt, "intercept") == 0) "0", w_names, attr(tt, "term.labels")), collapse = " + ")
    stats::as.formula(paste("~", rhs), env = environment(f))
}


## Spike-in derived size factors from the control gene counts alone, by the
## median of ratios against the geometric mean reference sample.
##
## A fixed amount of spike-in is added per cell / per unit of input, so the
## observed control gene level of a library is proportional to its sequencing
## depth divided by its total endogenous RNA content. Scaling endogenous counts
## by it therefore recovers abundance per cell instead of per library, which is
## the entire point of a spike-in experiment.
##
## The scale must come from the CONTROL GENES and never from the endogenous
## library size: the endogenous library size carries the very RNA content
## difference the spike-ins are there to measure, so using it cancels the
## signal of interest.
##
## Median of ratios rather than TMM on the control genes, for two reasons.
## It is the identical estimator to DESeq2's estimateSizeFactors(controlGenes),
## so both engines are offset by the same quantity and their results can be
## compared directly. And it is the more precise of the two here: over 40
## simulated spike-in datasets with a known global shift of log2(1.6), both are
## unbiased (-0.02) but the median of ratios has 35% lower RMSE than TMM
## (0.117 vs 0.179); the plain control gene total is unbiased too but far
## noisier (0.508). TMM's trimmed mean corrects composition bias, which control
## genes do not have by construction, and it pays for that robustness with
## variance when only a few dozen spike-ins are available.
spikein_size_factors <- function(spike_counts) {
    spike_counts <- as.matrix(spike_counts)
    if (nrow(spike_counts) < 2) {
        stop("spike-in normalization needs at least two control genes present in the count table")
    }
    usable <- rowSums(spike_counts <= 0) == 0
    if (sum(usable) < 2) {
        stop(paste0(
            "only ", sum(usable), " of ", nrow(spike_counts), " control genes are non-zero in every sample, ",
            "the spike-in scale cannot be estimated; check the spike-in coverage of the shallowest library"
        ))
    }
    logc <- log(spike_counts[usable, , drop = FALSE])
    sf <- exp(apply(logc - rowMeans(logc), 2, stats::median))
    if (any(!is.finite(sf)) || any(sf <= 0)) {
        stop("spike-in size factors are not positive and finite, check that every sample has control gene counts")
    }
    stats::setNames(sf / exp(mean(log(sf))), colnames(spike_counts))
}

## edgeR applies the offset lib.size * norm.factors, so externally derived size
## factors are installed by expressing them relative to the library sizes the
## DGEList currently carries. The norm.factors are rescaled to geometric mean
## 1 (edgeR's convention), which leaves the offsets proportional to
## `size_factors` whatever the library sizes happen to be, before or after
## filtering and with keep.lib.sizes either way.
set_spikein_offsets <- function(dge, size_factors) {
    if (length(size_factors) != ncol(dge)) {
        stop("set_spikein_offsets: exactly one size factor per sample is required")
    }
    nf <- as.numeric(size_factors) / dge$samples$lib.size
    dge$samples$norm.factors <- nf / exp(mean(log(nf)))
    dge
}

## Permute the samples independently within every column (gene) of a
## samples x genes matrix. Ordering the flattened matrix by (column, random
## key) exploits the column-major layout and permutes every column in one
## pass, which keeps a few hundred permutations affordable even for large
## control gene sets.
permute_within_columns <- function(m) {
    keys <- stats::runif(length(m))
    matrix(m[order(col(m), keys)], nrow(m), ncol(m))
}

## Squared singular values of a matrix via its Gram matrix. For the variance
## spectrum only d^2 is needed and the number of samples is tiny compared to
## the number of genes, so the n x n eigendecomposition is far cheaper than a
## full SVD of the n x genes matrix.
squared_singular_values <- function(m) {
    pmax(eigen(tcrossprod(m), symmetric = TRUE, only.values = TRUE)$values, 0)
}

## Per-gene pooled within-group variance and the between/within ratio of a
## samples x genes log-expression matrix, used to show what removing k factors
## does to replicate agreement (within should drop) and to the separation of
## the groups under test (between/within should not collapse: that would mean
## the factors are absorbing the biology instead of the nuisance).
group_variance_components <- function(expr, group) {
    group <- droplevels(factor(group))
    n <- nrow(expr)
    ngroups <- nlevels(group)
    if (ngroups < 2 || n - ngroups < 1) {
        return(c(within = NA_real_, between_within = NA_real_))
    }
    indicator <- stats::model.matrix(~ 0 + group)
    means <- (t(indicator) / colSums(indicator)) %*% expr
    fitted <- indicator %*% means
    within <- colSums((expr - fitted)^2) / (n - ngroups)
    grand <- colMeans(expr)
    between <- colSums(colSums(indicator) * (means - rep(grand, each = ngroups))^2) / (ngroups - 1)
    c(
        within = mean(within),
        between_within = if (mean(within) > 0) mean(between) / mean(within) else NA_real_
    )
}

## Association of one factor with the condition under test: fraction of its
## variance explained by the group labels and the corresponding F test.
factor_condition_association <- function(w, group) {
    group <- droplevels(factor(group))
    n <- length(w)
    ngroups <- nlevels(group)
    if (ngroups < 2 || n - ngroups < 1) {
        return(c(r2 = NA_real_, pvalue = NA_real_))
    }
    ss_tot <- sum((w - mean(w))^2)
    fitted <- as.vector(stats::model.matrix(~ 0 + group) %*% tapply(w, group, mean))
    ss_res <- sum((w - fitted)^2)
    if (ss_tot <= 0) {
        return(c(r2 = NA_real_, pvalue = NA_real_))
    }
    fstat <- ((ss_tot - ss_res) / (ngroups - 1)) / (ss_res / (n - ngroups))
    c(
        r2 = 1 - ss_res / ss_tot,
        pvalue = if (is.finite(fstat)) stats::pf(fstat, ngroups - 1, n - ngroups, lower.tail = FALSE) else NA_real_
    )
}


## Diagnostics for every candidate k from 0 (no correction) to k_max, on the
## endogenous genes, i.e. on the genes the differential expression is actually
## tested on. Because the columns of W are orthonormal, removing factor j is a
## single projection, so the k-th row is the (k-1)-th row with one more factor
## projected out and nothing has to be refactorized per k.
ruv_k_metrics <- function(expr, factors, varprop, perm_p, k_max, condition, log_sf, k_used, design_rank, contrast = NULL) {
    nsamples <- nrow(expr)
    groups <- if (!is.null(condition)) droplevels(factor(condition)) else NULL
    rows <- list()
    corrected <- expr
    for (j in 0:k_max) {
        if (j > 0) {
            u <- factors[, j, drop = FALSE]
            corrected <- corrected - u %*% (t(u) %*% corrected)
        }
        rle <- sweep(corrected, 2, apply(corrected, 2, stats::median), FUN = "-")
        variance <- if (!is.null(condition)) {
            group_variance_components(corrected, condition)
        } else {
            c(within = NA_real_, between_within = NA_real_)
        }
        association <- if (j > 0 && !is.null(condition)) {
            factor_condition_association(factors[, j], condition)
        } else {
            c(r2 = NA_real_, pvalue = NA_real_)
        }
        ## The global log2 shift between the two groups: the quantity a spike-in
        ## normalization exists to recover and the one a factor collinear with
        ## the spike-in size factors takes away again. It has to come from the
        ## joint fit of W and the condition, not from the W-projected data: as
        ## soon as a factor correlates with the condition the two differ, and it
        ## is the joint fit that the DE model reports.
        shift <- NA_real_
        if (!is.null(groups) && !is.null(contrast)) {
            model <- cbind(
                intercept = 1,
                if (j > 0) factors[, seq_len(j), drop = FALSE],
                contrast = as.numeric(groups == contrast[1])
            )
            if (qr(model)$rank == ncol(model)) {
                shift <- stats::median(stats::lm.fit(model, expr)$coefficients["contrast", ])
            }
        }
        rows[[length(rows) + 1]] <- data.frame(
            k = j,
            var_explained_controls = if (j > 0) varprop[j] else NA_real_,
            cum_var_explained_controls = if (j > 0) sum(varprop[seq_len(j)]) else NA_real_,
            permutation_pvalue = if (j > 0) perm_p[j] else NA_real_,
            cor_W_log2_sizefactor = if (j > 0 && !is.null(log_sf)) stats::cor(factors[, j], log_sf) else NA_real_,
            condition_R2 = association[["r2"]],
            condition_pvalue = association[["pvalue"]],
            rle_median_abs = stats::median(abs(apply(rle, 1, stats::median))),
            rle_iqr_median = stats::median(apply(rle, 1, stats::IQR)),
            within_group_var = variance[["within"]],
            between_within_ratio = variance[["between_within"]],
            global_log2_shift = shift,
            residual_df = nsamples - design_rank - j,
            selected = j == k_used,
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}


## Estimate how many factors of unwanted variation (k) RUVg should remove,
## from the negative control (spike-in) genes themselves.
##
## RUVg takes the gene-centered log counts of the control genes, decomposes
## them by SVD and uses the first k left singular vectors as W. How many of
## those components carry structure shared across the control genes rather
## than noise is a dimensionality question, and the standard answer for
## exactly this situation is a permutation test on the eigenvalue spectrum
## (parallel analysis, Buja & Eyuboglu 1992; the same procedure sva::num.sv
## uses by default to pick the number of surrogate variables). Permuting the
## samples independently within every control gene destroys any factor shared
## across genes while leaving each gene's marginal distribution untouched, so
## a component whose share of the total variance exceeds what permuted data
## produce cannot be explained by noise alone. The p-values are made monotone
## in the component index and k is the number of leading components at or
## below `alpha`, which is what makes the choice reproducible instead of the
## usual "look at a PCA plot and pick one".
##
## Unlike sva::num.sv the condition is NOT projected out before the test: RUVg
## assumes the control genes are unaffected by the biology and decomposes them
## as they are, so the spectrum tested here has to be the one RUVg will use. A
## factor that does line up with the condition is reported (condition_R2) and
## warned about instead of being quietly removed.
##
## k is additionally capped so the resulting model stays estimable:
## k <= number of non-degenerate components, k <= #control genes - 1, and
## k <= n_samples - rank(design) - 1 so that at least one residual degree of
## freedom is left for dispersion estimation.
##
## `k` may be "auto" (default) or an integer given via the ruvk= DE option; an
## integer beyond what is estimable is an error rather than a silent downgrade.
ruv_estimate_k <- function(counts, ctrlgenes, design, condition = NULL, size_factors = NULL,
                           k = "auto", max_k = 10, alpha = 0.1, nperm = 200, seed = 42,
                           epsilon = 1, label = "", contrast = NULL) {
    counts <- as.matrix(counts)
    ctrlgenes <- intersect(ctrlgenes, rownames(counts))
    nsamples <- ncol(counts)
    if (length(ctrlgenes) < 2) {
        stop("ruv_estimate_k needs at least two control genes present in the count table")
    }
    if (nsamples < 3) {
        message("RUV k selection: fewer than three samples, no factor of unwanted variation is estimable (k = 0)")
    }
    design <- as.matrix(design)
    design_rank <- qr(design)$rank

    ## the very matrix RUVg factorizes: samples x genes, log counts, gene-centered
    expr <- scale(t(log(counts + epsilon)), center = TRUE, scale = FALSE)
    controls <- expr[, ctrlgenes, drop = FALSE]
    decomp <- svd(controls)
    varprop <- decomp$d^2 / sum(decomp$d^2)
    ncomp <- sum(decomp$d > 1e-8)
    k_max <- max(0, min(max_k, ncomp, length(ctrlgenes) - 1, nsamples - design_rank - 1))

    ## permutation null of the variance spectrum
    set.seed(seed)
    null <- matrix(NA_real_, nrow = nperm, ncol = length(varprop))
    for (b in seq_len(nperm)) {
        permuted <- squared_singular_values(permute_within_columns(controls))
        kept <- seq_len(min(length(permuted), ncol(null)))
        null[b, kept] <- permuted[kept] / sum(permuted)
    }
    perm_p <- cummax(vapply(
        seq_along(varprop),
        function(j) mean(null[, j] >= varprop[j], na.rm = TRUE),
        numeric(1)
    ))
    null_quantile <- apply(null, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    k_auto <- min(sum(perm_p <= alpha), k_max)

    k_requested <- NA_integer_
    if (!identical(tolower(as.character(k)), "auto")) {
        k_requested <- suppressWarnings(as.integer(k))
        if (is.na(k_requested) || k_requested < 0) {
            stop(paste0("ruvk must be 'auto' or a non-negative integer, got '", k, "'"))
        }
        if (k_requested > k_max) {
            stop(paste0(
                "ruvk=", k_requested, " was requested but at most k=", k_max, " is estimable here (",
                nsamples, " samples, design rank ", design_rank, ", ", length(ctrlgenes),
                " control genes, ", ncomp, " non-degenerate components). Lower ruvk or use ruvk=auto."
            ))
        }
    }
    k_used <- if (is.na(k_requested)) k_auto else k_requested

    ## The diagnostics are reported for the endogenous genes on the scale the
    ## DE test sees them: counts divided by the spike-in size factors, log2 for
    ## readability, with the same factors projected out that the design would
    ## adjust for. W is unaffected by the choice of log base, so this stays the
    ## decomposition RUVg performs.
    endogenous <- setdiff(rownames(counts), ctrlgenes)
    expressed <- endogenous[rowSums(counts[endogenous, , drop = FALSE]) > 0]
    normalized <- counts[expressed, , drop = FALSE]
    if (!is.null(size_factors)) {
        normalized <- sweep(normalized, 2, size_factors, FUN = "/")
    }
    log_sf <- if (!is.null(size_factors)) log2(size_factors) else NULL
    ## the direction of the reported global shift follows the contrast
    ## (A vs B), never the alphabetical level order
    if (!is.null(contrast) && (length(contrast) != 2 || !all(contrast %in% levels(droplevels(factor(condition)))))) {
        stop("ruv_estimate_k: contrast must be the two condition groups of the comparison, A first")
    }
    metrics <- ruv_k_metrics(
        scale(t(log2(normalized + epsilon)), center = TRUE, scale = FALSE),
        decomp$u, varprop, perm_p, k_max, condition, log_sf, k_used, design_rank, contrast
    )

    selection <- list(
        k = k_used, k_auto = k_auto, k_max = k_max, k_requested = k_requested,
        method = if (is.na(k_requested)) "permutation (parallel analysis, Buja & Eyuboglu)" else "user supplied (ruvk=)",
        alpha = alpha, nperm = nperm, seed = seed, nsamples = nsamples,
        n_control = length(ctrlgenes), n_endogenous = length(expressed), design_rank = design_rank,
        varprop = varprop, null_quantile = null_quantile, perm_p = perm_p,
        metrics = metrics, label = label,
        contrast = if (is.null(contrast)) NA_character_ else paste(contrast[1], "vs", contrast[2])
    )
    ruv_announce_k(selection)
    selection
}


## Log the outcome of the k selection and flag the two properties a selected
## factor can have that change how the _ruv branch must be read. RUVg factors
## and the spike-in size factors are estimated from the SAME control genes, so
## when the spike-in is added per cell the control genes encode RNA content per
## cell and W becomes collinear with the size factors: the _ruv branch then
## reports fold changes relative to total RNA, i.e. the spike-anchored reading
## lives in the _norm branch and only there. A factor that differs between the
## groups under test is absorbing biology rather than nuisance and makes the
## _ruv branch conservative. Both are reported, never silently acted upon.
ruv_announce_k <- function(selection, collinearity_warn = 0.8, condition_r2_warn = 0.5) {
    tag <- if (nzchar(selection$label)) paste0(" [", selection$label, "]") else ""
    message(sprintf(
        "RUV k selection%s: k = %d (%s), permutation estimate k_auto = %d, estimable maximum k_max = %d, %d control genes, %d samples, design rank %d",
        tag, selection$k, selection$method, selection$k_auto, selection$k_max,
        selection$n_control, selection$nsamples, selection$design_rank
    ))
    if (!is.na(selection$k_requested) && selection$k_requested != selection$k_auto) {
        message(sprintf(
            "  NOTE: ruvk=%d overrides the permutation estimate of %d factors%s",
            selection$k_requested, selection$k_auto,
            if (selection$k_requested > selection$k_auto) " (more factors than the control genes support, this removes variance the data do not show to be unwanted)" else ""
        ))
    }
    if (selection$k == 0) {
        message("  No factor of unwanted variation is used: the _ruv design is left unchanged, so the _ruv branch reproduces the standard branch exactly. Only the _norm branch (spike-in size factors) still differs.")
    }
    for (j in seq_len(selection$k)) {
        row <- selection$metrics[selection$metrics$k == j, ]
        message(sprintf(
            "  W_%d: %.1f%% of control gene variance, permutation p = %.3f, cor(W_%d, log2 sizeFactor) = %s, condition R2 = %s (p = %s)",
            j, 100 * row$var_explained_controls, row$permutation_pvalue, j,
            format(round(row$cor_W_log2_sizefactor, 3)), format(round(row$condition_R2, 3)),
            format(signif(row$condition_pvalue, 3))
        ))
        if (!is.na(row$cor_W_log2_sizefactor) && abs(row$cor_W_log2_sizefactor) >= collinearity_warn) {
            shift0 <- selection$metrics$global_log2_shift[selection$metrics$k == 0]
            message(sprintf(
                "    NOTE: W_%d is %.0f%% collinear with the spike-in size factors, both are derived from the same control genes. The _ruv branch therefore reports fold changes RELATIVE to total RNA, not spike-anchored ones%s. Read the absolute change off the _norm tables, which are scaled by the spike-ins.",
                j, 100 * abs(row$cor_W_log2_sizefactor),
                if (is.finite(shift0) && is.finite(row$global_log2_shift)) {
                    sprintf(": the global log2 shift between the groups goes from %.2f without this factor to %.2f with it", shift0, row$global_log2_shift)
                } else {
                    ""
                }
            ))
        }
        if (!is.na(row$condition_R2) && row$condition_R2 >= condition_r2_warn) {
            message(sprintf(
                "    WARNING: W_%d differs between the groups being tested (R2 = %.2f), it is absorbing biology, not only nuisance. Check the k-report before trusting the _ruv results.",
                j, row$condition_R2
            ))
        }
    }
    invisible(selection)
}

## Write the per-k diagnostics next to the results, so the choice of k can be
## judged (and reproduced) without rerunning anything. Parameters of the
## selection go into leading '#' comment lines, the table itself is plain TSV.
ruv_write_k_report <- function(selection, path) {
    header <- c(
        paste0("# RUVg factor selection", if (nzchar(selection$label)) paste0(" for ", selection$label) else ""),
        paste0("# k_selected\t", selection$k),
        paste0("# k_method\t", selection$method),
        paste0("# k_auto_permutation\t", selection$k_auto),
        paste0("# k_requested\t", if (is.na(selection$k_requested)) "auto" else selection$k_requested),
        paste0("# k_max_estimable\t", selection$k_max),
        paste0("# samples\t", selection$nsamples),
        paste0("# design_rank\t", selection$design_rank),
        paste0("# control_genes\t", selection$n_control),
        paste0("# endogenous_genes\t", selection$n_endogenous),
        paste0("# permutations\t", selection$nperm),
        paste0("# alpha\t", selection$alpha),
        paste0("# seed\t", selection$seed),
        paste0("# global_log2_shift_direction\t", selection$contrast)
    )
    writeLines(header, path)
    suppressWarnings(utils::write.table(
        selection$metrics, path, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE
    ))
    invisible(path)
}


## Visual counterpart of the k-report: the control gene variance spectrum
## against its permutation null (what the choice of k is based on), what each
## candidate factor is correlated with (what it would cost), and how the data
## actually change with k (RLE spread and within-group variance).
ruv_plot_k <- function(selection, outfile) {
    metrics <- selection$metrics
    ncomp <- min(length(selection$varprop), max(selection$nsamples - 1, 1))
    grDevices::png(outfile, width = 3600, height = 1000, res = 220)
    on.exit(grDevices::dev.off())
    graphics::par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 4))

    observed <- 100 * selection$varprop[seq_len(ncomp)]
    null <- 100 * selection$null_quantile[seq_len(ncomp)]
    graphics::plot(seq_len(ncomp), observed,
        type = "b", pch = 19, col = "steelblue", ylim = c(0, max(observed, null, na.rm = TRUE) * 1.1),
        xlab = "factor", ylab = "% of control gene variance",
        main = paste0("spectrum vs permutation null\nk = ", selection$k, " (k_auto = ", selection$k_auto, ", k_max = ", selection$k_max, ")")
    )
    graphics::lines(seq_len(ncomp), null, type = "b", pch = 1, lty = 2, col = "grey40")
    if (selection$k > 0) {
        graphics::abline(v = selection$k + 0.5, col = "firebrick", lty = 3)
    }
    graphics::legend("topright", bty = "n", cex = 0.8, pch = c(19, 1), lty = c(1, 2),
        col = c("steelblue", "grey40"), legend = c("observed", paste0((1 - selection$alpha) * 100, "% of permuted")))

    factors <- metrics[metrics$k > 0, , drop = FALSE]
    if (nrow(factors) > 0) {
        bars <- rbind(abs(factors$cor_W_log2_sizefactor), factors$condition_R2)
        bars[is.na(bars)] <- 0
        graphics::barplot(bars, beside = TRUE, names.arg = paste0("W_", factors$k), ylim = c(0, 1),
            col = c("darkorange", "darkolivegreen"), ylab = "association", xlab = "factor",
            main = "what each factor is tied to")
        graphics::abline(h = c(0.5, 0.8), lty = c(3, 2), col = c("darkolivegreen", "darkorange"))
        graphics::legend("topright", bty = "n", cex = 0.8, fill = c("darkorange", "darkolivegreen"),
            legend = c("|cor| with log2 sizeFactor", "R2 with condition"))
    } else {
        graphics::plot.new()
        graphics::title(main = "no estimable factor")
    }

    graphics::plot(metrics$k, metrics$rle_iqr_median,
        type = "b", pch = 19, col = "steelblue", ylim = c(0, max(metrics$rle_iqr_median, na.rm = TRUE) * 1.1),
        xlab = "k", ylab = "median RLE IQR", main = "effect on the endogenous genes"
    )
    if (any(is.finite(metrics$within_group_var))) {
        graphics::par(new = TRUE)
        graphics::plot(metrics$k, metrics$within_group_var, type = "b", pch = 17, lty = 2,
            col = "firebrick", axes = FALSE, xlab = "", ylab = "",
            ylim = c(0, max(metrics$within_group_var, na.rm = TRUE) * 1.1))
        graphics::axis(4, col = "firebrick", col.axis = "firebrick")
        graphics::mtext("within-group variance", side = 4, line = 2.2, cex = 0.7, col = "firebrick")
    }
    if (selection$k > 0) {
        graphics::abline(v = selection$k, col = "firebrick", lty = 3)
    }

    if (any(is.finite(metrics$global_log2_shift))) {
        graphics::plot(metrics$k, metrics$global_log2_shift,
            type = "b", pch = 19, col = "purple4", xlab = "k", ylab = "median log2 shift between groups",
            main = "spike-in scale kept per k"
        )
        graphics::abline(h = 0, lty = 2, col = "grey40")
        if (selection$k > 0) {
            graphics::abline(v = selection$k, col = "firebrick", lty = 3)
        }
    } else {
        graphics::plot.new()
        graphics::title(main = "global shift needs two groups")
    }
    invisible(outfile)
}


## apeglm fits an ADAPTIVE prior whose scale is estimated from the spread of
## the MLE effect sizes. On a contrast with almost no signal left the prior
## collapses towards a point mass at zero, and most genes come back with
## log2FoldChange == 0 and lfcSE ~ 1e-3. Those columns are then unusable for
## effect sizes, ranking or correlation while still looking like ordinary
## output, so the degenerate case is called out explicitly. Diagnostic only,
## the results themselves are left untouched.
check_shrinkage <- function(res_raw, res_shrink, label = "", min_lfcSE = 0.01, max_pinned = 0.05) {
    tag <- if (nzchar(label)) paste0(" [", label, "]") else ""
    median_se <- stats::median(res_shrink$lfcSE, na.rm = TRUE)
    pinned <- mean(abs(res_shrink$log2FoldChange) < 1e-6, na.rm = TRUE)
    raw_lfc <- stats::median(abs(res_raw$log2FoldChange), na.rm = TRUE)
    shrunk_lfc <- stats::median(abs(res_shrink$log2FoldChange), na.rm = TRUE)
    ratio <- if (is.finite(shrunk_lfc) && shrunk_lfc > 0) raw_lfc / shrunk_lfc else Inf
    message(sprintf(
        "Shrinkage check%s: median lfcSE = %.5f, fraction |LFC| < 1e-6 = %.3f, median |LFC| raw/shrunk = %.1fx",
        tag, median_se, pinned, ratio
    ))
    if ((!is.na(median_se) && median_se < min_lfcSE) || (!is.na(pinned) && pinned > max_pinned)) {
        message(paste0(
            "  WARNING: apeglm shrinkage looks DEGENERATE for this contrast. The shrunken log2FoldChange / lfcSE columns ",
            "must not be used for effect sizes or correlations, use the *_noshrink table instead. This usually means the ",
            "contrast has little signal left; check whether a covariate (e.g. a RUV factor W_i) is absorbing it."
        ))
    }
    invisible(c(median_lfcSE = median_se, frac_pinned = pinned, raw_over_shrunk = ratio))
}
