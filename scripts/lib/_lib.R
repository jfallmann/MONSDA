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
    if (length(name_list) == 0) {
        message(paste("WARNING: gene id not found in annotation, using id as name: ", id))
        return(id)
    }
    if (length(name_list) == 1) {
        return(name_list)
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
    name_list <- unique(df$gene_name[df["type"] == "exon" & df["gene_id"] == id])
    if (length(name_list) == 0) {
        message(paste("WARNING: gene id not found in annotation, using id as name: ", id))
        return(id)
    }
    if (length(name_list) == 1) {
        return(name_list)
    }
    message(paste("WARNING: ambigous gene id: ", id))
    return(paste(name_list, collapse = "|"))
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
    res$Gene <- unlist(lapply(rownames(res), function(x) {
        get_gene_name(x, gtf_gene)
    }))
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
        Gene = unlist(lapply(gene_ids, function(x) {
            get_gene_name(x, gtf_gene)
        })),
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
