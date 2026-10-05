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
    name_list <- df$gene_name[df["type"] == "gene" & df["gene_id"] == id]
    if (length(unique(name_list)) == 1) {
        return(name_list[1])
    } else {
        message(paste("WARNING: ambigous gene id: ", id))
        return(paste(unique(name_list), sep = "|"))
    }
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

## Format a DESeq2 results object for export: add gene name and ID, select the
## canonical column order by name and append genomic coordinates after Gene_ID.
## Shrunk tables carry no stat column; raw (unshrunk) tables append stat last.
format_deseq2_results <- function(res, gtf_gene, shrink = TRUE) {
    res$Gene <- unlist(lapply(rownames(res), function(x) {
        get_gene_name(x, gtf_gene)
    }))
    res$Gene_ID <- rownames(res)
    if (shrink) {
        res <- res[, c("Gene_ID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
    } else {
        res <- res[, c("Gene_ID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat")]
    }
    add_gene_coordinates(res, res$Gene_ID, gtf_gene, after = "Gene_ID")
}
