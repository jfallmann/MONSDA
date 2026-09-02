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
