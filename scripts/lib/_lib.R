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
