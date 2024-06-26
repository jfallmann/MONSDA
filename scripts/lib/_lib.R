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
