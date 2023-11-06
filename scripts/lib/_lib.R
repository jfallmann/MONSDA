## FUNCS
get_gene_name <- function(id, df) {
    name_list <- df$gene[df["type"] == "gene" & df["gene_id"] == id]
    if (length(unique(name_list)) == 1) {
        return(name_list[1])
    } else {
        message(paste("WARNING: ambigous gene id: ", id))
        return(paste(unique(name_list), sep = "|"))
    }
}