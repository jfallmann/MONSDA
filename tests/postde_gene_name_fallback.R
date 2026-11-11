args <- commandArgs(trailingOnly = TRUE)
src <- readLines(file.path(args[1], "regulatory.R"))
start <- grep("^postde_expression_for_targets", src)
end <- grep("^postde_run_decoupler", src) - 1
eval(parse(text = src[start:end]))

expr <- matrix(c(1:4, 5:8, 9:12, 2:5, 3:6), nrow = 5, byrow = TRUE,
               dimnames = list(paste0("ENSG", 1:5), paste0("s", 1:4)))
gn <- c(ENSG1 = "TP53", ENSG2 = "MYC", ENSG3 = "MYC", ENSG4 = "A|B", ENSG5 = "")
entry <- list(expression = expr, gene_names = gn)

summary <- list()
r <- postde_expression_for_targets(entry, c("TP53", "MYC", "STAT3"))
summary$name_fallback_used <- r$namespace == "gene_name"
summary$name_fallback_rows <- paste(sort(rownames(r$expression)), collapse = ",")
summary$name_fallback_dup_keeps_highest_mean <- all(r$expression["MYC", ] == expr["ENSG3", ])
r2 <- postde_expression_for_targets(entry, c("ENSG1", "TP53"))
summary$id_match_wins <- r2$namespace == "id" && identical(r2$expression, expr)
r3 <- postde_expression_for_targets(entry, c("STAT3", "FOO"))
summary$no_match_keeps_ids <- r3$namespace == "id" && identical(r3$expression, expr)
r4 <- postde_expression_for_targets(list(expression = expr), "TP53")
summary$no_names_keeps_ids <- r4$namespace == "id"
cat(jsonlite::toJSON(summary, auto_unbox = TRUE), "\n")
