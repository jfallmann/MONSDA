#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else tempfile("postde_fixture_")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])))
postde_dir <- file.path(script_dir, "..", "scripts", "Analysis", "PostDE")
source(file.path(postde_dir, "common.R"))
source(file.path(postde_dir, "export.R"))
source(file.path(postde_dir, "enrichment.R"))
source(file.path(postde_dir, "regulatory.R"))

suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(GSVA)
    library(BiocParallel)
    library(decoupleR)
    library(variancePartition)
    library(clusterProfiler)
    library(GO.db)
    library(jsonlite)
})

set.seed(42)
summary <- list()

## synthetic cohort: 300 genes, 8 samples, 2 conditions
ngenes <- 300
samples <- paste0("s", 1:8)
cond <- factor(rep(c("B", "A"), each = 4), levels = c("B", "A"))
metadata <- data.frame(condition = cond, row.names = samples)
base <- exp(rnorm(ngenes, mean = 6, sd = 1.5))
effect <- rnorm(ngenes, mean = 0, sd = 1)
counts <- matrix(rnbinom(ngenes * 8, mu = base, size = 5), nrow = ngenes, dimnames = list(paste0("g", 1:ngenes), samples))
counts[, cond == "A"] <- matrix(rnbinom(ngenes * 4, mu = base * exp(effect), size = 5), nrow = ngenes)
counts <- round(counts)
counts[counts < 0] <- 0

## real DESeq2 fit
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "B")
dds <- DESeq(dds, betaPrior = FALSE)
res <- results(dds, contrast = c("condition", "A", "B"))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
res_df <- postde_results_deseq2(res)
summary$deseq2_rownames_aligned <- identical(rownames(res_df), rownames(res))
summary$deseq2_gene_id_matches <- all(rownames(res_df) == res_df$gene_id)

## real edgeR fit
dge <- DGEList(counts = counts, group = cond)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~condition, data = metadata)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, contrast = c(0, 1))
qlf_df <- postde_results_edger(qlf)
summary$edger_rownames_aligned <- identical(rownames(qlf_df), rownames(qlf$table))
summary$edger_gene_id_matches <- all(rownames(qlf_df) == qlf_df$gene_id)
summary$edger_stat_finite <- all(is.finite(qlf_df$stat[is.finite(qlf$table$F)]))
summary$edger_stat_na_noF <- all(is.na(qlf_df$stat[!is.finite(qlf$table$F)]))

## capture bundles (deseq2 and edger separately)
postde_bundles <- list()
postde_capture(engine = "deseq2", id = "A_vs_B", A = "A", B = "B", normalized = FALSE, metadata = metadata, counts = counts, expression = assay(vsd), formula = ~condition, results = res_df, mean_scale = "baseMean")
postde_write(outdir, "test", "deseq2")
deseq2_bundle <- readRDS(file.path(outdir, "DE_deseq2_test_postde.rds"))
deseq2_entry <- deseq2_bundle$contrasts[["A_vs_B"]]

postde_bundles <- list()
postde_capture(engine = "edger", id = "A_vs_B", A = "A", B = "B", normalized = FALSE, metadata = metadata, counts = counts, expression = cpm(dge, log = TRUE), formula = ~condition, results = qlf_df, mean_scale = "logCPM")
postde_write(outdir, "test", "edger")
edger_bundle <- readRDS(file.path(outdir, "DE_edger_test_postde.rds"))
edger_entry <- edger_bundle$contrasts[["A_vs_B"]]

## capture must reject misaligned rownames (known defect guard)
bad_res <- res_df
rownames(bad_res) <- rev(rownames(bad_res))
summary$capture_rejects_misaligned <- inherits(tryCatch(
    postde_capture(engine = "deseq2", id = "bad", A = "A", B = "B", normalized = FALSE, metadata = metadata, counts = counts, expression = assay(vsd), formula = ~condition, results = bad_res, mean_scale = "baseMean"),
    error = function(e) e
), "error")

## counts superset allowed: counts with extra genes
extra_counts <- rbind(counts, matrix(1, nrow = 5, ncol = 8, dimnames = list(paste0("extra", 1:5), samples)))
postde_bundles <- list()
postde_capture(engine = "deseq2", id = "superset", A = "A", B = "B", normalized = FALSE, metadata = metadata, counts = extra_counts, expression = assay(vsd), formula = ~condition, results = res_df, mean_scale = "baseMean")
summary$capture_counts_superset_ok <- length(postde_bundles) == 1

## GSVA scores vs direct call
gene_sets <- list(
    set1 = paste0("g", 1:20),
    set2 = paste0("g", 50:80),
    set3 = paste0("g", 100:130),
    set4 = paste0("g", 200:230)
)
expr <- assay(vsd)
scores <- postde_gsva_scores(expr, gene_sets, method = "gsva", min_size = 5, max_size = 100)
param <- GSVA::gsvaParam(exprData = expr, geneSets = gene_sets, minSize = 5, maxSize = 100, kcdf = "Gaussian")
scores_direct <- GSVA::gsva(param, BPPARAM = BiocParallel::SerialParam())
summary$gsva_max_diff <- max(abs(scores - scores_direct))

## differential scores vs direct limma
contrast <- deseq2_entry$contrast
diff <- postde_differential_scores(scores, deseq2_entry$design, contrast)
fit <- lmFit(scores, deseq2_entry$design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, number = nrow(scores), sort.by = "none")
summary$diff_score_difference_max_diff <- max(abs(diff$results$score_difference - tt$logFC))
summary$diff_stat_max_diff <- max(abs(diff$results$stat - tt$t))
summary$diff_columns <- paste(colnames(diff$results), collapse = ",")
summary$diff_feature_names_preserved <- identical(diff$results$feature_id, rownames(scores))

## decoupler ULM/MLM vs direct calls
net <- data.frame(
    source = c("TF1", "TF1", "TF1", "TF2", "TF2", "TF2", "TF3", "TF3", "TF3"),
    target = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"),
    mor = c(1, 1, -1, 1, -1, 1, -1, 1, 1)
)
net_path <- file.path(outdir, "test_network.tsv")
write.table(net, net_path, sep = "\t", row.names = FALSE, quote = FALSE)
dc_config <- list(enabled = TRUE, network = net_path, min_size = 1, methods = list("ulm", "mlm"))
dc_out <- postde_run_decoupler(deseq2_entry, list(min_size = 1, decoupler = dc_config), outdir)
acts_ulm <- decoupleR::run_ulm(expr, net, minsize = 1)
acts_mlm <- decoupleR::run_mlm(expr, net, minsize = 1)
ulm_mat <- postde_activity_matrix(acts_ulm, rownames(deseq2_entry$design))
mlm_mat <- postde_activity_matrix(acts_mlm, rownames(deseq2_entry$design))
summary$ulm_max_diff <- max(abs(read.delim(dc_out$ulm$activities, row.names = 1, check.names = FALSE) - ulm_mat))
summary$mlm_max_diff <- max(abs(read.delim(dc_out$mlm$activities, row.names = 1, check.names = FALSE) - mlm_mat))
summary$decoupler_network_meta_exists <- file.exists(file.path(outdir, "decoupler_network_meta.rds"))
summary$decoupler_filtered_tsv_exists <- file.exists(file.path(outdir, "decoupler_network_filtered.tsv"))

## decoupler gating: no allow_network -> error
dc_net_config <- list(enabled = TRUE, resource = "collectri", organism = "human", min_size = 1, methods = list("ulm"))
summary$decoupler_gated_without_network <- inherits(tryCatch(
    postde_run_decoupler(deseq2_entry, list(min_size = 1, decoupler = dc_net_config), outdir),
    error = function(e) e
), "error")

## decoupler MLM rank deficiency -> error (collinear sources)
collinear_net <- data.frame(
    source = c("TF1", "TF1", "TF2", "TF2"),
    target = c("g1", "g2", "g1", "g2"),
    mor = c(1, 1, 2, 2)
)
collinear_path <- file.path(outdir, "collinear_network.tsv")
write.table(collinear_net, collinear_path, sep = "\t", row.names = FALSE, quote = FALSE)
dc_collinear <- list(enabled = TRUE, network = collinear_path, min_size = 1, methods = list("mlm"))
summary$mlm_rank_deficiency_rejected <- inherits(tryCatch(
    postde_run_decoupler(deseq2_entry, list(min_size = 1, decoupler = dc_collinear), outdir),
    error = function(e) e
), "error")

## GO expansion vs clusterProfiler::buildGOmap
t2g <- data.frame(
    term = c("GO:0006915", "GO:0006915", "GO:0005737", "GO:0003677"),
    gene = c("g1", "g2", "g3", "g4")
)
gomap <- postde_build_gomap(t2g)
gomap_direct <- clusterProfiler::buildGOmap(t2g)
summary$gomap_identical <- identical(gomap, gomap_direct)
summary$gomap_expanded_rows <- nrow(gomap)
summary$gomap_unknown_term_rejected <- inherits(tryCatch(
    postde_build_gomap(data.frame(term = "GO:FAKE0001", gene = "g1")),
    error = function(e) e
), "error")

## term2gene/term2name strict validation
t2g_ora <- data.frame(
    term = c(
        rep("GO:0006915", 25), rep("GO:0005737", 25),
        rep("GO:0003677", 25), rep("GO:0005575", 25)
    ),
    gene = c(paste0("g", 1:25), paste0("g", 26:50), paste0("g", 51:75), paste0("g", 76:100))
)
t2g_path <- file.path(outdir, "term2gene.tsv")
write.table(t2g_ora, t2g_path, sep = "\t", row.names = FALSE, quote = FALSE)
t2n_path <- file.path(outdir, "term2name.tsv")
write.table(data.frame(term = c("GO:0006915", "GO:0005737", "GO:0003677", "GO:0005575"), name = c("apoptotic process", "cytoplasm", "DNA binding", "cellular_component")), t2n_path, sep = "\t", row.names = FALSE, quote = FALSE)
summary$term2gene_strict <- inherits(tryCatch(postde_read_term2gene(t2g_path), error = function(e) e), "data.frame")
summary$term2name_strict <- inherits(tryCatch(postde_read_term2name(t2n_path), error = function(e) e), "data.frame")
summary$term2gene_extra_col_rejected <- inherits(tryCatch(
    postde_read_term2gene(file.path(outdir, "term2gene_extra.tsv")),
    error = function(e) e
), "error")
write.table(cbind(t2g, extra = 1), file.path(outdir, "term2gene_extra.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## ORA all/up/down separate files via postde_run_enrichment
cp_config <- list(
    enabled = TRUE,
    term2gene = t2g_path,
    term2name = t2n_path,
    go_expand = FALSE
)
enr_out <- postde_run_enrichment(deseq2_entry, list(padj = 0.99, lfc = 0, min_size = 1, max_size = 500, seed = 1, clusterprofiler = cp_config), outdir)
summary$ora_separate_files <- all(file.exists(file.path(outdir, c("enricher_all_result.tsv", "enricher_up_result.tsv", "enricher_down_result.tsv"))))
summary$enrichment_provenance_exists <- file.exists(file.path(outdir, "enrichment_provenance.rds"))

## GSEA rank unavailable vs no hits
no_stat <- deseq2_entry$results
no_stat$stat <- NULL
gsea_no_stat <- postde_run_gsea(list(results = no_stat), t2g, NULL, 1, 500, 1, outdir)
summary$gsea_rank_unavailable_status <- gsea_no_stat$status
gsea_ok <- postde_run_gsea(deseq2_entry, t2g, NULL, 1, 500, 1, outdir)
summary$gsea_ranked_saved <- file.exists(file.path(outdir, "gsea_ranked.rds"))
summary$gsea_status <- gsea_ok$status

## dream: repeated subjects fixture (8 subjects, 16 samples, 300 genes)
subjects <- paste0("subj", 1:8)
dream_samples <- paste0("d", 1:16)
dream_cond <- factor(rep(rep(c("B", "A"), each = 1), 8), levels = c("B", "A"))
dream_subject <- factor(rep(subjects, each = 2))
dream_meta <- data.frame(condition = dream_cond, subject = dream_subject, row.names = dream_samples)
dream_base <- exp(rnorm(ngenes, mean = 6, sd = 1.5))
dream_effect <- rnorm(ngenes, mean = 0, sd = 1)
dream_counts <- matrix(rnbinom(ngenes * 16, mu = dream_base, size = 5), nrow = ngenes, dimnames = list(paste0("g", 1:ngenes), dream_samples))
dream_counts[, dream_cond == "A"] <- matrix(rnbinom(ngenes * 8, mu = dream_base * exp(dream_effect), size = 5), nrow = ngenes)
dream_counts <- round(dream_counts)
dream_meta_path <- file.path(outdir, "dream_metadata.tsv")
write.table(data.frame(sample = dream_samples, subject = dream_subject), dream_meta_path, sep = "\t", row.names = FALSE, quote = FALSE)
dream_entry <- list(
    id = "dream_test",
    A = "A",
    B = "B",
    normalized = FALSE,
    metadata = dream_meta,
    counts = dream_counts,
    expression = log2(dream_counts + 1),
    design = model.matrix(~condition, data = dream_meta),
    contrast = postde_contrast(~condition, dream_meta, "A", "B"),
    results = NULL,
    mean_scale = "baseMean"
)
dm_config <- list(enabled = TRUE, metadata = dream_meta_path, formula = "~ condition + (1|subject)")
dream_out <- postde_run_dream(dream_entry, list(dream = dm_config), outdir)
dream_res <- read.delim(file.path(outdir, "dream_result.tsv"), check.names = FALSE)
fixed_design <- model.matrix(~condition, data = dream_meta)
dge_d <- DGEList(counts = dream_counts)
keep_d <- filterByExpr(dge_d, design = fixed_design)
dge_d <- dge_d[keep_d, , keep.lib.sizes = FALSE]
dge_d <- calcNormFactors(dge_d, method = "TMM")
vobj_d <- voom(dge_d, fixed_design)
fit_d <- lmFit(vobj_d, fixed_design)
fit_d <- eBayes(fit_d)
tt_d <- topTable(fit_d, number = Inf, sort.by = "none")
common_genes <- intersect(dream_res$gene_id, rownames(tt_d))
summary$dream_max_logfc_diff <- max(abs(dream_res$logFC[match(common_genes, dream_res$gene_id)] - tt_d$logFC[match(common_genes, rownames(tt_d))]))
summary$dream_varpart_status <- dream_out$varpart$status
summary$dream_varpart_status_file <- file.exists(file.path(outdir, "dream_varpart_status.txt"))
vp_ok <- tryCatch(
    variancePartition::fitExtractVarPartModel(vobj_d, ~ (1 | condition) + (1 | subject), dream_meta),
    error = function(e) e
)
summary$dream_varpart_ok_path <- !inherits(vp_ok, "error")
summary$dream_skipped_normalized <- postde_run_dream(list(normalized = TRUE), list(dream = dm_config), outdir)$status

## dream safe formula/contrast validation
summary$dream_safe_formula_ok <- inherits(tryCatch(postde_safe_formula(~condition + (1|subject)), error = function(e) e), "logical")
summary$dream_unsafe_formula_rejected <- inherits(tryCatch(
    postde_safe_formula(~condition + system("echo hi")),
    error = function(e) e
), "error")
summary$dream_safe_contrast_ok <- inherits(tryCatch(postde_safe_contrast_expr("conditionA - conditionB"), error = function(e) e), "logical")
summary$dream_unsafe_contrast_rejected <- inherits(tryCatch(
    postde_safe_contrast_expr("conditionA; system('echo hi')"),
    error = function(e) e
), "error")

## empty TSV header-first with status sidecar
empty_path <- file.path(outdir, "empty_result.tsv")
postde_write_empty(empty_path, "no significant genes")
empty_lines <- readLines(empty_path)
summary$empty_tsv_header_first <- identical(empty_lines[1], "gene_id\tlogFC\tpvalue\tpadj\tstat\tmean")
summary$empty_tsv_status_sidecar <- file.exists(paste0(empty_path, ".status.json"))

## significant thresholds strict padj< and |lfc|>=
sig_res <- data.frame(
    gene_id = c("a", "b", "c", "d"),
    logFC = c(1, 1, 2, 2),
    pvalue = c(0.01, 0.01, 0.01, 0.01),
    padj = c(0.05, 0.049, 0.05, 0.049),
    stringsAsFactors = FALSE
)
sig_all <- postde_significant(sig_res, direction = "all", padj = 0.05, lfc = 1)
summary$sig_threshold_strict <- identical(sig_all, c("b", "d"))
sig_na <- sig_res
sig_na$logFC[1] <- NA
summary$sig_never_select_na_effect <- !any(is.na(postde_significant(sig_na, direction = "all", padj = 0.05, lfc = 1)))

## differential scores singleton handling: 0 residual df -> clean error
singleton_meta <- data.frame(condition = factor(c("B", "A")), row.names = c("x1", "x2"))
singleton_design <- model.matrix(~condition, data = singleton_meta)
singleton_contrast <- postde_contrast(~condition, singleton_meta, "A", "B")
singleton_scores <- matrix(rnorm(4), nrow = 2, dimnames = list(c("f1", "f2"), c("x1", "x2")))
summary$singleton_no_residual_df_rejected <- inherits(tryCatch(
    postde_differential_scores(singleton_scores, singleton_design, singleton_contrast),
    error = function(e) e
), "error")

## design validation: named contrast aligned by name, zero intercept rejected
named_contrast <- setNames(c(0, -1), c("(Intercept)", "conditionA"))
summary$design_named_contrast_aligned <- inherits(tryCatch(
    postde_validate_design(design, named_contrast),
    error = function(e) e
), "logical")
summary$design_zero_intercept_rejected <- inherits(tryCatch(
    postde_validate_design(design, c(1, 0)),
    error = function(e) e
), "error")

## safe ids
summary$sanitize_dot_rejected <- inherits(tryCatch(postde_sanitize_id(".."), error = function(e) e), "error")
summary$sanitize_collision_detected <- inherits(tryCatch(
    postde_check_sanitized_ids(c("a/b", "a_b")),
    error = function(e) e
), "error")

writeLines(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(outdir, "summary.json"))
cat("POSTDE_FIXTURE_OK", outdir, "\n")
