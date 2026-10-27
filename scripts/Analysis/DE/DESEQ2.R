# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

suppressPackageStartupMessages({
    require(utils)
    require(BiocParallel)
    require(DESeq2)
    require(ggrepel)
    require(ggplot2)
    require(rtracklayer)
    require(RUVSeq)
    require(dplyr)
    require(GenomeInfoDb)
    require(apeglm)
    require(EnhancedVolcano)
})

options(echo = TRUE)

## ARGS
args <- commandArgs(trailingOnly = TRUE)
argsLen <- length(args)
anname <- args[1]
countfile <- args[2]
gtf <- args[3]
outdir <- args[4]
cmp <- args[5]
combi <- args[6]
availablecores <- as.integer(args[7])
## everything beyond the fixed arguments is the free-form DE option string
## ("controlgenes=<file> ruvk=auto|<int>"), which snakemake passes as one
## argument and nextflow as one argument per whitespace separated token
deopts_raw <- if (argsLen > 7) args[8:argsLen] else character(0)

print(args)

## FUNCS
libp <- paste0(gsub("/bin/conda", "/envs/monsda", Sys.getenv("CONDA_EXE")), "/share/MONSDA/scripts/lib/_lib.R")
if (!file.exists(libp)) {
    postde_args <- commandArgs(trailingOnly = FALSE)
    postde_file <- sub("^--file=", "", postde_args[grep("^--file=", postde_args)])
    if (length(postde_file) == 0) {
        stop("Rscript --file path not found in commandArgs")
    }
    libp <- file.path(dirname(normalizePath(postde_file)), "..", "..", "lib", "_lib.R")
}
source(libp)

de_opts <- parse_de_options(deopts_raw)
spike <- de_option(de_opts, "controlgenes")
## number of RUVg factors of unwanted variation: "auto" estimates it from the
## control genes (permutation test on the variance spectrum), an integer pins
## it, 0 switches the W covariates off and normalizes on the spike-in scale only
ruvk <- de_option(de_opts, "ruvk", "auto")
## model formula for the per-contrast fit, over condition/type/batch, e.g.
## design=batch+condition. Empty means MONSDA derives it from the annotation.
design_spec <- de_option(de_opts, "design")

postde_enabled <- Sys.getenv("MONSDA_POSTDE") == "1"
if (postde_enabled) {
    postde_args <- commandArgs(trailingOnly = FALSE)
    postde_file <- sub("^--file=", "", postde_args[grep("^--file=", postde_args)])
    if (length(postde_file) == 0) {
        stop("MONSDA_POSTDE=1 but Rscript --file path not found in commandArgs")
    }
    postde_lib <- file.path(dirname(normalizePath(postde_file)), "..", "PostDE", "export.R")
    if (!file.exists(postde_lib)) {
        stop(paste0("MONSDA_POSTDE=1 but export.R not found at ", postde_lib))
    }
    source(postde_lib)
    postde_bundles <- list()
    postde_entry_ids <- character(0)
    postde_cache_dir <- file.path(outdir, ".postde_cache")
    dir.create(postde_cache_dir, recursive = TRUE, showWarnings = FALSE)
    postde_cache_dir <- normalizePath(postde_cache_dir)
}

## set thread-usage
BPPARAM <- MulticoreParam(workers = availablecores)

### SCRIPT
## Annotation
sampleData_all <- as.data.frame(read.table(gzfile(anname), row.names = 1, check.names = FALSE, sep = "\t"))
colnames(sampleData_all) <- c("condition", "type", "batch")
sampleData_all$batch <- as.factor(sampleData_all$batch)
sampleData_all$type <- as.factor(sampleData_all$type)
sampleData_all$condition <- as.factor(sampleData_all$condition)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)
gtf_gene <- droplevels(subset(gtf.df, type == "gene"))

## Combinations of conditions
comparison <- parse_comparisons(cmp)
validate_comparisons(comparison, levels(sampleData_all$condition))

## check combi
if (combi == "none") {
    combi <- ""
}

## readin counttable
countData_all <- as.matrix(read.table(gzfile(countfile), header = T, row.names = 1, check.names = FALSE))

# Check if names are consistent
if (!all(rownames(sampleData_all) %in% colnames(countData_all))) {
    stop("Count file does not correspond to the annotation file")
}

comparison_objs <- list()

WD <- getwd()
setwd(outdir)
outdir_abs <- getwd()

# Calculate CPM and TPM
cpm_matrix <- calc_cpm(countData_all)
tpm_matrix <- calc_tpm(countData_all, gtf_gene)

# Write out CPM and TPM tables
write.table(add_gene_coordinates(cpm_matrix, rownames(cpm_matrix), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "cpm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
write.table(add_gene_coordinates(tpm_matrix, rownames(tpm_matrix), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "tpm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)

# Normalize by spike in if available
if (spike != "") {
    print("Spike-in used, data will be normalized to spike in separately")
    spiken <- spike
    setwd(WD)
    ctrlgenes <- trimws(readLines(spiken))
    ctrlgenes <- ctrlgenes[nzchar(ctrlgenes)]
    setwd(outdir)
    missing_ctrlgenes <- setdiff(ctrlgenes, rownames(countData_all))
    if (length(missing_ctrlgenes) > 0) {
        warning(paste0("Control genes not found in count table and will be ignored: ", paste(missing_ctrlgenes, collapse = ", ")))
    }
    ctrlgenes <- intersect(ctrlgenes, rownames(countData_all))
    if (length(ctrlgenes) == 0) {
        stop(paste0("None of the control genes in ", spiken, " are present in the count table, cannot normalize by spike-in"))
    }
    ctrl_idx_all <- rownames(countData_all) %in% ctrlgenes # for spike-in-derived size factors
    # Derive size factors from the spike-ins themselves (controlGenes); normCounts() from RUVg
    # only removes unwanted variation (W_i) and does not carry the spike-in scale.
    dds_sf <- DESeqDataSetFromMatrix(countData = as.matrix(countData_all), colData = sampleData_all, design = ~1)
    dds_sf <- estimateSizeFactors(dds_sf, controlGenes = ctrl_idx_all)
    # how many factors of unwanted variation the control genes actually support
    ruv_all <- ruv_estimate_k(countData_all, ctrlgenes, model.matrix(~condition, data = sampleData_all),
        condition = sampleData_all$condition, size_factors = sizeFactors(dds_sf), k = ruvk, label = "DataSet")
    ruv_write_k_report(ruv_all, paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "ruv_k_selection.tsv", sep = "_"))
    ruv_plot_k(ruv_all, paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "RUV_k_selection.png", sep = "_"))
    if (ruv_all$k > 0) {
        ruv_set_all <- RUVg(newSeqExpressionSet(as.matrix(countData_all)), ctrlgenes, k = ruv_all$k)
        sampleData_ruv <- cbind(sampleData_all, pData(ruv_set_all))
    } else {
        sampleData_ruv <- sampleData_all
    }
    # counts_norm is the spike-in size factor normalized matrix (the _norm branch); the RUVg
    # factors live in sampleData_ruv and enter the model as covariates, not as a count scale
    counts_norm <- as.data.frame(counts(dds_sf, normalized = TRUE))
    countData_clean <- countData_all %>% subset(!row.names(countData_all) %in% ctrlgenes) # removing spike-ins for standard analysis
    counts_norm_clean <- counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes) # removing spike-ins for standard analysis
     write.table(add_gene_coordinates(counts_norm, rownames(counts_norm), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "counts_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(counts_norm_clean, rownames(counts_norm_clean), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "counts_norm_clean.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(countData_clean, rownames(countData_clean), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "counts_clean.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(sampleData_ruv, gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "sampleData_ruv.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    # Calculate CPM and TPM
    cpm_norm_matrix <- calc_cpm(counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes))
    tpm_norm_matrix <- calc_tpm(counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes), gtf_gene)

    # Write out CPM and TPM tables
    write.table(add_gene_coordinates(cpm_norm_matrix, rownames(cpm_norm_matrix), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "cpm_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(tpm_norm_matrix, rownames(tpm_norm_matrix), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "tpm_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)

}
#### Now plot and print over-all comparisons

## Create design-table considering different types (paired, unpaired) and batches
design <- ~ 0 + condition
print(design)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData_all, colData = sampleData_all, design = design)

# filter low counts
smallestGroupSize <- min_group_size(sampleData_all$condition)
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# run for each pair of conditions
dds <- DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, betaPrior = FALSE)

rld <- rlogTransformation(dds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "PCA.png", sep = "_"), width=1900, height=1200, res=75)
print(DESeq2::plotPCA(vsd, intgroup = c("condition")) + geom_text_repel(aes(label = name), arrow = arrow(length = unit(0.01, "npc")), box.padding = 1, max.overlaps = 100, force = 5) + theme_bw()) # requires ggrepel)
dev.off()

# We also write the normalized counts to file
write.table(add_gene_coordinates(as.data.frame(assay(rld)), rownames(assay(rld)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "rld.tsv.gz", sep = "_")), sep = "\t", col.names = NA)
write.table(add_gene_coordinates(as.data.frame(assay(vsd)), rownames(assay(vsd)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "vsd.tsv.gz", sep = "_")), sep = "\t", col.names = NA)

# Here we choose blind so that the initial conditions setting does not influence the outcome, ie we want to see if the conditions cluster based purely on the individual datasets, in an unbiased way. According to the documentation, the rlogTransformation method that converts counts to log2 values is apparently better than the old varienceStabilisation method when the data size factors vary by large amounts.

par(mai = ifelse(1:4 <= 2, par("mai"), 0))
px <- counts(dds)[, 1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1, length(ord), length = 50)]
last <- ord[length(ord)]
vstcol <- c("blue", "black")
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "VST-and-log2.png", sep = "_"),width=400,height=350,res=45)
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type = "l", lty = 1, col = vstcol, xlab = "n", ylab = "f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n / s[1]))), fill = vstcol)
dev.off()

##############################
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hm_geom <- scaled_png_geometry(ncol(counts(dds)))
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap1.png", sep = "_"), width=hm_geom$width, height=hm_geom$height, res=hm_geom$res)
heatmap.2(counts(dds, normalized = TRUE)[select, ],
    col = hmcol,
    Rowv = FALSE, Colv = FALSE, scale = "none",
    dendrogram = "none", trace = "none", margin = c(hm_geom$margin, 8), cexCol = hm_geom$cex
)
dev.off()
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap2.png", sep = "_"), width=hm_geom$width, height=hm_geom$height, res=hm_geom$res)
heatmap.2(assay(rld)[select, ],
    col = hmcol,
    Rowv = FALSE, Colv = FALSE, scale = "none",
    dendrogram = "none", trace = "none", margin = c(hm_geom$margin, 8), cexCol = hm_geom$cex
)
dev.off()
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap3.png", sep = "_"), width=hm_geom$width, height=hm_geom$height, res=hm_geom$res)
heatmap.2(assay(vsd)[select, ],
    col = hmcol,
    Rowv = FALSE, Colv = FALSE, scale = "none",
    dendrogram = "none", trace = "none", margin = c(hm_geom$margin, 8), cexCol = hm_geom$cex
)
dev.off()

# The above shows heatmaps for 30 most highly expressed genes (not necessarily the biggest fold change). The data is of raw counts (left), regularized log transformation (center) and from variance stabilizing transformation (right) and you can clearly see the effect of the transformation has by shrinking the variance so that we don’t get the squish effect shown in the left hand graph.
##############################
# Sample-to-sample Pearson correlation of the rlog-transformed counts. A
# correlation matrix is directly interpretable (bounded -1..1) and, via the
# shared plot_sample_correlation_heatmap() helper, stays legible and keeps its
# color key from overlapping the plot regardless of sample count.
plot_sample_correlation_heatmap(
    assay(rld),
    paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap-samplebysample.png", sep = "_"),
    group = colData(dds)$condition,
    title = "Sample-to-sample Pearson correlation (rlog)"
)

### Now we run comparisons

for (contrast in comparison) {
    contrast_name <- contrast$name
    print(paste("Comparing ", contrast_name, sep = ""))

    # determine contrast
    A <- contrast$A
    B <- contrast$B

    # subset Datasets for pairwise comparison: metadata-based selection, B then A order
    sel <- select_contrast_samples(sampleData_all, countData_all, A, B)
    sampleData <- sel$sampleData
    countData <- sel$countData

    ## Design for this contrast: `design=` from the config if set, otherwise
    ## derived from the types and batches present in the two groups
    design <- build_de_design(design_spec, sampleData, A, B, contrast_name)
    print(design)

    ## Fit one spike-in branch (_norm or _ruv) and write its QC outputs. Identical to the
    ## plain branch below except for the suffix, the optional externally supplied size
    ## factors and the design, so the three readings stay directly comparable.
    fit_spike_branch <- function(counts_b, meta_b, design_b, size_factors, suffix) {
        d <- DESeqDataSetFromMatrix(countData = counts_b, colData = meta_b, design = design_b)
        if (!is.null(size_factors)) {
            sizeFactors(d) <- size_factors[colnames(d)]
        }
        d <- d[rowSums(counts(d) >= 10) >= min_group_size(meta_b$condition), ]
        d$condition <- droplevels(d$condition)
        d$condition <- relevel(d$condition, ref = B[[1]])
        d <- DESeq(d, parallel = TRUE, BPPARAM = BPPARAM, betaPrior = FALSE)
        rld_b <- rlogTransformation(d, blind = FALSE)
        vsd_b <- varianceStabilizingTransformation(d, blind = FALSE)
        png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", paste0("PCA", suffix, ".png"), sep = "_"), width=1900, height=1200, res=75)
        print(DESeq2::plotPCA(vsd_b, intgroup = c("condition")) + geom_text_repel(aes(label = name), arrow = arrow(length = unit(0.02, "npc")), box.padding = .5) + theme_bw()) # requires ggrepel)
        dev.off()
        write.table(add_gene_coordinates(as.data.frame(assay(rld_b)), rownames(assay(rld_b)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", paste0("rld", suffix, ".tsv.gz"), sep = "_")), sep = "\t", col.names = NA)
        write.table(add_gene_coordinates(as.data.frame(assay(vsd_b)), rownames(assay(vsd_b)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", paste0("vsd", suffix, ".tsv.gz"), sep = "_")), sep = "\t", col.names = NA)
        list(dds = d, vsd = vsd_b)
    }

    ## Results, tables and figures of one spike-in branch. Byte-identical in structure to
    ## the plain branch below, only the suffix and the fitted object differ.
    emit_spike_branch <- function(dds_b, meta_b, counts_b, vsd_b, design_b, suffix) {
        listname <- paste0(contrast_name, suffix)
        res <- results(dds_b, contrast = c("condition", A, B), parallel = TRUE, BPPARAM = BPPARAM)
        resn <- res
        res_shrink <- lfcShrink(dds = dds_b, coef = paste("condition", A, "vs", B, sep = "_"), res = res, type = "apeglm")
        check_shrinkage(res, res_shrink, label = listname)

        if (postde_enabled) {
            postde_capture(engine = "deseq2", id = listname, A = A, B = B, normalized = TRUE, metadata = meta_b, counts = counts_b, expression = assay(vsd_b), formula = design_b, results = postde_results_deseq2(res), mean_scale = "baseMean")
        }

        resOrdered <- format_deseq2_results(res_shrink[order(res_shrink$log2FoldChange), ], gtf_gene, shrink = TRUE, center = TRUE)
        pdf(file = paste("Figures/DE", "DESEQ2", combi, contrast_name, paste0("figure_Volcano", suffix, ".pdf"), sep = "_"), width = 15, height = 10)
        print(EnhancedVolcano(resOrdered,
            lab = resOrdered$Gene, x = "log2FoldChange", y = "padj",
            title = paste0(contrast_name, suffix, "_p005_lfc15"),
            pCutoff = 0.05, FCcutoff = 1.5, pointSize = 3.0, labSize = 5.0, colAlpha = .3,
            legendLabels = c("Not sig.", "Log (base 2) FC", "p-value", "p-value & Log (base 2) FC"),
            legendPosition = "right", legendLabSize = 10, legendIconSize = 5.0,
            drawConnectors = TRUE, widthConnectors = 0.75
        ))
        dev.off()
        write.table(as.data.frame(apply(resOrdered, 2, as.character)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", paste0("results", suffix, ".tsv.gz"), sep = "_")), sep = "\t", row.names = FALSE, quote = F)

        res <- resn[order(resn$log2FoldChange), ]
        res$Gene_ID <- rownames(res)
        res <- format_deseq2_results(res, gtf_gene, shrink = FALSE, center = TRUE)
        write.table(as.data.frame(apply(res, 2, as.character)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", paste0("results", suffix, "_noshrink.tsv.gz"), sep = "_")), sep = "\t", row.names = FALSE, quote = F)

        png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", paste0("MA", suffix, ".png"), sep = "_"), width=1900, height=1200, res=300)
        DESeq2::plotMA(res_shrink)
        dev.off()
        resn
    }

    # Normalize by spike in if available
    if (spike != "") {
        print("Spike-in used, data will be normalized to spike in separately")
        spiken <- spike
        setwd(WD)
        ctrlgenes <- trimws(readLines(spiken))
        ctrlgenes <- ctrlgenes[nzchar(ctrlgenes)]
        setwd(outdir)
        missing_ctrlgenes <- setdiff(ctrlgenes, rownames(countData))
        if (length(missing_ctrlgenes) > 0) {
            warning(paste0("Control genes not found in count table and will be ignored: ", paste(missing_ctrlgenes, collapse = ", ")))
        }
        ctrlgenes <- intersect(ctrlgenes, rownames(countData))
        if (length(ctrlgenes) == 0) {
            stop(paste0("None of the control genes in ", spiken, " are present in the count table, cannot normalize by spike-in"))
        }

        ctrl_idx <- rownames(countData) %in% ctrlgenes # spike-ins are never tested themselves
        counts_endo <- as.matrix(countData)[!ctrl_idx, , drop = FALSE]

        # Spike-in derived size factors (median of ratios over the control genes): the
        # normalization of the _norm branch, and the reference the k diagnostics compare the
        # factors of unwanted variation against.
        sf_spike <- spikein_size_factors(as.matrix(countData)[ctrl_idx, , drop = FALSE])
        write_scaling_log(sf_spike, colnames(countData), paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "scaling.log", sep = "_"))

        # how many factors of unwanted variation the control genes support for this contrast
        ruv_sel <- ruv_estimate_k(countData, ctrlgenes, stats::model.matrix(design, data = sampleData),
            condition = sampleData$condition, size_factors = sf_spike, k = ruvk, label = contrast_name, contrast = c(A, B))
        ruv_write_k_report(ruv_sel, paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "ruv_k_selection.tsv", sep = "_"))
        ruv_plot_k(ruv_sel, paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "RUV_k_selection.png", sep = "_"))

        if (ruv_sel$k > 0) {
            ruv_set <- RUVg(newSeqExpressionSet(as.matrix(countData)), ctrlgenes, k = ruv_sel$k)
            sampleData_ruv <- cbind(sampleData, pData(ruv_set))
            design_ruv <- ruv_add_w_terms(design, colnames(pData(ruv_set)))  # Last argument is variable of interest
        } else {
            sampleData_ruv <- sampleData
            design_ruv <- design
        }
        countData <- countData %>% subset(!row.names(countData) %in% ctrlgenes) # removing spike-ins for standard analysis
        print(paste0("Design of the _ruv branch: ", paste(all.vars(design_ruv), collapse = ", ")))

        # _norm: the spike-ins set the scale. Size factors from the control genes, design
        # unchanged, so its fold changes are absolute (per cell) rather than per library.
        print("Spike-in branch _norm: size factors from the control genes, no RUV covariates")
        branch_norm <- fit_spike_branch(counts_endo, sampleData, design, sf_spike, "_norm")
        dds_norm <- branch_norm$dds
        vsd_norm <- branch_norm$vsd

        # _ruv: RUVg exactly as the RUVSeq vignette prescribes it. The W_i covariates enter
        # the design and normalization is DESeq2's own median of ratios over the endogenous
        # genes; the spike-in size factors are deliberately NOT applied on top. Both are
        # estimated from the same control genes, so combining them would correct the same
        # quantity twice and leave the result depending on which of the two wins.
        print("Spike-in branch _ruv: RUVg factors as covariates, standard median-of-ratios normalization")
        branch_ruv <- fit_spike_branch(counts_endo, sampleData_ruv, design_ruv, NULL, "_ruv")
        dds_ruv <- branch_ruv$dds
        vsd_ruv <- branch_ruv$vsd
    }

    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = design)

    # filter low counts
    smallestGroupSize <- min_group_size(sampleData$condition)
    keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    dds <- dds[keep, ]

    # drop unused samples
    dds$condition <- droplevels(dds$condition)

    # relevel to base condition B
    dds$condition <- relevel(dds$condition, ref = B[[1]])

    # run for each pair of conditions
    dds <- DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, betaPrior = FALSE)
    print(resultsNames(dds))

    # Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

    rld <- rlogTransformation(dds, blind = FALSE)
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

    png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "PCA.png", sep = "_"), width=1900, height=1200, res=75)
    print(DESeq2::plotPCA(vsd, intgroup = c("condition")) + geom_text_repel(aes(label = name), arrow = arrow(length = unit(0.02, "npc")), box.padding = .5) + theme_bw()) # requires ggrepel)
    # DESeq2::plotPCA(rld, intgroup=c('condition')) + geom_text(aes(label = name), position = position_nudge(y = 2))
    dev.off()

    # We also write the normalized counts to file
    write.table(add_gene_coordinates(as.data.frame(assay(rld)), rownames(assay(rld)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "rld.tsv.gz", sep = "_")), sep = "\t", col.names = NA)
    write.table(add_gene_coordinates(as.data.frame(assay(vsd)), rownames(assay(vsd)), gtf_gene), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "vsd.tsv.gz", sep = "_")), sep = "\t", col.names = NA)

    tryCatch({
        # initialize empty objects
        res <- ""
        resOrdered <- ""
        res <- results(dds, contrast = c("condition", A, B), parallel = TRUE, BPPARAM = BPPARAM)
        resn <- res
        res_shrink <- lfcShrink(dds = dds, coef = paste("condition", A, "vs", B, sep = "_"), res = res, type = "apeglm")
        check_shrinkage(res, res_shrink, label = contrast_name)

        # add comp object to list for image
        if (postde_enabled) {
            postde_capture(engine = "deseq2", id = contrast_name, A = A, B = B, normalized = FALSE, metadata = sampleData, counts = countData, expression = assay(vsd), formula = design, results = postde_results_deseq2(res), mean_scale = "baseMean")
        }
        comparison_objs[[contrast_name]] <- res

        # sort and output
        resOrdered <- res_shrink[order(res_shrink$log2FoldChange), ]

        # # Add gene names  (check how gene_id col is named )
        resOrdered <- format_deseq2_results(resOrdered, gtf_gene, shrink = TRUE)

        # plotVolcano
        pdf(
            file = paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure_Volcano.pdf", sep = "_"), width = 15, height = 10
        )
        print(EnhancedVolcano(resOrdered,
            lab = resOrdered$Gene,
            x = "log2FoldChange",
            y = "padj",
            title = paste0(contrast_name, "_p005_lfc15", sep = ""),
            pCutoff = 0.05,
            FCcutoff = 1.5,
            pointSize = 3.0,
            labSize = 5.0,
            colAlpha = .3,
            legendLabels = c(
                "Not sig.", "Log (base 2) FC", "p-value",
                "p-value & Log (base 2) FC"
            ),
            legendPosition = "right",
            legendLabSize = 10,
            legendIconSize = 5.0,
            drawConnectors = TRUE,
            widthConnectors = 0.75
        ))
        dev.off()

        resOrdered <- as.data.frame(apply(resOrdered, 2, as.character))

        # write the table to a tsv file
        write.table(as.data.frame(resOrdered), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

        # sort and output
        res <- resn[order(resn$log2FoldChange), ]

        res$Gene <- vapply(rownames(res), function(x) {
            get_gene_name(x, gtf_gene)
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
        res$Gene_ID <- rownames(res)
        res <- format_deseq2_results(res, gtf_gene, shrink = FALSE)
        res <- as.data.frame(apply(res, 2, as.character))

        write.table(as.data.frame(res), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results_noshrink.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

        # plotMA
        png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "MA.png", sep = "_"), width=1900, height=1200, res=75)
        DESeq2::plotMA(res_shrink)
        dev.off()

        if (spike != "") { # DE runs for the two spike-in branches
            # _norm: spike-in size factors, plain design  -> absolute (per cell) fold changes
            comparison_objs[[paste0(contrast_name, "_norm")]] <- emit_spike_branch(dds_norm, sampleData, counts_endo, vsd_norm, design, "_norm")
            # _ruv: RUVg covariates, standard normalization -> unwanted variation removed
            comparison_objs[[paste0(contrast_name, "_ruv")]] <- emit_spike_branch(dds_ruv, sampleData_ruv, counts_endo, vsd_ruv, design_ruv, "_ruv")
        }
        # cleanup
        rm(res, resOrdered)
        print(paste("cleanup done for ", contrast_name, sep = ""))
    })
}


##############################

if (postde_enabled) {
    postde_write(outdir_abs, combi, "deseq2")
}

save.image(file = paste("DE", "DESEQ2", combi, "SESSION.gz", sep = "_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
