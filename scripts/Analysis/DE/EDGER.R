suppressPackageStartupMessages({
    require(utils)
    require(BiocParallel)
    require(edgeR)
    require(rtracklayer)
    require(RUVSeq)
    require(dplyr)
    require(GenomeInfoDb)
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

## plotMDS needs at least 3 samples/columns, fall back to a placeholder plot otherwise
safe_plotMDS <- function(dgeobj, ...) {
    if (ncol(dgeobj) < 3) {
        plot.new()
        text(0.5, 0.5, paste0("MDS plot skipped: only ", ncol(dgeobj), " sample(s), need at least 3"), cex = 0.8)
    } else {
        plotMDS(dgeobj, ...)
    }
}

## When a contrast has no residual degrees of freedom (e.g. two groups with a single sample each),
## edgeR cannot estimate dispersion from the data at all (estimateDisp/glmQLFit fail with
## "No residual df"/"NA dispersions not allowed"/"Could not estimate prior df"). Fall back to a
## fixed, typical BCV and edgeR's classic exactTest for that contrast only, same as EDGER_no_rep.R.
FALLBACK_BCV <- 0.4
no_residual_df <- function(design) {
    (nrow(design) - qr(design)$rank) <= 0
}

### SCRIPT
print(paste("Run EdgeR DE with ", availablecores, " cores", sep = ""))

# set thread-usage
BPPARAM <- MulticoreParam(workers = availablecores)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)
gtf_gene <- droplevels(subset(gtf.df, type == "gene"))

## Annotation
sampleData_all <- as.data.frame(read.table(gzfile(anname), row.names = 1, check.names = FALSE, sep = "\t"))
colnames(sampleData_all) <- c("condition", "type", "batch")
sampleData_all$condition <- as.factor(sampleData_all$condition)
sampleData_all$batch <- as.factor(sampleData_all$batch)
sampleData_all$type <- as.factor(sampleData_all$type)
samples <- rownames(sampleData_all)

## Combinations of conditions
comparison <- parse_comparisons(cmp)
validate_comparisons(comparison, levels(sampleData_all$condition))

## check combi
if (combi == "none") {
    combi <- ""
}

## readin counttable
countData_all <- read.table(countfile, header = TRUE, row.names = 1, check.names = FALSE)

# Check if names are consistent
if (!all(rownames(sampleData_all) %in% colnames(countData_all))) {
    stop("Count file does not correspond to the annotation file")
}

comparison_objs <- list()

WD <- getwd()
setwd(outdir)
outdir_abs <- getwd()

## Calculate CPM and TPM
cpm_matrix <- calc_cpm(countData_all)
tpm_matrix <- calc_tpm(countData_all, gtf_gene)

## Write out CPM and TPM tables
write.table(add_gene_coordinates(cpm_matrix, rownames(cpm_matrix), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "cpm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
write.table(add_gene_coordinates(tpm_matrix, rownames(tpm_matrix), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "tpm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)

## Normalize by spike in if available
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
    ctrl_idx_all <- rownames(countData_all) %in% ctrlgenes # for spike-in-derived normalization factors
    # Derive the scale from the spike-ins themselves (median of ratios over the control
    # genes); counts_norm below is that _norm scale. RUVg's factors are covariates, not a
    # count scale, and live in sampleData_ruv.
    size_factor_all <- spikein_size_factors(as.matrix(countData_all)[ctrl_idx_all, , drop = FALSE])
    # how many factors of unwanted variation the control genes actually support
    ruv_all <- ruv_estimate_k(countData_all, ctrlgenes, model.matrix(~condition, data = sampleData_all),
        condition = sampleData_all$condition, size_factors = size_factor_all, k = ruvk, label = "DataSet")
    ruv_write_k_report(ruv_all, paste("Tables/DE", "EDGER", combi, "DataSet", "table", "ruv_k_selection.tsv", sep = "_"))
    ruv_plot_k(ruv_all, paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "RUV_k_selection.png", sep = "_"))
    if (ruv_all$k > 0) {
        ruv_set_all <- RUVg(newSeqExpressionSet(as.matrix(countData_all)), ctrlgenes, k = ruv_all$k)
        sampleData_ruv <- cbind(sampleData_all, pData(ruv_set_all))
    } else {
        sampleData_ruv <- sampleData_all
    }
    counts_norm <- as.data.frame(sweep(as.matrix(countData_all), 2, size_factor_all, FUN = "/"))
    countData_clean <- countData_all %>% subset(!row.names(countData_all) %in% ctrlgenes) # removing spike-ins for standard analysis
    counts_norm_clean <- counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes) # removing spike-ins for standard analysis
    write.table(add_gene_coordinates(counts_norm, rownames(counts_norm), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "counts_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(counts_norm_clean, rownames(counts_norm_clean), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "counts_norm_clean.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(countData_clean, rownames(countData_clean), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "counts_clean.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(sampleData_ruv, gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "sampleData_ruv.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    # Calculate CPM and TPM
    cpm_norm_matrix <- calc_cpm(counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes))
    tpm_norm_matrix <- calc_tpm(counts_norm %>% subset(!row.names(counts_norm) %in% ctrlgenes), gtf_gene)

    # Write out CPM and TPM tables
    write.table(add_gene_coordinates(cpm_norm_matrix, rownames(cpm_norm_matrix), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "cpm_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
    write.table(add_gene_coordinates(tpm_norm_matrix, rownames(tpm_norm_matrix), gtf_gene), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "tpm_norm.tsv.gz", sep = "_")), sep = "\t", col.names = NA, quote = FALSE)
}

## Create design-table considering different types (paired, unpaired) and batches
des <- ~ 0 + condition
design <- model.matrix(des, data = sampleData_all)
colnames(design) <- levels(sampleData_all$condition)
print(design)

genes <- rownames(countData_all)
samples <- rownames(sampleData_all)

dge <- DGEList(counts = countData_all, group = sampleData_all$condition, samples = samples, genes = genes)

## filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

## normalize with TMM
dge <- calcNormFactors(dge, method = "TMM", BPPARAM = BPPARAM)

## create file normalized table
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm), 1:ncol(tmm) - 1)]

tmm <- add_gene_coordinates(tmm, tmm$ID, gtf_gene, after = "ID")
write.table(as.data.frame(tmm), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "AllConditionsNormalized.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

## geometry (canvas size/res/point&label cex) scaled to the number of samples,
## so DataSet figures stay legible whether there are a handful or hundreds
ds_geom <- scaled_png_geometry(ncol(dge))

## sample-to-sample Pearson correlation on log-CPM, consistent with the
## correlation plot produced by the DESEQ2 DataSet figures
plot_sample_correlation_heatmap(
    cpm(dge, log = TRUE),
    paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsCorrelation.png", sep = "_"),
    group = dge$samples$group,
    title = "Sample-to-sample Pearson correlation (log2 CPM)"
)

## create file MDS-plot with and without summarized replicates
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsMDS.png", sep = "_")
png(out, width=ds_geom$width, height=ds_geom$height, res=ds_geom$res)
safe_plotMDS(dge, col = as.numeric(dge$samples$group), pch = 16, cex = ds_geom$cex)
if (ncol(dge) >= 3) {
    legend("topright", legend = levels(as.factor(dge$samples$group)), col = seq_along(levels(as.factor(dge$samples$group))), pch = 16, cex = 0.8)
}
dev.off()

## estimate Dispersion
dge <- estimateDisp(dge, design, robust = TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsBCV.png", sep = "_")
png(out, width=ds_geom$width, height=ds_geom$height, res=ds_geom$res)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust = TRUE)

## create file quasi-likelihood-dispersion-plot
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsQLDisp.png", sep = "_")
png(out, width=ds_geom$width, height=ds_geom$height, res=ds_geom$res)
plotQLDisp(fit)
dev.off()

## Analyze according to comparison groups
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
    sampleData$condition <- relevel(sampleData$condition, ref = B)

    samples <- rownames(sampleData)
    ## name types and levels for design
    bl <- sapply("batch", paste0, levels(sampleData$batch)[1:length(levels(sampleData$batch)) - 1])
    tl <- sapply("type", paste0, levels(sampleData$type)[1:length(levels(sampleData$type)) - 1])

    ## Design for this contrast: `design=` from the config if set, otherwise
    ## derived from the types and batches present in the two groups
    des <- build_de_design(design_spec, sampleData, A, B, contrast_name)
    design <- model.matrix(des, data = sampleData)
    print(design)

    ## Results, tables and figures of one spike-in branch. Byte-identical in structure to
    ## the plain branch below, only the suffix and the fitted object differ.
    emit_spike_branch <- function(branch, meta_b, counts_b, des_b, suffix) {
        contrast_col <- paste0("condition", A)
        if (!contrast_col %in% colnames(branch$design)) {
            stop(paste0("Spike-in design ", suffix, " has no column '", contrast_col, "' (columns: ", paste(colnames(branch$design), collapse = ", "), ")"))
        }
        AvsB <- setNames(numeric(ncol(branch$design)), colnames(branch$design))
        AvsB[contrast_col] <- 1
        qlf <- glmQLFTest(branch$fit, contrast = AvsB) ## glm quasi-likelihood-F-Test

        if (postde_enabled) {
            postde_capture(engine = "edger", id = paste0(contrast_name, suffix), A = A, B = B, normalized = TRUE, metadata = meta_b, counts = counts_b, expression = cpm(branch$dge, log = TRUE), formula = des_b, results = postde_results_edger(qlf), mean_scale = "logCPM")
        }

        res <- format_edger_results(qlf$table, gtf_gene, center = TRUE)
        pdf(file = paste("Figures/DE", "EDGER", combi, contrast_name, paste0("figure_Volcano", suffix, ".pdf"), sep = "_"), width = 15, height = 10)
        print(EnhancedVolcano(res,
            lab = res$Gene, x = "log2FoldChange", y = "padj",
            title = paste0(contrast_name, suffix, "_p005_lfc15"),
            pCutoff = 0.05, FCcutoff = 1.5, pointSize = 3.0, labSize = 5.0, colAlpha = .3,
            legendLabels = c("Not sig.", "Log (base 2) FC", "p-value", "p-value & Log (base 2) FC"),
            legendPosition = "right", legendLabSize = 10, legendIconSize = 5.0,
            drawConnectors = TRUE, widthConnectors = 0.75
        ))
        dev.off()
        write.table(as.data.frame(apply(res, 2, as.character)), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", paste0("results", suffix, ".tsv.gz"), sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        for (sorting in c("logFC", "PValue")) {
            tops <- format_edger_results(topTags(qlf, n = nrow(qlf$table), sort.by = sorting)$table, gtf_gene, center = TRUE)
            label <- if (sorting == "logFC") "resultsLogFCsorted" else "resultsPValueSorted"
            write.table(as.data.frame(apply(tops, 2, as.character)), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", paste0(label, suffix, ".tsv.gz"), sep = "_")), sep = "\t", quote = F, row.names = FALSE)
        }

        png(paste("Figures/DE", "EDGER", combi, contrast_name, "figure", paste0("MD", suffix, ".png"), sep = "_"), width=1900, height=1200, res=300)
        print(plotMD(qlf, main = paste0(contrast_name, suffix)))
        abline(h = c(-1, 1), col = "blue")
        dev.off()
        qlf
    }

    ## Fit one spike-in branch (_norm or _ruv) and write its QC outputs. Identical to the
    ## plain branch below except for the suffix, the normalization and the design, so the
    ## three readings stay directly comparable.
    fit_spike_branch <- function(counts_b, meta_b, des_b, size_factors, suffix) {
        design_b <- model.matrix(des_b, data = meta_b)
        print(paste0("Design of the ", suffix, " branch: ", paste(colnames(design_b), collapse = ", ")))
        d <- DGEList(counts = counts_b, group = meta_b$condition, samples = samples, genes = rownames(counts_b))
        d <- d[filterByExpr(d), , keep.lib.sizes = FALSE]
        d$samples$group <- relevel(d$samples$group, ref = B[[1]])
        if (is.null(size_factors)) {
            d <- calcNormFactors(d, method = "TMM", BPPARAM = BPPARAM)
        } else {
            d <- set_spikein_offsets(d, size_factors[colnames(d)])
        }

        tmm_b <- as.data.frame(cpm(d))
        colnames(tmm_b) <- t(d$samples$samples)
        tmm_b$ID <- d$genes$genes
        tmm_b <- tmm_b[c(ncol(tmm_b), 1:ncol(tmm_b) - 1)]
        tmm_b <- add_gene_coordinates(tmm_b, tmm_b$ID, gtf_gene, after = "ID")
        write.table(as.data.frame(tmm_b), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", paste0("Normalized", suffix, ".tsv.gz"), sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        png(paste("Figures/DE", "EDGER", combi, contrast_name, "figure", paste0("MDS", suffix, ".png"), sep = "_"), width=1900, height=1200, res=300)
        print(safe_plotMDS(d, col = as.numeric(d$samples$group), cex = 1))
        dev.off()

        if (no_residual_df(design_b)) {
            stop(paste0("Spike-in design ", suffix, " for contrast ", contrast_name, " has no residual degrees of freedom (", nrow(design_b), " samples, design rank ", qr(design_b)$rank, "). The nuisance-adjusted model cannot be fit; refusing to fall back to an unadjusted test that ignores the covariates."))
        }
        d <- estimateDisp(d, design_b, robust = TRUE)

        png(paste("Figures/DE", "EDGER", combi, contrast_name, "figure", paste0("BCV", suffix, ".png"), sep = "_"), width=1900, height=1200, res=300)
        print(plotBCV(d))
        dev.off()

        f <- glmQLFit(d, design_b, robust = TRUE)

        png(paste("Figures/DE", "EDGER", combi, contrast_name, "figure", paste0("QLDisp", suffix, ".png"), sep = "_"), width=1900, height=1200, res=300)
        print(plotQLDisp(f))
        dev.off()
        list(dge = d, fit = f, design = design_b)
    }

    ## check genes and spike-ins
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

        ## Spike-in derived scale (median of ratios over the control genes): the normalization
        ## of the _norm branch, and the reference the k diagnostics compare the factors of
        ## unwanted variation against.
        size_factor_spike <- spikein_size_factors(as.matrix(countData)[ctrl_idx, , drop = FALSE])
        write_scaling_log(size_factor_spike, colnames(countData), paste("Tables/DE", "EDGER", combi, contrast_name, "table", "scaling.log", sep = "_"))

        ## how many factors of unwanted variation the control genes support for this contrast
        ruv_sel <- ruv_estimate_k(countData, ctrlgenes, design,
            condition = sampleData$condition, size_factors = size_factor_spike, k = ruvk, label = contrast_name, contrast = c(A, B))
        ruv_write_k_report(ruv_sel, paste("Tables/DE", "EDGER", combi, contrast_name, "table", "ruv_k_selection.tsv", sep = "_"))
        ruv_plot_k(ruv_sel, paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "RUV_k_selection.png", sep = "_"))

        if (ruv_sel$k > 0) {
            ruv_set <- RUVg(newSeqExpressionSet(as.matrix(countData)), ctrlgenes, k = ruv_sel$k)
            sampleData_ruv <- cbind(sampleData, pData(ruv_set))
            des_ruv <- ruv_add_w_terms(des, colnames(pData(ruv_set)))  # Last argument is variable of interest
        } else {
            sampleData_ruv <- sampleData
            des_ruv <- des
        }
        countData <- countData %>% subset(!row.names(countData) %in% ctrlgenes) # removing spike-ins for standard analysis

        ## _norm: the spike-ins set the scale. Offsets from the control genes, design
        ## unchanged, so its fold changes are absolute (per cell) rather than per library.
        print("Spike-in branch _norm: offsets from the control genes, no RUV covariates")
        branch_norm <- fit_spike_branch(counts_endo, sampleData, des, size_factor_spike, "_norm")
        dge_norm <- branch_norm$dge
        fit_norm <- branch_norm$fit
        design_norm <- branch_norm$design

        ## _ruv: RUVg exactly as the RUVSeq vignette prescribes it. The W_i covariates enter
        ## the design and normalization is edgeR's own TMM over the endogenous genes; the
        ## spike-in offsets are deliberately NOT applied on top. Both are estimated from the
        ## same control genes, so combining them would correct the same quantity twice and
        ## leave the result depending on which of the two wins.
        print("Spike-in branch _ruv: RUVg factors as covariates, standard TMM normalization")
        branch_ruv <- fit_spike_branch(counts_endo, sampleData_ruv, des_ruv, NULL, "_ruv")
        dge_ruv <- branch_ruv$dge
        fit_ruv <- branch_ruv$fit
        design_ruv <- branch_ruv$design
    }

    # Same without spike-in normalization
    genes <- rownames(countData)
    dge <- DGEList(counts = countData, group = sampleData$condition, samples = samples, genes = genes)

    ## filter low counts
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]

    # relevel to base condition B
    dge$samples$group <- relevel(dge$samples$group, ref = B[[1]])

    ## normalize with TMM
    dge <- calcNormFactors(dge, method = "TMM", BPPARAM = BPPARAM)

    ## create file normalized table
    tmm <- as.data.frame(cpm(dge))
    colnames(tmm) <- t(dge$samples$samples)
    tmm$ID <- dge$genes$genes
    tmm <- tmm[c(ncol(tmm), 1:ncol(tmm) - 1)]

    tmm <- add_gene_coordinates(tmm, tmm$ID, gtf_gene, after = "ID")
    write.table(as.data.frame(tmm), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "Normalized.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

    ## create file MDS-plot with and without summarized replicates
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "MDS.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
    print(safe_plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1))
    dev.off()

    no_rep <- no_residual_df(design)
    if (no_rep) {
        print(paste0("No residual degrees of freedom for contrast ", contrast_name, ": at least one group has no replicates, dispersion cannot be estimated from the data. Falling back to a fixed BCV of ", FALLBACK_BCV, " (edgeR exactTest) for this contrast."))
    } else {
        ## estimate Dispersion
        dge <- estimateDisp(dge, design, robust = TRUE)
    }

    ## create file BCV-plot - visualizing estimated dispersions
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "BCV.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
    if (no_rep) {
        plot.new()
        text(0.5, 0.5, paste0("BCV plot skipped: no replicates, using fixed BCV = ", FALLBACK_BCV), cex = 0.8)
    } else {
        print(plotBCV(dge))
    }
    dev.off()

    if (!no_rep) {
        ## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
        fit <- glmQLFit(dge, design, robust = TRUE)
    }

    ## create file quasi-likelihood-dispersion-plot
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "QLDisp.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
    if (no_rep) {
        plot.new()
        text(0.5, 0.5, "QLDisp plot skipped: no replicates", cex = 0.8)
    } else {
        print(plotQLDisp(fit))
    }
    dev.off()

    tryCatch({
        # determine contrast, only for complex cases, not needed for our pairwise comparisons now
        # A <- strsplit(contrast_groups[[1]][1], "\\+")
        # B <- strsplit(contrast_groups[[1]][2], "\\+")
        # minus <- 1/length(A[[1]])*(-1)
        # plus <- 1/length(B[[1]])
        # contrast <- cbind(integer(dim(design)[2]), colnames(design))
        # for(i in A[[1]]){
        #    contrast[which(contrast[,2]==i)]<- minus
        # }
        # for(i in B[[1]]){
        #    contrast[which(contrast[,2]==i)]<- plus
        # }
        # contrast <- as.numeric(contrast[,1])
        
        if (no_rep) {
            qlf <- exactTest(dge, pair = c(B[1], A[1]), dispersion = FALLBACK_BCV^2, prior.count = 2) ## no replicates: fixed-BCV exact test
            qlf$table$F <- NA
            qlf$table <- qlf$table[, c("logFC", "logCPM", "F", "PValue")]
        } else {
            AvsB <- makeContrasts(TreatvsUntreat = paste("condition", A, sep = ""), levels = design)
            qlf <- glmQLFTest(fit, contrast = AvsB) ## glm quasi-likelihood-F-Test
        }
        # add comp object to list for image
        if (postde_enabled) {
            postde_capture(engine = "edger", id = contrast_name, A = A, B = B, normalized = FALSE, metadata = sampleData, counts = countData, expression = cpm(dge, log = TRUE), formula = des, results = postde_results_edger(qlf), mean_scale = "logCPM")
        }
        comparison_objs[[contrast_name]] <- qlf

        # # Add gene names  (check how gene_id col is named )
        res <- format_edger_results(qlf$table, gtf_gene)

        # plotVolcano
        pdf(
            file = paste("Figures/DE", "EDGER", combi, contrast_name, "figure_Volcano.pdf", sep = "_"), width = 15, height = 10
        )
        print(EnhancedVolcano(res,
            lab = res$Gene,
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

        res <- as.data.frame(apply(res, 2, as.character))

        # create results table
        write.table(as.data.frame(res), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        # create sorted results Tables
        tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
        tops <- format_edger_results(tops$table, gtf_gene)
        tops <- as.data.frame(apply(tops, 2, as.character))
        write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsLogFCsorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "PValue")
        tops <- format_edger_results(tops$table, gtf_gene)
        tops <- as.data.frame(apply(tops, 2, as.character))
        write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsPValueSorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        ## plot lFC vs CPM
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "MD.png", sep = "_")
        png(out, width=1900, height=1200, res=300)
        print(plotMD(qlf, main = contrast_name))
        abline(h = c(-1, 1), col = "blue")
        dev.off()

        if (spike != "") { # DE runs for the two spike-in branches
            # _norm: spike-in offsets, plain design   -> absolute (per cell) fold changes
            comparison_objs[[paste0(contrast_name, "_norm")]] <- emit_spike_branch(branch_norm, sampleData, counts_endo, des, "_norm")
            # _ruv: RUVg covariates, TMM normalization -> unwanted variation removed
            comparison_objs[[paste0(contrast_name, "_ruv")]] <- emit_spike_branch(branch_ruv, sampleData_ruv, counts_endo, des_ruv, "_ruv")
        }

        # cleanup
        rm(qlf, res, tops)
        print(paste("cleanup done for ", contrast_name, sep = ""))
    }, error = function(e) {
        message(paste0("Error while processing contrast ", contrast_name, ": ", conditionMessage(e), ". Skipping this contrast."))
    })
}


if (postde_enabled) {
    postde_write(outdir_abs, combi, "edger")
}

save.image(file = paste("DE_EDGER", combi, "SESSION.gz", sep = "_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
