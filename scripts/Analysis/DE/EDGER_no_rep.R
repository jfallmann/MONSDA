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
bcv <- if (argsLen > 7) as.numeric(args[8]) else 0.2
spike <- if (argsLen > 8) args[9] else ""
if (!exists(quote(bcv))){
    bcv <- 0.2
}

print(args)
print(paste0("Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates. You selected ", bcv, sep = ""))

## FUNCS
libp <- paste0(gsub("/bin/conda", "/envs/monsda", Sys.getenv("CONDA_EXE")), "/share/MONSDA/scripts/lib/_lib.R")
source(libp)

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

# Same for all samples without design specific normalization
## name types and levels for design

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

## create dummy file MDS-plot with and without summarized replicates
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsMDS.png", sep = "_")
png::writePNG(array(0, dim = c(1,1,4)), out)

## create dummy file BCV-plot - visualizing estimated dispersions
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsBCV.png", sep = "_")
png::writePNG(array(0, dim = c(1,1,4)), out)

## create dummy file quasi-likelihood-dispersion-plot
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsQLDisp.png", sep = "_")
png::writePNG(array(0, dim = c(1,1,4)), out)

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

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == condition)$type)) > 1 | length(unique(subset(sampleData, B == condition)$type)) > 1) {
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1) {
            des <- ~ type + batch + condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), tl, bl)
        } else {
            des <- ~ type + condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(condition), tl)
        }
    } else {
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1) {
            des <- ~ batch + condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), bl)
        } else {
            des <- ~condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- levels(sampleData$condition)
        }
    }
    print(design)

    ## check genes and spike-ins
    if (spike != "") {
        print("Spike-in used, data will be normalized to spike in separately")
        spiken <- strsplit(spike, "=")[[1]][2]
        setwd(WD)
        ctrlgenes <- readLines(spiken)
        setwd(outdir)
        counts_norm <- RUVg(newSeqExpressionSet(as.matrix(countData)), ctrlgenes, k = 1)
        ctrl_idx <- rownames(counts(counts_norm)) %in% ctrlgenes # for spike-in-derived normalization factors
        counts_norm_mat <- counts(counts_norm)[!ctrl_idx, , drop = FALSE] # removing spike-ins for actual DE testing
        genes <- rownames(counts_norm_mat)
        countData <- countData %>% subset(!row.names(countData) %in% ctrlgenes) # removing spike-ins for standard analysis

        sampleData_norm <- cbind(sampleData, pData(counts_norm))
        design_norm <- model.matrix(as.formula(paste(deparse(des), colnames(pData(counts_norm))[1], sep = " + ")), data = sampleData_norm)
        # colnames(design_norm) <- c(colnames(design),"W_1")

        dge_norm <- DGEList(counts = counts_norm_mat, group = sampleData$condition, samples = samples, genes = genes)

        ## filter low counts; keep original (pre-filter) lib sizes since the spike-in-derived
        ## normalization factors below are only valid relative to them
        keep <- filterByExpr(dge_norm)
        dge_norm <- dge_norm[keep, , keep.lib.sizes = TRUE]

        # relevel to base condition B
        dge_norm$samples$group <- relevel(dge_norm$samples$group, ref = B[[1]])

        ## normalize using the spike-in (control gene) counts rather than TMM on the endogenous
        ## genes: this is what actually puts the spike-in scale into the normalized results, the
        ## W_1 covariate alone only adjusts for unwanted variation, not scale.
        spike_counts <- counts(counts_norm)[ctrl_idx, , drop = FALSE]
        dge_spike <- calcNormFactors(DGEList(counts = spike_counts), method = "TMM")
        dge_norm$samples$norm.factors <- dge_spike$samples$norm.factors

        eff_lib_norm <- dge_norm$samples$lib.size * dge_norm$samples$norm.factors
        scaling_factors <- eff_lib_norm / exp(mean(log(eff_lib_norm)))
        write_scaling_log(scaling_factors, dge_norm$samples$samples, paste("Tables/DE", "EDGER", combi, contrast_name, "table", "scaling.log", sep = "_"))

        ## create file normalized table
        tmm_norm <- as.data.frame(cpm(dge_norm))
        colnames(tmm_norm) <- t(dge_norm$samples$samples)
        tmm_norm$ID <- dge_norm$genes$genes
        tmm_norm <- tmm_norm[c(ncol(tmm_norm), 1:ncol(tmm_norm) - 1)]

        tmm_norm <- add_gene_coordinates(tmm_norm, tmm_norm$ID, gtf_gene, after = "ID")
        write.table(as.data.frame(tmm_norm), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "DataSet", "table", "Normalized_norm.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

        ## create dummy file MDS-plot with and without summarized replicates
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "MDS_norm.png", sep = "_")
        png(out, width=1900, height=1200, res=300)        
        dev.off()

        ## create dummy file BCV-plot - visualizing estimated dispersions
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "BCV_norm.png", sep = "_")
        png(out, width=1900, height=1200, res=300)
        dev.off()

        ## create dummy file quasi-likelihood-dispersion-plot
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "QLDisp_norm.png", sep = "_")
        png(out, width=1900, height=1200, res=300)
        dev.off()
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
    write.table(as.data.frame(tmm), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "DataSet", "table", "Normalized.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

    ## create dummy file MDS-plot with and without summarized replicates
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "MDS.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
    dev.off()

    ## create dummy file BCV-plot - visualizing estimated dispersions
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "BCV.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
    dev.off()

    ## create dummy file quasi-likelihood-dispersion-plot
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "QLDisp.png", sep = "_")
    png(out, width=1900, height=1200, res=300)
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

        ## Testing
        # qlf <- glmQLFTest(fit, contrast=contrast) ## glm quasi-likelihood-F-Test
        #AvsB <- makeContrasts(TreatvsUntreat = paste("condition", A, sep = ""), levels = design)
        ## estimate Dispersion, THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES
        # dge <- estimateDisp(dge, design, robust = TRUE)
        # AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES => no glmQLFTest possible
        # qlf <- glmQLFTest(fit, contrast = AvsB) ## glm quasi-likelihood-F-Test
        qlf <- exactTest(dge, pair = c(B, A), dispersion = bcv^2, prior.count = 2)

        # add comp object to list for image
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

        if (spike != "") { # Same for spike-in normalized
            rm(qlf, tops, res)

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

            ## Testing
            # qlf <- glmQLFTest(fit, contrast=contrast) ## glm quasi-likelihood-F-Test
            #AvsB <- makeContrasts(TreatvsUntreat = paste("condition", A, sep = ""), levels = design)
            # THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES, no glmQLFTest possible
            #qlf <- glmQLFTest(fit, contrast = AvsB) ## glm quasi-likelihood-F-Test
            qlf <- exactTest(dge_norm, pair = c(B, A), dispersion = bcv^2, prior.count = 2)
            # add comp object to list for image
            comparison_objs[[paste0(contrast_name, "_norm")]] <- qlf

            # # Add gene names  (check how gene_id col is named )
            res <- format_edger_results(qlf$table, gtf_gene, center = TRUE)

            # plotVolcano
            pdf(
                file = paste("Figures/DE", "EDGER", combi, contrast_name, "figure_Volcano_norm.pdf", sep = "_"), width = 15, height = 10
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
            write.table(as.data.frame(res), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "results_norm.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

            # create sorted results Tables
            tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
            tops <- format_edger_results(tops$table, gtf_gene, center = TRUE)
            tops <- as.data.frame(apply(tops, 2, as.character))
            write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsLogFCsorted_norm.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

            tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "PValue")
            tops <- format_edger_results(tops$table, gtf_gene, center = TRUE)
            tops <- as.data.frame(apply(tops, 2, as.character))
            write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsPValueSorted_norm.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

            ## plot lFC vs CPM
            out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "MD_norm.png", sep = "_")
            png(out, width=1900, height=1200, res=300)
            print(plotMD(qlf, main = contrast_name))
            abline(h = c(-1, 1), col = "blue")
            dev.off()
        }

        # cleanup
        rm(qlf, res, tops)
        print(paste("cleanup done for ", contrast_name, sep = ""))
    }, error = function(e) {
        message(paste0("Error while processing contrast ", contrast_name, ": ", conditionMessage(e), ". Skipping this contrast."))
    })
}


save.image(file = paste("DE_EDGER", combi, "SESSION.gz", sep = "_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
