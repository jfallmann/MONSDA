suppressPackageStartupMessages({
    require(BiocParallel)
    library(tximport)
    library(GenomicFeatures)
    library(DRIMSeq)
    library(rtracklayer)
    library(readr)
})

options(echo = TRUE)

### ARGS
args <- commandArgs(trailingOnly = TRUE)
argsLen <- length(args)
anname <- args[1]
gtf <- args[2]
outdir <- args[3]
cmp <- args[4]
combi <- args[5]
cores <- as.integer(args[6])
filter <- if (argsLen > 6) args[7] else "min_samps_feature_expr = 1, min_feature_expr = .1, min_samps_gene_expr = 1, min_gene_expr = 1"
print(args)

## FUNCS
libp <- paste0(gsub("/bin/conda", "/envs/monsda", Sys.getenv("CONDA_EXE")), "/share/MONSDA/scripts/lib/_lib.R")
source(libp)

### SCRIPT
# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)

BPPARAM <- MulticoreParam(workers = cores)

# Importing counts
sampleData_all <- read.table(file = gzfile(anname), header = TRUE, row.names = NULL)
sampleData_all$original_sample_id <- sampleData_all$sample_id
sampleData_all$sample_id <- paste("sample", sampleData_all$sample_id, sampleData_all$condition, sep = "_")
sampleData_all$sample_id <- make.names(sampleData_all$sample_id)
sampleData_all$condition <- as.factor(sampleData_all$condition)
sampleData_all$type <- as.factor(sampleData_all$type)
sampleData_all$batch <- as.factor(sampleData_all$batch)

if (file_test("-d", sampleData_all$path[1])) {
    files <- file.path(sampleData_all$path, "quant.sf.gz")
} else {
    files <- sampleData_all$path
}
names(files) <- sampleData_all$sample_id

# Transcript-to-gene mapping
txdb.filename <- file.path(paste(gtf, "sqlite", sep = "."))
txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# get counts
txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM") 
cts_all <- txi$counts[gsub("\\.[0-9]+$", "", rownames(txi$counts)) %in% gsub("\\.[0-9]+$", "", txdf$TXNAME), ]
cts_all <- cts_all[rowSums(cts_all) > 0, ]
rownames(cts_all) <- gsub("\\.[0-9]+$", "", rownames(cts_all))

comparisons <- strsplit(cmp, ",")

if (combi == "none") {
    combi <- ""
}

setwd(outdir)
comparison_objs <- list()

for (contrast in comparisons[[1]]) {
    contrast_name <- strsplit(contrast, ":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast, ":")[[1]][2], "-vs-")
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)
    conds <- c(A, B)

    # Subset sampleData and counts for this comparison
    sampleData <- droplevels(subset(sampleData_all, condition %in% conds))
    cts <- cts_all[, sampleData$sample_id, drop = FALSE]

    # Subset txdf to only transcripts present in cts
    txdf_sub <- txdf[match(rownames(cts), txdf$TXNAME), ]
    counts <- data.frame(gene_id = txdf_sub$GENEID, feature_id = txdf_sub$TXNAME, cts, check.names = FALSE)
    counts <- counts[!is.na(counts$gene_id), ]

    #  QC: Library size plot 
    suppressWarnings(suppressMessages(library(ggplot2)))
    lib_sizes <- data.frame(
        sample = sampleData$sample_id,
        total_counts = colSums(cts)
    )
    p_libsize <- ggplot(lib_sizes, aes(x = sample, y = total_counts)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        xlab("Sample") + ylab("Library Size") +
        ggtitle("Library Size per Sample") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(
        filename = paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "LibrarySize.png", sep = "_"),
        plot = p_libsize, width = 10, height = 6, dpi = 150
    )

    #  QC: Features per gene table 
    features_per_gene <- as.data.frame(table(counts$gene_id))
    colnames(features_per_gene) <- c("gene_id", "num_features")
    write.table(
        features_per_gene,
        file = paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "features_per_gene.tsv", sep = "_"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    # Set up design matrix
    bl <- sapply("batch", paste0, levels(sampleData$batch)[-1])
    tl <- sapply("type", paste0, levels(sampleData$type)[-1])
    if (length(levels(sampleData$type)) > 1) {
        if (length(levels(sampleData$batch)) > 1) {
            design <- model.matrix(~ 0 + type + batch + condition, data = sampleData)
            colnames(design) <- c(tl, bl, levels(sampleData$condition))
        } else {
            design <- model.matrix(~ 0 + type + condition, data = sampleData)
            colnames(design) <- c(tl, levels(sampleData$condition))
        }
    } else {
        if (length(levels(sampleData$batch)) > 1) {
            design <- model.matrix(~ 0 + batch + condition, data = sampleData)
            colnames(design) <- c(bl, levels(sampleData$condition))
        } else {
            design <- model.matrix(~ 0 + condition, data = sampleData)
            colnames(design) <- levels(sampleData$condition)
        }
    }

    # Create dmDSdata object for this comparison
    d <- dmDSdata(counts = counts, samples = sampleData)

    # Filter (use DRIMSeq defaults scaled to group size)
    n <- nrow(sampleData)
    n.small <- max(1, ceiling(n / length(levels(sampleData$condition))))
    default_filter <- paste(
        "min_samps_feature_expr =", n.small,
        ", min_feature_expr = 10",
        ", min_samps_feature_prop =", n.small,
        ", min_feature_prop = 0.1",
        ", min_samps_gene_expr =", n.small,
        ", min_gene_expr = 10"
    )
    filter_string <- ifelse(is.null(filter) || filter == "", default_filter, filter)
    d <- eval(parse(text = paste0("dmFilter(d, ", filter_string, ")")))

    #  QC: Filtering summary table 
    counts_filtered <- counts(d)
    filter_summary <- data.frame(
        n_genes_before = length(unique(txdf_sub$GENEID)),
        n_genes_after = length(unique(counts_filtered$gene_id)),
        n_features_before = nrow(txdf_sub),
        n_features_after = nrow(counts_filtered)
    )
    write.table(
        filter_summary,
        file = paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "filtering_summary.tsv", sep = "_"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    #  QC: Mean-variance plot 
    mean_expr <- rowMeans(cts)
    var_expr <- apply(cts, 1, var)
    meanvar_df <- data.frame(mean = mean_expr, variance = var_expr)
    meanvar_df <- meanvar_df[meanvar_df$mean > 0 & meanvar_df$variance > 0, ]
    p_meanvar <- ggplot(meanvar_df, aes(x = mean, y = variance)) +
        geom_point(alpha = 0.5) +
        scale_x_log10() + scale_y_log10() +
        theme_bw() +
        xlab("Mean Expression") + ylab("Variance") +
        ggtitle("Mean-Variance Plot")
    ggsave(
        filename = paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "MeanVariance.png", sep = "_"),
        plot = p_meanvar, width = 8, height = 6, dpi = 150
    )

    # Set up contrast vector
    # Direction: positive for group A (first side of contrast), negative for group B
    plus <- 1 / length(A)
    minus <- -1 / length(B)
    contrast_vec <- cbind(integer(dim(design)[2]), colnames(design))
    for (i in A) {
        contrast_vec[which(contrast_vec[, 2] == i)] <- plus
    }
    for (i in B) {
        contrast_vec[which(contrast_vec[, 2] == i)] <- minus
    }
    contrast_vec <- as.numeric(contrast_vec[, 1])

    # Run DRIMSeq workflow
    set.seed(1)
    d <- dmPrecision(d, design = design, BPPARAM = BPPARAM)
    d <- dmFit(d, design = design)
    d <- dmTest(d, contrast = contrast_vec)

    #  QC: PCA of transcript proportions (robust to constant/zero columns) 
    proportions_tmp <- tryCatch({ DRIMSeq::proportions(d) }, error = function(e) NULL)
    if (!is.null(proportions_tmp) && ncol(proportions_tmp) > 3) {
        prop_mat <- as.matrix(proportions_tmp[, grep("^sample", colnames(proportions_tmp)), drop = FALSE])
        # Remove columns (samples) that are constant or all zero
        is_const_col <- apply(prop_mat, 2, function(x) length(unique(x)) == 1)
        prop_mat <- prop_mat[, !is_const_col, drop = FALSE]
        # Remove rows (transcripts) that are constant or all zero
        is_const_row <- apply(prop_mat, 1, function(x) length(unique(x)) == 1)
        prop_mat <- prop_mat[!is_const_row, , drop = FALSE]
        if (nrow(prop_mat) > 2 && ncol(prop_mat) > 2) {
            pca <- prcomp(t(prop_mat), scale. = TRUE)
            pca_df <- data.frame(
                PC1 = pca$x[, 1],
                PC2 = pca$x[, 2],
                sample = colnames(prop_mat),
                condition = sampleData$condition[match(colnames(prop_mat), sampleData$sample_id)]
            )
            p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
                geom_point(size = 3) +
                geom_text(vjust = 1.5, hjust = 1.2, size = 2) +
                theme_bw() +
                ggtitle("PCA of Transcript Proportions")
            ggsave(
                filename = paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "PCA_Proportions.png", sep = "_"),
                plot = p_pca, width = 8, height = 6, dpi = 150
            )
        }
    }

    # add comp object to list for image
    comparison_objs[[contrast_name]] <- d

    # calculate LFC from proportions table
    proportions <- DRIMSeq::proportions(d)
    A_lvls <- A[[1]]
    B_lvls <- B[[1]]
    samples_of_group_A <- subset(samples(d), condition %in% A_lvls)$sample_id
    samples_of_group_B <- subset(samples(d), condition %in% B_lvls)$sample_id
    proportions[paste(A[[1]], "mean", sep = "_")] <- rowMeans(proportions[as.vector(samples_of_group_A)], na.rm = TRUE)
    proportions[paste(B[[1]], "mean", sep = "_")] <- rowMeans(proportions[as.vector(samples_of_group_B)], na.rm = TRUE)
    proportions[[paste(A[[1]], "mean", sep = "_")]][is.nan(proportions[[paste(A[[1]], "mean", sep = "_")]])] <- 0
    proportions[[paste(B[[1]], "mean", sep = "_")]][is.nan(proportions[[paste(B[[1]], "mean", sep = "_")]])] <- 0

    # Avoid -Inf when a group mean is zero: add a tiny pseudocount derived from the data scale
    nonzero_props <- as.numeric(unlist(proportions[grep("^sample", colnames(proportions))]))
    nonzero_props <- nonzero_props[nonzero_props > 0]
    eps <- if (length(nonzero_props) > 0) {
        max(1e-10, min(nonzero_props) / 10)
    } else {
        1e-10
    }
    a_mean <- proportions[[paste(A[[1]], "mean", sep = "_")]]
    b_mean <- proportions[[paste(B[[1]], "mean", sep = "_")]]
    proportions[["lfc"]] <- log2((a_mean + eps) / (b_mean + eps))
    props_transcripts <- proportions[c("feature_id", "lfc")]
    props_genes <- aggregate(proportions, list(proportions$gene_id), mean)[c("Group.1", "lfc")]
    colnames(props_genes) <- c("gene_id", "lfc")

    # generate a single p-value per gene and transcript
    res <- DRIMSeq::results(d)
    res.txp <- DRIMSeq::results(d, level = "feature")

    # filter out NA's
    no.na <- function(x) ifelse(is.na(x), 1, x)
    res$pvalue <- no.na(res$pvalue)
    res.txp$pvalue <- no.na(res.txp$pvalue)
    res$adj_pvalue <- no.na(res$adj_pvalue)
    res.txp$adj_pvalue <- no.na(res.txp$adj_pvalue)

    res <- merge(props_genes, res)
    res.txp <- merge(props_transcripts, res.txp)
    res <- res[, c(1, 6, 2, 3, 4, 5)]
    res.txp <- res.txp[, c(1, 7, 2, 3, 4, 5, 6)]

    proportions$Gene <- lapply(proportions$gene_id, function(x) {
        get_gene_name(x, gtf.df)
    })
    res$Gene <- lapply(res$gene_id, function(x) {
        get_gene_name(x, gtf.df)
    })
    res.txp$Gene <- lapply(res.txp$gene_id, function(x) {
        get_gene_name(x, gtf.df)
    })

    res <- res[order(res$adj_pvalue), ]
    res.txp <- res.txp[order(res$adj_pvalue), ]
    proportions <- proportions[order(proportions$lfc), ]

    proportions.print <- as.data.frame(apply(proportions, 2, as.character))
    for (c in 3:(length(colnames(proportions.print)) - 4)) {
        colnames(proportions.print)[c] <- sampleData$original_sample_id[sampleData$sample_id == colnames(proportions.print)[c]]
    }
    res.print <- as.data.frame(apply(res, 2, as.character))
    res.txp.print <- as.data.frame(apply(res.txp, 2, as.character))

    # CREATE RESULTS TABLES
    write.table(as.data.frame(proportions.print), gzfile(paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "proportions.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    write.table(as.data.frame(res.print), gzfile(paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "genes.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    write.table(as.data.frame(res.txp.print), gzfile(paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "transcripts.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    write.table(as.data.frame(genewise_precision(d)), gzfile(paste("Tables/DTU", "DRIMSEQ", combi, contrast_name, "table", "genewise-precision.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

    ## Plot feature per gene histogram
    tryCatch(
        {
            png(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "FeatPerGene.png", sep = "_"), width=1900, height=1200, res=300)
            print(plotData(d))
            dev.off()
        },
        warning = function(w) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "FeatPerGene.png", sep = "_"))
            message("Caught a warning!")
            print(w)
        },
        error = function(e) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "FeatPerGene.png", sep = "_"))
        }
    )
    ## Plot precision
    png(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "Precision.png", sep = "_"), width=1900, height=1200, res=300)
    print(plotPrecision(d))
    dev.off()

    ## Plot gene-level p-values
    png(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "figure", "PValues.png", sep = "_"), width=1900, height=1200, res=300)
    print(plotPValues(d))
    dev.off()

    # plot proportions
    figures <- data.frame()
    sigs <- which(res$adj_pvalue < 0.05)

    limit <- 10
    counter <- 1
    message("create proportions plots")
    for (gene in sigs) {
        if (counter > limit) {
            break
        }
        if (is.na(gene)) {
            next
        }
        name1 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res$Gene[gene], "figure", "plotProportions", "props.png", sep = "_")
        png(name1, width=1900, height=1200, res=300)
        print(plotProportions(d, res$gene_id[gene], group_variable = "condition"))
        dev.off()

        name2 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res$Gene[gene], "figure", "lineplot.png", sep = "_")
        png(name2, width=1900, height=1200, res=300)
        print(plotProportions(d, res$gene_id[gene], group_variable = "condition", plot_type = "lineplot"))
        dev.off()

        name3 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res$Gene[gene], "figure", "ribbonplot.png", sep = "_")
        png(name3, width=1900, height=1200, res=300)
        print(plotProportions(d, res$gene_id[gene], group_variable = "condition", plot_type = "ribbonplot"))
        dev.off()

        figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir, name1, sep = "/")))
        figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir, name2, sep = "/")))
        figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir, name3, sep = "/")))

        counter <- counter + 1
    }
    tryCatch(
        {
            colnames(figures) <- c("geneID", "geneName", "file")
            write.table(figures, paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigGeneFigures.tsv", sep = "_"), sep = "\t", quote = F, row.names = FALSE, col.names = TRUE)
        },
        warning = function(w) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigGeneFigures.tsv", sep = "_"))
            message("Caught a warning!")
            print(w)
        },
        error = function(e) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigGeneFigures.tsv", sep = "_"))
        }
    )

    sigt <- which(res.txp$adj_pvalue < 0.05)
    figures <- data.frame()
    limit <- 10
    counter <- 1
    message("create proportions plots for transcripts")
    for (trans in sigt) {
        if (counter > limit) {
            break
        }
        if (is.na(trans)) {
            next
        }
        name1 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res.txp$Gene[trans], res.txp$feature_id[trans], "figure", "plotProportions_transcript", "props.png", sep = "_")
        png(name1, width=1900, height=1200, res=300)
        print(plotProportions(d, res.txp$gene_id[trans], group_variable = "condition"))
        dev.off()

        name2 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res.txp$Gene[trans], res.txp$feature_id[trans], "figure", "lineplot_transcript.png", sep = "_")
        png(name2, width=1900, height=1200, res=300)
        print(plotProportions(d, res.txp$gene_id[trans], group_variable = "condition", plot_type = "lineplot"))
        dev.off()

        name3 <- paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, res.txp$Gene[trans], res.txp$feature_id[trans], "figure", "ribbonplot_transcript.png", sep = "_")
        png(name3, width=1900, height=1200, res=300)
        print(plotProportions(d, res.txp$gene_id[trans], group_variable = "condition", plot_type = "ribbonplot"))
        dev.off()

        figures <- rbind(figures, c(res.txp$feature_id[trans], res.txp$Gene[trans], paste(outdir, name1, sep = "/")))
        figures <- rbind(figures, c(res.txp$feature_id[trans], res.txp$Gene[trans], paste(outdir, name2, sep = "/")))
        figures <- rbind(figures, c(res.txp$feature_id[trans], res.txp$Gene[trans], paste(outdir, name3, sep = "/")))

        counter <- counter + 1
    }
    tryCatch(
        {
            colnames(figures) <- c("transcriptID", "geneName", "file")
            write.table(figures, paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigTranscriptFigures.tsv", sep = "_"), sep = "\t", quote = F, row.names = FALSE, col.names = TRUE)
        },
        warning = function(w) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigTranscriptFigures.tsv", sep = "_"))
            message("Caught a warning!")
            print(w)
        },
        error = function(e) {
            file.create(paste("Figures/DTU", "DRIMSEQ", combi, contrast_name, "list", "sigTranscriptFigures.tsv", sep = "_"))
        }
    )

    # cleanup
    rm(res, res.txp, proportions)
    print(paste("cleanup done for ", contrast_name, sep = ""))
}

save.image(file = paste("DTU_DRIMSEQ", combi, "SESSION.gz", sep = "_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
