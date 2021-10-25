#https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

suppressPackageStartupMessages({
    require(utils)
    require(BiocParallel)
    require(DESeq2)
    require(ggrepel)
    require(rtracklayer)
    require(RUVSeq)
    require(dplyr)
    require(GenomeInfoDb)
    require(apeglm)
})

options(echo=TRUE)

## ARGS
args            <- commandArgs(trailingOnly = TRUE)
argsLen         <- length(args);
anname          <- args[1]
countfile       <- args[2]
gtf             <- args[3]
outdir          <- args[4]
cmp             <- args[5]
combi           <- args[6]
availablecores  <- as.integer(args[7])
spike           <- if (argsLen > 7) args[8] else ''

print(args)

### FUNCS
get_gene_name <- function(id, df){
    name_list <- df$gene_name[df['gene_id'] == id]
    if(length(unique(name_list)) == 1){
        return(name_list[1])
    }else{
        message(paste("WARNING: ambigous gene id: ", id))
        return (paste(unique(name_list), sep="|"))
    }
}

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

### SCRIPT
## Annotation
sampleData_all <- as.data.frame(read.table(gzfile(anname), row.names=1))
colnames(sampleData_all) <- c("condition", "type", "batch")
sampleData_all$batch <- as.factor(sampleData_all$batch)
sampleData_all$type <- as.factor(sampleData_all$type)
sampleData_all$condition <- as.factor(sampleData_all$condition)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)
gtf_gene <- droplevels(subset(gtf.df, type == "gene"))

## Combinations of conditions
comparison <- strsplit(cmp, ",")

## check combi
if (combi == "none"){
    combi <- ''
}

## readin counttable
countData_all <- as.matrix(read.table(gzfile(countfile), header=T, row.names=1))

#Check if names are consistent
if (!all(rownames(sampleData_all) %in% colnames(countData_all))){
    stop("Count file does not correspond to the annotation file")
}

comparison_objs <- list()

WD <- getwd()
setwd(outdir)

for(contrast in comparison[[1]]){

    contrast_name <- strsplit(contrast, ":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast, ":")[[1]][2], "-vs-")
    print(paste("Comparing ", contrast_name, sep=""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names=FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names=FALSE)

    #subset Datasets for pairwise comparison
    countData <- cbind(countData_all[ , grepl( paste(B, '_', sep='') , colnames( countData_all ) )], countData_all[ ,grepl(paste(A, '_', sep='') , colnames( countData_all )) ])
    sampleData <- droplevels(rbind(subset(sampleData_all, B == condition), subset(sampleData_all, A == condition)))

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == condition)$type)) > 1 | length(unique(subset(sampleData, B == condition)$type)) > 1){
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1){
            design <- ~ type + batch + condition
        } else{
            design <- ~ type + condition
        }
    } else{
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1){
            design <- ~ batch + condition
        } else{
            design <- ~ condition
        }
    }
    print(design)

    # Normalize by spike in if available
    if (spike != ''){
        print("Spike-in used, data will be normalized to spike in separately")
        spike = strsplit(spike, "=")[[1]][2]
        setwd(WD)
        ctrlgenes <- readLines(spike)
        setwd(outdir)

        counts_norm <- RUVg(newSeqExpressionSet(as.matrix(countData)), ctrlgenes, k=1)
        countData <- countData %>% subset(!row.names(countData) %in% ctrlgenes)  # removing spike-ins for standard analysis

        sampleData_norm <- cbind(sampleData, pData(counts_norm))
        design_norm <- as.formula(paste(deparse(design), colnames(pData(counts_norm))[1], sep=" + "))

        dds_norm <- DESeqDataSetFromMatrix(countData = counts(counts_norm), colData = sampleData_norm, design= design)

        #filter low counts
        keep_norm <- rowSums(counts(dds_norm)) >= 10
        dds_norm <- dds_norm[keep_norm,]

        #drop unused samples
        dds_norm$condition <- droplevels(dds_norm$condition)

        #relevel to base condition B
        dds_norm$condition <- relevel(dds_norm$condition, ref = B[[1]])

        dds_norm <- DESeq(dds_norm, parallel=TRUE, BPPARAM=BPPARAM, betaPrior=FALSE)
        rld_norm <- rlogTransformation(dds_norm, blind=FALSE)
        vsd_norm <-varianceStabilizingTransformation(dds_norm, blind=FALSE)

        png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "PCA_norm.png", sep="_"))
        DESeq2::plotPCA(rld_norm, intgroup=c('condition')) + geom_text_repel(aes(label = name), arrow = arrow(length = unit(0.02, "npc")), box.padding = .5)  # requires ggrepel
        # DESeq2::plotPCA(rld_norm, intgroup=c('condition')) + geom_text(aes(label = name), position = position_nudge(y = 2))
        dev.off()

        #We also write the normalized counts to file
        write.table(as.data.frame(assay(rld_norm)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "rld_norm.tsv.gz", sep="_")), sep="\t", col.names=NA)
        write.table(as.data.frame(assay(vsd_norm)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "vsd_norm.tsv.gz", sep="_")), sep="\t", col.names=NA)
    }

    #Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = design)

    #filter low counts
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]

    #drop unused samples
    dds$condition <- droplevels(dds$condition)

    #relevel to base condition B
    dds$condition <- relevel(dds$condition, ref = B[[1]])

    #run for each pair of conditions
    dds <- DESeq(dds, parallel=TRUE, BPPARAM=BPPARAM, betaPrior=FALSE)
    print(resultsNames(dds))

    #Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

    rld<- rlogTransformation(dds, blind=FALSE)
    vsd<-varianceStabilizingTransformation(dds, blind=FALSE)

    png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "PCA.png", sep="_"))
    DESeq2::plotPCA(rld, intgroup=c('condition')) + geom_text_repel(aes(label = name), arrow = arrow(length = unit(0.02, "npc")), box.padding = .5)  # requires ggrepel
    # DESeq2::plotPCA(rld, intgroup=c('condition')) + geom_text(aes(label = name), position = position_nudge(y = 2))
    dev.off()

    #We also write the normalized counts to file
    write.table(as.data.frame(assay(rld)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "rld.tsv.gz", sep="_")), sep="\t", col.names=NA)
    write.table(as.data.frame(assay(vsd)), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "vsd.tsv.gz", sep="_")), sep="\t", col.names=NA)

    tryCatch({

        # initialize empty objects
        res=""
        resOrdered=""
        res <- results(dds, contrast=c('condition', A, B), parallel=TRUE, BPPARAM=BPPARAM)
        res_shrink <- lfcShrink(dds=dds, coef=paste("condition", A, "vs", B, sep="_"), res=res, type='apeglm')

        # add comp object to list for image
        comparison_objs[[contrast_name]] <- res

        # sort and output
        resOrdered <- res_shrink[order(res_shrink$log2FoldChange),]

        # # Add gene names  (check how gene_id col is named )
        resOrdered$Gene  <- lapply(rownames(resOrdered) , function(x){get_gene_name(x, gtf_gene)})
        resOrdered$Gene_ID <- rownames(resOrdered)
        resOrdered <- resOrdered[, c(7,6,1,2,3,4,5)]
        resOrdered <- as.data.frame(apply(resOrdered, 2, as.character))

        # write the table to a tsv file
        write.table(as.data.frame(resOrdered), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results.tsv.gz", sep="_")), sep="\t", row.names=FALSE, quote=F)

        # sort and output
        res <- res[order(res$log2FoldChange),]

        res$Gene  <- lapply(rownames(res) , function(x){get_gene_name(x, gtf_gene)})
        res$Gene_ID <- rownames(res)
        res <- res[, c(7,6,1,2,3,4,5)]
        res <- as.data.frame(apply(res, 2, as.character))

        write.table(as.data.frame(res), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results_noshrink.tsv.gz", sep="_")), sep="\t", row.names=FALSE, quote=F)

        # plotMA
        png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "MA.png", sep="_"))
        DESeq2::plotMA(res_shrink)
        dev.off()

        if (spike != ''){  # DE run for spike-in normalized data

            # initialize empty objects
            res=""
            resOrdered=""
            res <- results(dds_norm, contrast=c('condition', A, B), parallel=TRUE, BPPARAM=BPPARAM)
            res_shrink <- lfcShrink(dds=dds_norm, coef=paste("condition", A, "vs",B,sep="_"), res=res, type='apeglm')

            # add comp object to list for image
            listname <- paste(contrast_name, "_norm",sep="")
            comparison_objs[[listname]] <- res

            # sort and output
            resOrdered <- res_shrink[order(res_shrink$log2FoldChange),]

            # # Add gene names  (check how gene_id col is named )
            resOrdered$Gene  <- lapply(rownames(resOrdered), function(x){get_gene_name(x, gtf_gene)})
            resOrdered$Gene_ID <- rownames(resOrdered)
            resOrdered <- resOrdered[,c(7,6,1,2,3,4,5)]
            resOrdered <- as.data.frame(apply(resOrdered, 2, as.character))

            # write the table to a tsv file
            write.table(as.data.frame(resOrdered), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results_norm.tsv.gz", sep="_")), sep="\t", row.names=FALSE, quote=F)


            # sort and output
            res <- res[order(res$log2FoldChange),]

            res$Gene  <- lapply(rownames(res) , function(x){get_gene_name(x, gtf_gene)})
            res$Gene_ID <- rownames(res)
            res <- res[, c(7,6,1,2,3,4,5)]
            res <- as.data.frame(apply(res, 2, as.character))

            write.table(as.data.frame(res), gzfile(paste("Tables/DE", "DESEQ2", combi, contrast_name, "table", "results_norm_noshrink.tsv.gz", sep="_")), sep="\t", row.names=FALSE, quote=F)

            # plotMA
            png(paste("Figures/DE", "DESEQ2", combi, contrast_name, "figure", "MA_norm.png", sep="_"))
            DESeq2::plotMA(res)
            dev.off()
        }

        # cleanup
        rm(res, resOrdered)
        print(paste('cleanup done for ', contrast_name, sep=''))
    })
}

#### Now plot and print over-all comparisons

## Create design-table considering different types (paired, unpaired) and batches
design <- ~0 + condition
print(design)

#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData_all, colData = sampleData_all, design = design)

#filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#run for each pair of conditions
dds <- DESeq(dds, parallel=TRUE, BPPARAM=BPPARAM, betaPrior=FALSE)

rld<- rlogTransformation(dds, blind=FALSE)
vsd<-varianceStabilizingTransformation(dds, blind=FALSE)

png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "PCA.png", sep="_"))
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

#We also write the normalized counts to file
write.table(as.data.frame(assay(rld)), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "rld.tsv.gz", sep="_")), sep="\t", col.names=NA)
write.table(as.data.frame(assay(vsd)), gzfile(paste("Tables/DE", "DESEQ2", combi, "DataSet", "table", "vsd.tsv.gz", sep="_")), sep="\t", col.names=NA)


# Here we choose blind so that the initial conditions setting does not influence the outcome, ie we want to see if the conditions cluster based purely on the individual datasets, in an unbiased way. According to the documentation, the rlogTransformation method that converts counts to log2 values is apparently better than the old varienceStabilisation method when the data size factors vary by large amounts.

par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'black')
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "VST-and-log2.png", sep="_"))
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type='l', lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

##############################
library('RColorBrewer')
library('gplots')
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap1.png", sep="_"), width=800, height=750)
heatmap.2(counts(dds, normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(12,8))
dev.off()
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap2.png", sep="_"), width=800, height=750)
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(12, 8))
dev.off()
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap3.png", sep="_"), width=800, height=750)
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(12, 8))
dev.off()

#The above shows heatmaps for 30 most highly expressed genes (not necessarily the biggest fold change). The data is of raw counts (left), regularized log transformation (center) and from variance stabilizing transformation (right) and you can clearly see the effect of the transformation has by shrinking the variance so that we don’t get the squish effect shown in the left hand graph.
##############################
#Now we calculate sample to sample distances so we can make a dendrogram to look at the clustering of samples.
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), condition)
#updated in latest vignette (See comment by Michael Love)
#this line was incorrect
#heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(16, 16))
#From the Apr 2015 vignette
hc <- hclust(distsRL)
png(paste("Figures/DE", "DESEQ2", combi, "DataSet", "figure", "heatmap-samplebysample.png", sep="_"), width=800, height=800)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.off()

##############################

save.image(file = paste("DE", "DESEQ2", combi, "SESSION.gz", sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
