suppressPackageStartupMessages({
    require(utils)
    require(BiocParallel)
    require(DESeq2)
    require(rtracklayer)
    require(RUVSeq)
    require(dplyr)
    require(GenomeInfoDb)
    require(apeglm)
})

setwd("/scr/k70san3/berni/Kidneys/fall/")

anname          <- "DE/fastqc-trimgalore-hisat2/deseq2_DE/Tables.bak/fastqc-trimgalore-hisat2_ANNOTATION.gz"
countfile       <- "DE/fastqc-trimgalore-hisat2/deseq2_DE/Tables.bak/fastqc-trimgalore-hisat2_COUNTS.gz"
gtf             <- "GENOMES/mm10/gencode.vM25.annotation.gtf.gz"
outdir          <- "/scr/k70san2/robin/fall/deseq_DE"
cmp             <-  "T2DMvsCTRL:T2DM-vs-CTRL,T2DMvsT2DMGC:T2DM-vs-T2DMGC,T2DMGCvsCTRL:T2DMGC-vs-CTRL,T1DMvsCTRL:T1DM-vs-CTRL,T1DMvsT1DMGC:T1DM-vs-T1DMGC,T1DMGCvsCTRL:T1DMGC-vs-CTRL"
combi           <- "fastqc-trimgalore-hisat2"
availablecores  <- 16

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

### SCRIPT
## Annotation
sampleData <- as.data.frame(read.table(gzfile(anname), row.names=1))
colnames(sampleData) <- c("condition", "type", "batch")
sampleData$batch <- as.factor(sampleData$batch)
sampleData$type <- as.factor(sampleData$type)
sampleData$condition <- as.factor(sampleData$condition)

## readin counttable
countData <- as.matrix(read.table(gzfile(countfile), header=T, row.names=1))

## Combinations of conditions
comparison <- strsplit(cmp, ",")

## Create design-table considering different types (paired, unpaired) and batches
design <- ~ batch + condition

#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design= design)

#filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#run for each pair of conditions
dds <- DESeq(dds, parallel=TRUE, BPPARAM=BPPARAM)  #, betaPrior=TRUE)
resultsNames(dds)

A <- "T1DM"
B <- "CTRL"

res <- results(dds, contrast=c('condition', A, B), parallel=TRUE, BPPARAM=BPPARAM)
res
shrink <- lfcShrink(dds=dds, coef=paste("condition",A,"vs",B,sep="_"), res=res, type='apeglm')
shrink

setwd(outdir)

png("MA_res.png")
plotMA(res)
dev.off()

png("MA_shrink.png")
plotMA(shrink)
dev.off()
