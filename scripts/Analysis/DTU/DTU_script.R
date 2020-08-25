suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(DRIMSeq)
  library(stageR)
  library(DEXSeq)
  library(DESeq2)
  library(edgeR)
})


# args <- commandArgs(trailingOnly = TRUE)
# 
# anname          <- args[1]
# countfile       <- args[2]
# outdir          <- args[3]
# cmp             <- args[4]
# availablecores  <- as.integer(args[5])

wd <- "/scr/k70san2/robin/dtu_testing"     

# anname          <- 
gtf             <- file.path(wd,"GENOMES/gencode.v35.annotation.gtf.gz")
# countfile       <- args[2]
# outdir          <- args[3]
# cmp             <- args[4]
# availablecores  <- 15

sample_id <- c("WT_d01","WT_d02","WT_d03","WT_d04","WT_d06","WT_d07","WT_d08","KO_d02","KO_d03","KO_d04","KO_d06","KO_d07","KO_d08")
condition <- c("WT","WT","WT","WT","WT","WT","WT","KO","KO","KO","KO","KO","KO")
samps <- data.frame(sample_id, condition)
samps$condition <- factor(samps$condition)

files <- file.path(wd,"DTU/quant", samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
files[1] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d01/quant.sf" 
files[2] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d02/quant.sf" 
files[3] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d03/quant.sf" 
files[4] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d04/quant.sf" 
files[5] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d06/quant.sf" 
files[6] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d07/quant.sf" 
files[7] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/WT/d08/quant.sf" 
files[8] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d02/quant.sf" 
files[9] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d03/quant.sf" 
files[10] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d04/quant.sf" 
files[11] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d06/quant.sf" 
files[12] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d07/quant.sf" 
files[13] <- "/scr/k70san2/robin/dtu_testing/DTU/quant/KO/d08/quant.sf" 

txi <- tximport(files, type='salmon', txOut=TRUE, countsFromAbundance = "scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

gtf <- gtf
txdb.filename <- file.path(wd,"GENOMES/gencode.v35.annotation.sqlite")

txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

# DRIMSEQ
counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)
d <- dmDSdata(counts=counts, samples=samps)

n <- 13
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
table(table(counts(d)$gene_id))
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)

# 1 estimate the precision, 
# 2 fit regression coefficients and perform null hypothesis testing on the coefficient of interest, 
# 3 test the coefficient associated with the difference between condition 2 and condition 1
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="conditionWT")
})

res <- DRIMSeq::results(d)
res.txp <- DRIMSeq::results(d, level="feature")

# filter out NA's
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

# plot the estimated proportions for one of the significant genes
idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
plotProportions(d, res$gene_id[idx], "condition")

nrow(subset(res, adj_pvalue < 0.05))



## two-stage testing procedure
#vector of p-values for the screening stage
pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)

# column matrix of the confirmation p-values
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)

# two column data.frame with the transcript and gene identifiers
tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# stageR analysis following DRIMSEQ
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

head(drim.padj)
nrow(subset(drim.padj, transcript < 0.05))


# Post-hoc filtering on the standard deviation in proportions
res.txp.filt <- DRIMSeq::results(d, level="feature")
smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}
filt <- smallProportionSD(d)
res.txp.filt$pvalue[filt] <- 1 
res.txp.filt$adj_pvalue[filt] <- 1




# DEXSEQ analysis
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)

# stageR analysis following DEXSEQ
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
})

head(dex.padj)





# DESEQ2 analysis
txi.g <- tximport(files, type="salmon", tx2gene=txdf[,2:1])

dds <- DESeqDataSetFromTximport(txi.g, samps, ~condition)
dds <- DESeq(dds)
dres <- DESeq2::results(dds)
all(dxr.g$gene %in% rownames(dres))
dres <- dres[dxr.g$gene,]
col <- rep(8, nrow(dres))
col[rownames(dres) %in% dge.genes] <- 1
col[rownames(dres) %in% dte.genes] <- 2
col[rownames(dres) %in% dtu.genes] <- 3

bigpar()
# here cap the smallest DESeq2 adj p-value
cap.padj <- pmin(-log10(dres$padj), 100)
# this vector only used for plotting
jitter.padj <- -log10(dxr.g$qval + 1e-20)
jp.idx <- jitter.padj == 20
jitter.padj[jp.idx] <- rnorm(sum(jp.idx),20,.25)
plot(cap.padj, jitter.padj, col=col,
     xlab="Gene expression",
     ylab="Transcript usage")
legend("topright",
       c("DGE","DTE","DTU","null"),
       col=c(1:3,8), pch=20, bty="n")


# EDGER analysis
cts.g <- txi.g$counts
normMat <- txi.g$length
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts.g/normMat)) + log(colSums(cts.g/normMat))
y <- DGEList(cts.g)
y <- scaleOffset(y, t(t(log(normMat)) + o))
keep <- filterByExpr(y)
y <- y[keep,]

y <- estimateDisp(y, design_full)
fit <- glmFit(y, design_full)
lrt <- glmLRT(fit)
tt <- topTags(lrt, n=nrow(y), sort="none")[[1]]

common <- intersect(res$gene_id, rownames(tt))
tt <- tt[common,]
res.sub <- res[match(common, res$gene_id),]

plot(-log10(tt$FDR), -log10(res.sub$adj_pvalue), col=col,
     xlab="Gene expression",
     ylab="Transcript usage")
legend("topright",
       c("DGE","DTE","DTU","null"),
       col=c(1:3,8), pch=20, bty="n")
