require(ggplot2)
require(ggpubr)
require(magrittr)
require(edgeR)
require(devtools)
require(jsonlite)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
smp <- args[2]
grp <- args[3]
exp <- args[4]
cond <- args[5]

###READ IN
data <- read.table(file,h=F,sep="\t")
samples <- as.vector(t(read.table(smp,h=F,sep="\t"))[,1])
groups <- as.vector(t(read.table(grp,h=F,sep="\t"))[,1])
genes <- data$V1

###GET RID OF ID COLUMN AND CONVERT NA TO 0
counts <- data[,-c(1)]
counts[is.na(counts)] <- 0

summary(counts)

###CREATE DE-GENE-LIST

dge <- DGEList(counts=counts, group=groups, samples=samples, genes=genes)

###REMOVE LOW EXPRESSION

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge$samples
dge$genes

###Normalize by TMM

dge <- calcNormFactors(dge)
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]
write.table(tmm, file=paste(file,"_normalized.tsv",sep=''), sep="\t", quote=F, row.names=FALSE)

#dge
#colnames(dge)

###plot MDS
out <- paste(file,"_MDS",".png",sep="")
png(out, width = 350, height = 350)
plotMDS(dge)
dev.off()

###Design table
Experiments <- as.vector(t(read.table(exp,h=F,sep="\t"))[,1])
Conditions <- as.vector(t(read.table(cond,h=F,sep="\t"))[,1])

design <- model.matrix(~0+Conditions)
rownames(design) <- colnames(dge)

#design
###Estimate dispersion

dge <- estimateDisp(dge, design, robust=TRUE)

###plot Dispersion

out <- paste(file,"_BCV",".png",sep="")
png(out, width = 350, height = 350)
plotBCV(dge)
dev.off()

###Estimate DE
####fit glms
fit <- glmFit(dge, design)
####likelihood ratiotest
lrt <- glmLRT(fit)
####top DE
topTags(lrt)

###DE distribution
summary(decideTests(lrt))

###plot lFC vs CPM

out <- paste(file,"_MD",".png",sep="")
png(out, width = 350, height = 350)
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

###Save session info
session_info() %>%
    write_json(paste(file,"_session_info.json",sep=''), force = TRUE)
