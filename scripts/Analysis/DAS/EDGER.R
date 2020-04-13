
suppressPackageStartupMessages({
  require(edgeR)
  require(dplyr)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)

anname          <- args[1]
countfile       <- args[2]
outdir          <- args[3]
cmp             <- args[4]
availablecores  <- as.integer(args[5])

## for manual use
# wd <-"/home/roberto/Rscripts/DAS3/"
# setwd(wd)
# anname          <- "data/ANNOTATION.gz"
# countfile       <- "data/COUNTS.gz"
# outdir          <- wd
# cmp             <- "AD.allvsCTRL.all:AD_total+AD_3+AD_5_6+AD_white-vs-CTRL_total+CTRL_3+CTRL_5_6+CTRL_white,AD.totalvsAD.single:AD_total-vs-AD_3+AD_5_6+AD_white,CTRL.totalvsCTRL.single:CTRL_total-vs-CTRL_3+CTRL_5_6+CTRL_white,3:AD_3-vs-CTRL_3,5_6:AD_5_6-vs-CTRL_5_6,total:AD_total-vs-CTRL_total,white:AD_white-vs-CTRL_white"
# availablecores  <- as.integer("1")


## Gives Colors for MDS Plot
RainbowColor <- function(groups){
  groupsAsNumbers <- as.numeric(groups)
  spektrum <- rainbow(max(groupsAsNumbers),alpha=1)
  cl <- c()
  for(i in groupsAsNumbers){
    cl <- c(cl,spektrum[i])
  }
  return(cl)
}


### MAIN ###
############

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

message(paste('Will run EdgeR DAS with ',availablecores,' cores',sep=''))

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type")
sampleData <- as.data.frame(sampleData)
groups <- factor(sampleData$condition)
samples <- rownames(sampleData)
types <- factor(sampleData$type)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

## readin counttable
read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
colnames(dcounts) <- c("GeneID", rownames(sampleData))

## create ExonID's
id <- as.character(dcounts[,1])
n <- id
split(n,id) <- lapply(split(n ,id), seq_along )
rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
dcounts <- dcounts[,2:ncol(dcounts)]

## get genes and exon names out
splitted <- strsplit(rownames(dcounts), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)

dge <- DGEList(counts=dcounts, group=groups, samples=samples, genes=genesrle)

## Addd exons to dge
dge$genes$exons <- exons

## filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

## normalize with TMM
dge <- calcNormFactors(dge, method = "TMM", BPPARAM=BPPARAM)

## create file normalized table
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]
write.table(tmm, file=paste(outdir,"All_Conditions_normalized_table.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)

## create file MDS-plot with and without sumarized replicates
out <- paste(outdir,"All_Conditions_MDS.png",sep="")
png(out, width = 400, height = 400)
colors <- RainbowColor(dge$samples$group)
plotMDS(dge, col=colors)
dev.off()
DGEsum <- sumTechReps(dge, ID=groups)
out <- paste(outdir,"All_Conditions_sum_MDS.png", sep="")
png(out, width = 400, height = 400)
colors <- RainbowColor(DGEsum$samples$group)
plotMDS(DGEsum, col=colors)
dev.off()


## Create design-table considering different types (paired, unpaired)
if (length(levels(types))>1){
  design <- model.matrix(~0+groups+types, data=sampleData)
} else {
  design <- model.matrix(~0+groups, data=sampleData)
}
colnames(design) <- levels(groups)

## estimate Dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- paste(outdir,"All_Conditions_BCV.png",sep="")
png(out, width = 400, height = 400)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust=TRUE)

## create file quasi-likelihood-dispersion-plot
out <- paste(outdir,"All_Conditions_QLDisp.png",sep="")
png(out, width = 400, height = 400)
plotQLDisp(fit)
dev.off()

## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

  contrast_name <- strsplit(contrast,":")[[1]][1]
  contrast_groups <- strsplit(strsplit(contrast,":")[[1]][2], "-vs-")

  message(paste("Comparing ",contrast_name, sep=""))

  # determine contrast
  A <- strsplit(contrast_groups[[1]][1], "+")
  B <- strsplit(contrast_groups[[1]][2], "+")
  minus <- 1/length(A[[1]])*(-1)
  plus <- 1/length(B[[1]])
  contrast <- cbind(integer(dim(design)[2]), colnames(design))
  for(i in A[[1]]){
    contrast[which(contrast[,2]==i)]<- minus
  }
  for(i in B[[1]]){
    contrast[which(contrast[,2]==i)]<- plus
  }
  contrast <- as.numeric(contrast[,1])

  # create files topSpliced by gene, simes and exon method
  sp <- diffSpliceDGE(fit, contrast=contrast, geneid="genes", exonid="exons")
  tops <- topSpliceDGE(sp, test="gene", n=length(fit$counts))
  write.table(tops, file=paste(outdir,contrast_name,"_diffSplice_geneTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)
  tops <- topSpliceDGE(sp, test="simes", n=length(fit$counts))
  write.table(tops, file=paste(outdir,contrast_name,"_diffSplice_simesTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)
  tops <- topSpliceDGE(sp, test="exon", n=length(fit$counts))
  write.table(tops, file=paste(outdir,contrast_name,"_diffSplice_exonTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)

  # create files diffSplicePlots
  tops <- topSpliceDGE(sp, test="simes", n=10)
  for(i in 1:10){
    geneID <- tops$genes[i]
    out <- paste(outdir,contrast_name,"_topSplice_simes_",i,".png",sep="")
    png(out, width = 800, height = 400)
    plotSpliceDGE(sp, geneid=geneID, genecol="genes")
    dsdev.off()
  }
}

save.image(file = paste(outdir,"EDGER_DAS_SESSION.gz",sep=""), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)