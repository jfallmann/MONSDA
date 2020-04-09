
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
# wd <-"/home/roberto/Rscripts/DAS2/"
# setwd(wd)
# anname          <- "data/ANNOTATION.gz"
# countfile       <- "data/COUNTS.gz"
# outdir          <- wd
# cmp             <- "AD_3-vs-CTRL_3,AD_5_6-vs-CTRL_5_6,AD_total-vs-CTRL_total,AD_white-vs-CTRL_white"
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

DGEListFromFeatureCounts <- function (countfile, sampleData, comp){

  ## readin counttable
  read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
  colnames(dcounts) <- c("GeneID", rownames(sampleData))

  ## get subset of annotation and counts by comparisonpair
  if (!is.null(comp)) {
    sampleData <- sampleData[which(sampleData$condition==comp[[1]][1] | sampleData$condition==comp[[1]][2]),]
    dcounts <- dcounts[c(1,grep(comp[[1]][1],colnames(dcounts)), grep(comp[[1]][2],colnames(dcounts)))]
  }

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

  dge <- DGEList(counts=dcounts, group=factor(sampleData$condition), samples=rownames(sampleData), genes=genesrle)

  ## Addd exons to dge
  dge$genes$exons <- exons

  ## return DGEList and possibly subsettet sampleData
  return(list("dge"=dge, "SD"=sampleData))
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

## Combinations of conditions
comparison <- strsplit(cmp, ",")

## create file colored MDS-plot with all conditions
res = DGEListFromFeatureCounts(countfile, sampleData, NULL)
DGE <- res$dge
conditions <- res$SD$condition
keep <- filterByExpr(DGE)
DGE <- DGE[keep, , keep.lib.sizes=FALSE]
DGE <- calcNormFactors(DGE, method = "TMM", BPPARAM=BPPARAM)
DGE <- sumTechReps(DGE, ID=conditions)
out <- "All_Conditions_MDS.png"
png(out, width = 400, height = 400)
colors <- RainbowColor(DGE$samples$group)
plotMDS(DGE, col=colors)
dev.off()


## read in count table and normalize, pairwise
for(pair in comparison[[1]]){

  # pair <- comparison[[1]][3]

  message(paste("Comparing ",pair,sep=""))
  comp<- strsplit(pair,"-vs-")

  res <- DGEListFromFeatureCounts(countfile, sampleData, comp)
  DGE <- res$dge
  conditions <- res$SD$condition

  keep <- filterByExpr(DGE)
  table(keep)
  DGE <- DGE[keep, , keep.lib.sizes=FALSE]
  DGE <- calcNormFactors(DGE, method = "TMM", BPPARAM=BPPARAM)

  # check for different types (paired, unpaired)
  if (length(levels(res$SD$type))>1){
    design <- model.matrix(~factor(conditions)+res$SD$type)
  } else {
    design <- model.matrix(~factor(conditions))
  }

  # create file normalized table
  tmm <- as.data.frame(cpm(DGE))
  colnames(tmm) <- t(DGE$samples$samples)
  tmm$ID <- DGE$genes$genes
  tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]
  write.table(tmm, file=paste(pair,"_normalized_table.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)

  # create file colored MDS-plot
  out <- paste(pair,"_MDS.png",sep="")
  png(out, width = 400, height = 400)
  colors <- RainbowColor(DGE$samples$group)
  plotMDS(DGE, col=colors)
  dev.off()

  # create file BCV-plot
  DGE <- estimateDisp(DGE, design, robust=TRUE)
  out <- paste(pair,"_BCV.png",sep="")
  png(out, width = 400, height = 400)
  plotBCV(DGE)
  dev.off()

  # create file QLDisp-plot
  fit <- glmQLFit(DGE, design, robust=TRUE)
  out <- paste(pair,"_QLDisp.png",sep="")
  png(out, width = 400, height = 400)
  plotQLDisp(fit)
  dev.off()

  # create files topSpliced by gene, simes and exon method
  sp <- diffSpliceDGE(fit, geneid="genes", exonid="exons")
  tops <- topSpliceDGE(sp, test="gene", n=length(fit$counts))
  write.table(tops, file=paste(pair,"_diffSplice_geneTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)
  tops <- topSpliceDGE(sp, test="simes", n=length(fit$counts))
  write.table(tops, file=paste(pair,"_diffSplice_simesTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)
  tops <- topSpliceDGE(sp, test="exon", n=length(fit$counts))
  write.table(tops, file=paste(pair,"_diffSplice_exonTest.tsv",sep=""), sep="\t", quote=F, row.names=FALSE)

  # create files diffSplicePlots
  tops <- topSpliceDGE(sp, test="simes", n=10)
  for(i in 1:10){
    geneID <- tops$genes[i]
    out <- paste(pair,"_topSplice_simes_",i,".png",sep="")
    png(out, width = 800, height = 400)
    plotSpliceDGE(sp, geneid=geneID, genecol="genes")
    dev.off()
  }
}

save.image(file = "EDGER_DAS_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
