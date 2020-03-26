
suppressPackageStartupMessages({
  require(edgeR)
  require(dplyr)
  require(GenomicRanges)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)

anname          <- args[1]
countfile       <- args[2]
flatanno        <- args[3]
outdir          <- args[4]
cmp             <- args[5]
availablecores  <- as.integer(args[6])

### for manual use
#wd <-"/path/to/directory/"
#setwd(wd)
#anname          <- "RUN_DEU_Analysis.anno.gz"
#countfile       <- "RUN_DEU_Analysis.tbl.gz"
#flatanno        <- "Homo_Anno_dexseq.gtf.gz"
#outdir          <- wd
#cmp             <- "AD_3-vs-CTRL_3,AD_5_6-vs-CTRL_5_6,AD_total-vs-CTRL_total,AD_white-vs-CTRL_white"
#availablecores  <- as.integer("1")


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

## Read Fcount output and convert to DGEList
DGEListFromFeatureCounts <- function (countfile, sampleData, flattenedfile=NULL, comp=NULL){

  ##  Take a fcount file and convert it to dcounts for dexseq
  message("Read Fcount output and convert to DGEList")

  if (!is.null(comp)) {
    samples   <- c(as.character(rownames(sampleData[which(sampleData$condition==comp[[1]][1]),])),as.character(rownames(sampleData[which(sampleData$condition==comp[[1]][2]),])))
    conditions   <- c(as.character(sampleData[which(sampleData$condition==comp[[1]][1]),]$condition),as.character(sampleData[which(sampleData$condition==comp[[1]][2]),]$condition))
    read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
    colnames(dcounts) <- c("GeneID", rownames(sampleData) )
    dcounts   <- dcounts[c(1,grep(comp[[1]][1],colnames(dcounts)), grep(comp[[1]][2],colnames(dcounts)))]
  }else{
    samples <-rownames(sampleData)
    conditions <- sampleData$condition
    read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
    colnames(dcounts) <- c("GeneID", rownames(sampleData) )
  }

  id <- as.character(dcounts[,1])
  n <- id
  split(n,id) <- lapply(split(n ,id), seq_along )
  rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
  dcounts <- dcounts[,2:ncol(dcounts)]

  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name
  ##filter low counts
  keep <- rowSums(dcounts) >= 10
  dcounts <- dcounts[keep,]

  ## get genes and exon names out
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply(splitted, "[[", 1)

  ## parse the flattened file
  if (!is.null(flattenedfile)) {
    aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
    colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
    aggregates$strand <- gsub("\\.", "*", aggregates$strand)
    aggregates <- aggregates[which(aggregates$class == "exonic_part"),]  # exonic_part
    aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
    aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
    # trim the gene_ids to 255 chars in order to match with featurecounts
    longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
    warning(paste0(longIDs, " aggregate geneIDs were found truncated in featureCounts output"), call. = FALSE)
    aggregates$gene_id <- substr(aggregates$gene_id,1,255)

    transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
    transcripts <- strsplit(transcripts, "\\+")
    exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr) # exonic_part_number
    exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
    names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")

    names(transcripts) <- names(exoninfo)
    if (!all(rownames(dcounts) %in% names(exoninfo))) {
      stop("Count files do not correspond to the flattened annotation file")
    }
    matching <- match(rownames(dcounts), names(exoninfo))
    stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
    stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))

    dge <- DGEList(counts=dcounts, group=conditions, samples=samples, genes=exoninfo[matching])
    return(dge)
  }else{
    dge <- DGEList(counts=dcounts, group=conditions, samples=samples, genes=genesrle)
    return(dge)
  }
}

### MAIN ###
############

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type")
sampleData <- as.data.frame(sampleData)

## Combinations of conditions
comparison <- strsplit(cmp, ",")

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)
message(paste('Will run EdgeR with ',availablecores,' cores',sep=''))

## create file colored MDS-plot with all conditions
DGE = DGEListFromFeatureCounts(countfile, sampleData, flattenedfile = flatanno)
keep <- filterByExpr(DGE)
DGE <- DGE[keep, , keep.lib.sizes=FALSE]
DGE <- calcNormFactors(DGE, method = "TMM", BPPARAM=BPPARAM)
tmm <- as.data.frame(cpm(DGE))
colnames(tmm) <- t(DGE$samples$samples)
tmm$ID <- DGE$genes$genes
tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]
out <- "All_Conditions_MDS.png"
png(out, width = 700, height = 700)
colors <- RainbowColor(DGE$samples$group)
plotMDS(DGE, col=colors)
dev.off()

## read in count table and normalizez, pairwise
for(pair in comparison[[1]]){
  message(paste("Comparing ",pair,sep=""))
  comp<- strsplit(pair[[1]],"-vs-")

  DGE = DGEListFromFeatureCounts(countfile, sampleData, flattenedfile=flatanno, comp=comp)

  keep <- filterByExpr(DGE)
  DGE <- DGE[keep, , keep.lib.sizes=FALSE]

  DGE <- calcNormFactors(DGE, method = "TMM", BPPARAM=BPPARAM)

  tmm <- as.data.frame(cpm(DGE))
  colnames(tmm) <- t(DGE$samples$samples)
  tmm$ID <- DGE$genes$genes
  tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]

  # create file normalized table
  write.table(tmm, file=paste("normalized_table_",pair,".tsv",sep=""), sep="\t", quote=F, row.names=FALSE)

  # create file colored MDS-plot
  out <- paste(pair[[1]],"_MDS.png",sep="")
  png(out, width = 400, height = 400)
  colors <- RainbowColor(DGE$samples$group)
  plotMDS(DGE, col=colors)
  dev.off()

  #
  design <- model.matrix(~0+group, data=DGE$samples)
  colnames(design) <- levels(DGE$samples$group)
  DGE <- estimateDisp(DGE, design, robust=TRUE, BPPARAM=BPPARAM)

  # create file BCV-plot
  out <- paste(pair[[1]],"_BCV.png",sep="")
  png(out, width = 400, height = 400)
  plotBCV(DGE)
  dev.off()

  # create file QLDisp-plot
  fit <- glmQLFit(DGE, design)
  out <- paste(pair[[1]],"_QLDisp.png",sep="")
  png(out, width = 400, height = 400)
  plotQLDisp(fit)
  dev.off()

  qlf <- glmQLFTest(fit, contrast = c(-1,1))
  topTags(qlf)
  is.de <- decideTests(qlf, p.value=0.05)
  summary(is.de)

  # create file MD-plot
  out <- paste(pair[[1]],"_MD.png",sep="")
  png(out, width = 400, height = 400)
  plotMD(qlf, main=pair[[1]])
  dev.off()

}
