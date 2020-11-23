suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(DRIMSeq)
  library(stageR)
  library(DEXSeq)
  library(DESeq2)
  library(edgeR)
})


args <- commandArgs(trailingOnly = TRUE)
anname  <- args[1]
gtf     <- args[2]
outdir  <- args[3]
cmp     <- args[4]
cores   <- as.integer(args[5])

#anname  <- file.path("/home/roberto/PROJECTS/DTU_dev/DTU/LoveSonesonPatro/Tables/ANNOTATION.gz")
#gtf     <- file.path("/home/roberto/PROJECTS/DTU_dev/GENOMES/hg38/gencode.v35.annotation.gtf.gz")
#outdir  <- "/home/roberto/PROJECTS/DTU_dev/DTU"
#cmp     <- "1_vs_2:group1-vs-group2"
#cores   <- 15





# Importing counts
samps <- read.table(file = gzfile(anname), header=TRUE, row.names=NULL)
samps$sample_id <- paste(samps$sample_id, samps$condition, sep="_")
samps$condition <- factor(samps$condition)
files <- file.path(samps$path, "quant.sf")
names(files) <- samps$sample_id

groups <- factor(samps$condition)
types <- factor(samps$type)
batches <- factor(samps$batch)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

txi <- tximport(files, type='salmon', txOut=TRUE, countsFromAbundance = "scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
cts.original <- cts
rownames(cts) <- sub("\\|.*", "", rownames(cts))

# Transcript-to-gene mapping
txdb.filename <- file.path("GENOMES/gencode.v35.annotation.sqlite")
txdb <- makeTxDbFromGFF(gtf, format="gtf")
saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# DRIMSEQ
#check for integrity -> define exception..
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)
d <- dmDSdata(counts=counts, samples=samps)

# Filter before running procedures:
#   (1) it has a count of at least 10 in at least n.small samples
#   (2) it has a relative abundance proportion of at least 0.1 in at least n.small samples
#   (3) the total count of the corresponding gene is at least 10 in all n samples

n <- nrow(samps)
n.small <- n/length(levels(samps$condition))  # its not really the smallest group, needs to be improved
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)

# shows how many of the remaining genes have N isoforms
table(table(counts(d)$gene_id))


# create designmatrix

## original code for simple model
# design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
# colnames(design_full)

## Create design-table considering different types (paired, unpaired) and batches
if (length(levels(types)) > 1){
  if (length(levels(batches)) > 1){
    design <- model.matrix(~0+groups+types+batches, data=samps)
    colnames(design) <- c(levels(groups),tl,bl)
  } else{
    design <- model.matrix(~0+groups+types, data=samps)
    colnames(design) <- c(levels(groups),tl)
  }
} else{
  if (length(levels(batches)) > 1){
    design <- model.matrix(~0+groups+batches, data=samps)
    colnames(design) <- c(levels(groups),bl)
  } else{
    design <- model.matrix(~0+groups, data=samps)
    colnames(design) <- levels(groups)
  }
}


## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

  setwd(outdir)

  contrast_name <- strsplit(contrast,":")[[1]][1]
  contrast_groups <- strsplit(strsplit(contrast,":")[[1]][2], "-vs-")

  print(paste("Comparing ",contrast_name, sep=""))



  # determine contrast
  A <- strsplit(contrast_groups[[1]][1], "\\+")
  B <- strsplit(contrast_groups[[1]][2], "\\+")
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

  d <- d[1:250,]

  # 1 estimate the precision,
  # 2 fit regression coefficients and perform null hypothesis testing on the coefficient of interest,
  # 3 test the coefficient associated with the difference between condition 2 and condition 1
  set.seed(1)
  system.time({
    d <- dmPrecision(d, design=design)
    d <- dmFit(d, design=design)
    d <- dmTest(d, contrast=contrast)
  })

  # generate a single p-value per gene and transcript
  res <- DRIMSeq::results(d)
  res.txp <- DRIMSeq::results(d, level="feature")

  # filter out NA's
  no.na <- function(x) ifelse(is.na(x), 1, x)
  res$pvalue <- no.na(res$pvalue)
  res.txp$pvalue <- no.na(res.txp$pvalue)


  # plot the estimated proportions for one of the significant genes
  idx <- which(res$adj_pvalue < 0.05)[1]
  res[idx,]
  out <- paste("DTU",contrast_name,res$gene_id[idx],"estimated_proportions.png",sep="_")
  png(out, width = 400, height = 400)
  plotProportions(d, res$gene_id[idx], "condition")
  dev.off()

  write.table(as.data.frame(res.txp), gzfile(paste("DTU",contrast_name,"results.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)
}

save.image(file = "DTU_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
