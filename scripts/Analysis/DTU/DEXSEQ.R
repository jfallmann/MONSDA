suppressPackageStartupMessages({
    require(BiocParallel)
    library(tximport)
    library(GenomicFeatures)
    library(DRIMSeq)
    library(DEXSeq)
    library(rtracklayer)
})

options(echo=TRUE)

### ARGS
args <- commandArgs(trailingOnly = TRUE)
anname  <- args[1]
gtf     <- args[2]
outdir  <- args[3]
combi   <- args[4]
cmp     <- args[5]
cores   <- as.integer(args[6])
print(args)


### FUNCS
get_gene_name <- function(id, df){
    name_list <- df$gene_name[df['gene_id'] == id]
    if(length(unique(name_list)) == 1){
        return(name_list[1])
    }else{
        message(paste("WARNING: ambigous gene id: ",id))
        return (paste("ambigous",id,sep="_"))
    }
}

### SCRIPT
# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)

BPPARAM = MulticoreParam(workers=cores)

# Importing counts
samps <- read.table(file = gzfile(anname), header=TRUE, row.names=NULL)
samps$sample_id <- paste(samps$sample_id, samps$condition, sep="_")
samps$condition <- factor(samps$condition)
files <- file.path(samps$path, "quant.sf.gz")
names(files) <- samps$sample_id

groups <- factor(samps$condition)
types <- factor(samps$type)
batches <- factor(samps$batch)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

## check combi
if (combi == "none"){
    combi <- ''
}

txi <- tximport(files, type='salmon', txOut=TRUE, countsFromAbundance = "scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
cts.original <- cts
rownames(cts) <- sub("\\|.*", "", rownames(cts))

# Transcript-to-gene mapping
txdb.filename <- file.path(paste(gtf, "sqlite", sep="."))
txdb <- makeTxDbFromGFF(gtf, format="gtf")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# DRIMSEQ
#check for integrity -> define exception..
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts, check.names=FALSE)
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

## shows how many of the remaining genes have N isoforms
#table(table(counts(d)$gene_id))

# create designmatrix
#   original code for simple model
design <- model.matrix(~groups, data=DRIMSeq::samples(d))

## name types and levels for design
# bl <- sapply("batch", paste0, levels(batches)[-1])
# tl <- sapply("type", paste0, levels(types)[-1])

## Create design-table considering different types (paired, unpaired) and batches
# if (length(levels(types)) > 1){
#    if (length(levels(batches)) > 1){
#        design <- model.matrix(~0+groups+types+batches, data=samps)
#        colnames(design) <- c(levels(groups),tl,bl)
#    } else{
#        design <- model.matrix(~0+groups+types, data=samps)
#        colnames(design) <- c(levels(groups),tl)
#    }
# } else{
#    if (length(levels(batches)) > 1){
#        design <- model.matrix(~0+groups+batches, data=samps)
#        colnames(design) <- c(levels(groups),bl)
#    } else{
#        design <- model.matrix(~0+groups, data=samps)
#        colnames(design) <- levels(groups)
#    }
# }

comparison_objs <- list()
setwd(outdir)

## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

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

    # DEXSeq on DRIMSeq filtering
    sample.data <- DRIMSeq::samples(d)
    count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
    # The design formula of the DEXSeqDataSet here uses the language “exon” but this should be read as “transcript” for our analysis.
    dxd <- DEXSeqDataSet(countData=count.data,
                sampleData=sample.data,
                design=~sample + exon + condition:exon,
                featureID=counts(d)$feature_id,
                groupID=counts(d)$gene_id
            )
    dxd <- DEXSeq::estimateSizeFactors(dxd) 
    dxd <- DEXSeq::estimateDispersions(dxd, quiet=TRUE)
    dxd <- testForDEU(dxd, reducedModel=~sample + exon)

    dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
    qval <- perGeneQValue(dxr)
    dxr.g <- data.frame(gene=names(qval),qval)

    columns <- c("featureID","groupID","pvalue","padj")
    dxr <- as.data.frame(dxr[,columns])

    dxr$Gene  <- lapply(dxr$groupID, function(x){get_gene_name(x,gtf.df)})
    dxr <- dxr[,c(1,2,5,3,4)]
    dxr <- apply(dxr,2,as.character)

    write.table(as.data.frame(dxr), gzfile(paste("Tables/DTU","DEXSEQ",combi,contrast_name,"table","results.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)
    
    # # stageR following DEXSeq
    # strp <- function(x) substr(x,1,15)
    # pConfirmation <- matrix(dxr$pvalue,ncol=1)
    # dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
    # pScreen <- qval
    # names(pScreen) <- strp(names(pScreen))
    # tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
    # for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
    #
    # stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
    #                       pScreenAdjusted=TRUE, tx2gene=tx2gene)
    # stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
    # suppressWarnings({
    #   dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
    #                                  onlySignificantGenes=TRUE)
    # })

    # write.table(as.data.frame(dex.padj), gzfile(paste("DTU_DEXSEQ",contrast_name,"stageR-filtered.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)
}

save.image(file = paste("DTU_DEXSEQ", combi, "SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
