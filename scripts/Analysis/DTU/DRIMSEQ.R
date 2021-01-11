suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(DRIMSeq)
  library(stageR)
  library(DEXSeq)
  library(DESeq2)
  library(edgeR)
  library(rtracklayer)
})

options(echo=TRUE)

### ARGS
args <- commandArgs(trailingOnly = TRUE)
anname  <- args[1]
gtf     <- args[2]
outdir  <- args[3]
cmp     <- args[4]
cores   <- as.integer(args[5])
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

# define cutoffs
cutoffs <- strsplit(cutts, '-')[[1]]
pv_cut   <- as.numeric(sub("pval:", "", cutoffs[1]))
lfc_cut  <- as.numeric(sub("lfc:", "", cutoffs[2]))

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
#   original code for simple model
#   design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
#   colnames(design_full)

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

comparison_objs <- list()
setwd(outdir)

## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

  contrast_name <- strsplit(contrast,":")[[1]][1]
  contrast_groups <- strsplit(strsplit(contrast,":")[[1]][2], "-vs-")

  print(paste("Comparing ",contrast_name, sep=""))
  tryCatch({

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

    BPPARAM = MulticoreParam(workers=availablecores)

    # # reduce data for testing
    # d <- d[1:250,]

    # 1 estimate the precision,
    # 2 fit regression coefficients and perform null hypothesis testing on the coefficient of interest,
    # 3 test the coefficient associated with the difference between condition 2 and condition 1
    set.seed(1)
    system.time({
      d <- dmPrecision(d, design=design)
      d <- dmFit(d, design=design)
      d <- dmTest(d, contrast=contrast)
    })

    # add comp object to list for image
    comparison_objs <- append(comparison_objs, d)

    # calculate LFC from proportions table
    proportions <- DRIMSeq::proportions(d)
    samples_of_group_A <- subset(samples(d), condition==A)$sample_id
    samples_of_group_B <- subset(samples(d), condition==B)$sample_id
    proportions[paste(A[[1]],"mean",sep='_')] <- rowMeans(proportions[as.vector(samples_of_group_A)])
    proportions[paste(B[[1]],"mean",sep='_')] <- rowMeans(proportions[as.vector(samples_of_group_B)])
    proportions["lfc"] <- log2(proportions[paste(A[[1]],"mean",sep='_')]) - log2(proportions[paste(B[[1]],"mean",sep='_')])
    props_transcripts<-proportions[c("feature_id","lfc")]
    props_genes <- aggregate(proportions, list(proportions$gene_id), mean)[c("Group.1","lfc")]
    colnames(props_genes) <- c("gene_id", "lfc")

    # generate a single p-value per gene and transcript
    res <- DRIMSeq::results(d)
    res.txp <- DRIMSeq::results(d, level="feature")

    # filter out NA's
    no.na <- function(x) ifelse(is.na(x), 1, x)
    res$pvalue <- no.na(res$pvalue)
    res.txp$pvalue <- no.na(res.txp$pvalue)

    # put together results tables
    res <- merge(props_genes, res)
    res.txp <- merge(props_transcripts,res.txp )
    res <- res[,c(1,6,2,3,4,5)]
    res.txp <- res.txp[,c(1,7,2,3,4,5,6)]

    # Add gene names
    proportions$Gene  <- lapply(proportions$gene_id, function(x){get_gene_name(x,gtf.df)})
    res$Gene          <- lapply(res$gene_id, function(x){get_gene_name(x,gtf.df)})
    res.txp$Gene      <- lapply(res.txp$gene_id, function(x){get_gene_name(x,gtf.df)})

    proportions <- apply(proportions,2,as.character)
    res         <- apply(res,2,as.character)
    res.txp     <- apply(res.txp,2,as.character)

    # CREATE RESULTS TABLE
    write.table(as.data.frame(proportions), gzfile(paste("DTU","DRIMSEQ",contrast_name,"results_proportions.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)
    write.table(as.data.frame(res), gzfile(paste("DTU","DRIMSEQ",contrast_name,"results_genes.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)
    write.table(as.data.frame(res.txp), gzfile(paste("DTU","DRIMSEQ",contrast_name,"results_transcripts.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

    # cleanup
    rm(res,res.txp,proportions,BPPARAM)
    print(paste('cleanup done for ', contrast_name, sep=''))

    # #  stageR following DRIMSeq (included additionally in workflow by Love&Soneson&Patro)
    # pScreen <- res$pvalue
    # strp <- function(x) substr(x,1,15)
    # names(pScreen) <- strp(res$gene_id)
    #
    # pConfirmation <- matrix(res.txp$pvalue, ncol=1)
    # rownames(pConfirmation) <- strp(res.txp$feature_id)
    #
    # tx2gene <- res.txp[,c("feature_id", "gene_id")]
    # for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
    #
    # stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
    #                       pScreenAdjusted=FALSE, tx2gene=tx2gene)
    # stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=pv_cut)
    # suppressWarnings({
    #   drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
    #                                   onlySignificantGenes=TRUE)
    # })
    # # head(drim.padj)
    #
    # write.table(as.data.frame(drim.padj), gzfile(paste("DTU_DRIMSEQ",contrast_name,"results_stageR-filtered.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)
    #
    # # Post-hoc filtering
    # res.txp.filt <- DRIMSeq::results(d, level="feature")
    # smallProportionSD <- function(d, filter=0.1) {
    #   cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
    #   gene.cts <- rowsum(cts, counts(d)$gene_id)
    #   total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
    #   props <- cts/total.cts
    #   propSD <- sqrt(rowVars(props))
    #   propSD < filter
    # }
    # filt <- smallProportionSD(d)
    # res.txp.filt$pvalue[filt] <- 1
    # res.txp.filt$adj_pvalue[filt] <- 1
    #
    # write.table(as.data.frame(res.txp.filt), gzfile(paste("DTU_DRIMSEQ",contrast_name,"results_post-hoc-filtered-on-SD.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)

    }, error=function(e){
        print(warnings)
        file.create(paste("DTU","DRIMSEQ",contrast_name,"results_proportions.tsv.gz", sep="_"))
        file.create(paste("DTU","DRIMSEQ",contrast_name,"results_genes.tsv.gz", sep="_"))
        file.create(paste("DTU","DRIMSEQ",contrast_name,"results_transcripts.tsv.gz", sep="_"))
        cat("WARNING :",conditionMessage(e), "\n")
    } )
}

save.image(file = "DRIMSEQ_DTU_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
