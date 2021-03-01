suppressPackageStartupMessages({
    library(tximport)
    library(GenomicFeatures)
    library(DRIMSeq)
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


anname  <- "DTU/drimseq_DTU/Tables/_ANNOTATION.gz"
gtf     <- "GENOMES/hg38/gencode.v35.annotation.gtf.gz"
outdir  <- "DTU/drimseq_DTU"
combi   <- "none"
cmp     <- "1vs2:group1-vs-group2"
cores   <- 2



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
# cutoffs <- strsplit(cutts, '-')[[1]]
# pv_cut   <- as.numeric(sub("pval:", "", cutoffs[1]))
# lfc_cut  <- as.numeric(sub("lfc:", "", cutoffs[2]))

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
# dir.create(file.path(outdir,"Figures"))

## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

    contrast <- comparisons[[1]]

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

        # BPPARAM = MulticoreParam(workers=cores)

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
    comparison_objs <- c(comparison_objs, d)

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
    res$adj_pvalue <- no.na(res$adj_pvalue)
    res.txp$adj_pvalue <- no.na(res.txp$adj_pvalue)

    res <- merge(props_genes, res)
    res.txp <- merge(props_transcripts,res.txp )
    res <- res[,c(1,6,2,3,4,5)]
    res.txp <- res.txp[,c(1,7,2,3,4,5,6)]

    proportions$Gene  <- lapply(proportions$gene_id, function(x){get_gene_name(x,gtf.df)})
    res$Gene          <- lapply(res$gene_id, function(x){get_gene_name(x,gtf.df)})
    res.txp$Gene      <- lapply(res.txp$gene_id, function(x){get_gene_name(x,gtf.df)})

    res <- res[order(res$adj_pvalue),]
    res.txp <- res.txp[order(res$adj_pvalue),]
    proportions <- proportions[order(proportions$lfc),]

    proportions.print <- as.data.frame(apply(proportions,2,as.character))
    res.print         <- as.data.frame(apply(res,2,as.character))
    res.txp.print     <- as.data.frame(apply(res.txp,2,as.character))

    # CREATE RESULTS TABLES
    # setwd(file.path(outdir,"Tables"))
    write.table(as.data.frame(proportions.print), gzfile(paste("Tables/DTU","DRIMSEQ",combi,contrast_name,"table","proportions.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)
    write.table(as.data.frame(res.print), gzfile(paste("Tables/DTU","DRIMSEQ",combi,contrast_name,"table","genes.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)
    write.table(as.data.frame(res.txp.print), gzfile(paste("Tables/DTU","DRIMSEQ",combi,contrast_name,"table","transcripts.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)
    write.table(as.data.frame(genewise_precision(d)), gzfile(paste("Tables/DTU","DRIMSEQ",combi,contrast_name,"table","genewise-precision.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

    # setwd(file.path(outdir, "Figures"))

    ## Plot feature per gene histogram
    png(paste("Figures/DTU","DRIMSEQ",combi,contrast_name,"figure","FeatPerGene.png",sep="_"))
    print(plotData(d))
    dev.off()

    ## Plot precision
    png(paste("Figures/DTU","DRIMSEQ",combi,contrast_name,"figure","Precision.png",sep="_"))
    print(plotPrecision(d))
    dev.off()

    ## Plot gene-level p-values
    png(paste("Figures/DTU","DRIMSEQ",combi,contrast_name,"figure","PValues.png",sep="_"))
    print(plotPValues(d))
    dev.off()

    # plot proportions
    figures <- data.frame()
    sigs <- which(res$adj_pvalue < 0.05)

    limit   <- 10
    counter <- 1
    message("create proportions plots")
    for(gene in sigs){
        if(counter>limit){break}
        if(is.na(gene)){next}
        suppressMessages({
            name1 <- paste("Figures/DTU","DRIMSEQ",combi,contrast_name,res$Gene[gene],"figure","plotProportions","props.png",sep="_")
            png(name1)
            print(plotProportions(d, res$gene_id[gene], group_variable = "condition"))
            dev.off()

            name2 <- paste("Figures/DTU","DRIMSEQ",combi,contrast_name,res$Gene[gene],"figure","lineplot.png",sep="_")
            png(name2)
            print(plotProportions(d, res$gene_id[gene], group_variable = "condition", plot_type = "lineplot"))
            dev.off()

            name3 <- paste("Figures/DTU","DRIMSEQ",combi,contrast_name,res$Gene[gene],"figure","ribbonplot.png",sep="_")
            png(name3)
            print(plotProportions(d, res$gene_id[gene], group_variable = "condition", plot_type = "ribbonplot"))
            dev.off()

            figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir,name1, sep="/")))
            figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir,name2, sep="/")))
            figures <- rbind(figures, c(res$gene_id[gene], res$Gene[gene], paste(outdir,name3, sep="/")))

        counter <- counter+1
        })
    }
    colnames(figures) <- c("geneID","geneName","file")
    write.table(figures, paste("Figures/DTU","DRIMSEQ",combi,contrast_name,"list","sigGenesFigures.tsv", sep="_"), sep="\t", quote=F, row.names=FALSE, col.names=TRUE)

    # cleanup
    rm(res,res.txp,proportions)
    print(paste('cleanup done for ', contrast_name, sep=''))

    })
}

save.image(file = paste("DTU","DRIMSEQ",combi,"SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
