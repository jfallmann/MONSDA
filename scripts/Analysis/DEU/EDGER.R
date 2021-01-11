
suppressPackageStartupMessages({
    require(BiocParallel)
    require(dplyr)
    require(edgeR)
})

options(echo=TRUE)

## ARGS
args <- commandArgs(trailingOnly = TRUE)
anname          <- args[1]
countfile       <- args[2]
outdir          <- args[3]
cmp             <- args[4]
availablecores  <- as.integer(args[5])
print(args)

## FUNCS
get_gene_name <- function(id, df){
  name_list <- df$gene_name[df['gene_id'] == id]
  if(length(unique(name_list)) == 1){
    return(name_list[1])
  }else{
    message(paste("WARNING: ambigous gene id: ",id))
    return (paste("ambigous",id,sep="_"))
  }
}

RainbowColor <- function(groups){
  groupsAsNumbers <- as.numeric(groups)
  spektrum <- rainbow(max(groupsAsNumbers),alpha=1)
  cl <- c()
  for(i in groupsAsNumbers){
    cl <- c(cl,spektrum[i])
  }
  return(cl)
}


### SCRIPT
print(paste('Will run EdgeR DEU with ',availablecores,' cores',sep=''))

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type","batch")
sampleData <- as.data.frame(sampleData)
groups <- factor(sampleData$condition)
samples <- rownames(sampleData)
types <- factor(sampleData$type)
batches <- factor(sampleData$batch)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

## readin counttable
countData <- read.table(countfile,header = TRUE,row.names=1)

#Check if names are consistent
if (!all(rownames(sampleData) %in% colnames(countData))){
    stop("Count file does not correspond to the annotation file")
}

## get genes names out
genes <- rownames(countData)

splitted <- strsplit(rownames(countData), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)

dge <- DGEList(counts=countData, group=groups, samples=samples, genes=genesrle)

## Addd exons to dge
dge$genes$exons <- exons

## filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

## normalize with TMM
dge <- calcNormFactors(dge, method = "TMM", BPPARAM=BPPARAM)

##name types and levels for design
bl <- sapply("batch",paste0,levels(batches)[-1])
tl <- sapply("type",paste0,levels(types)[-1])

## Create design-table considering different types (paired, unpaired) and batches
if (length(levels(types)) > 1){
    if (length(levels(batches)) > 1){
        design <- model.matrix(~0+groups+types+batches, data=sampleData)
        colnames(design) <- c(levels(groups),tl,bl)
    } else{
        design <- model.matrix(~0+groups+types, data=sampleData)
        colnames(design) <- c(levels(groups),tl)
    }
} else{
    if (length(levels(batches)) > 1){
        design <- model.matrix(~0+groups+batches, data=sampleData)
        colnames(design) <- c(levels(groups),bl)
    } else{
        design <- model.matrix(~0+groups, data=sampleData)
        colnames(design) <- levels(groups)
    }
}

print(paste('FITTING DESIGN: ',design, sep=""))

## create file normalized table
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]

setwd(outdir)

write.table(as.data.frame(tmm), gzfile("EDGER_DEU_All_Conditions_normalized.tsv.gz"), sep="\t", quote=F, row.names=FALSE)

## create file MDS-plot with and without sumarized replicates
out <- "EDGER_DEU_All_Conditions_MDS.png"
png(out, width = 400, height = 400)
colors <- RainbowColor(dge$samples$group)
plotMDS(dge, col=colors)
dev.off()
if (length(levels(groups)) > 2){
    print("Will plot MDS for Count sums")
    DGEsum <- sumTechReps(dge, ID=groups)
    out <- "EDGER_DE_All_Conditions_sum_MDS.png"
    png(out, width = 400, height = 400)
    colors <- RainbowColor(DGEsum$samples$group)
    plotMDS(DGEsum, col=colors)
    dev.off()
}

## estimate Dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- "EDGER_DEU_All_Conditions_BCV.png"
png(out, width = 400, height = 400)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust=TRUE)

## create file quasi-likelihood-dispersion-plot
out <- "EDGER_DEU_All_Conditions_QLDisp.png"
png(out, width = 400, height = 400)
plotQLDisp(fit)
dev.off()

comparison_objs <- list()

## Analyze according to comparison groups
for(contrast in comparisons[[1]]){

    contrast_name <- strsplit(contrast,":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast,":")[[1]][2], "-vs-")

    BPPARAM = MulticoreParam(workers=availablecores)

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

        # quasi-likelihood F-Test
        qlf <- glmQLFTest(fit, contrast = contrast)
        is.de <- decideTests(qlf, p.value=0.05)
        summary(is.de)

        # add comp object to list for image
        comparison_objs <- append(comparison_objs, qlf)

        # create results table
        write.table(as.data.frame(qlf), gzfile(paste("DEU","EDGER",contrast_name,"results.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        # # create file MD-plot
        # out <- paste("EDGER_DEU_",contrast_name,"_MD.png",sep="")
        # png(out, width = 400, height = 400)
        # plotMD(qlf, main=contrast_name)
        # dev.off()

        # cleanup
        rm(qlf, BPPARAM)
        print(paste('cleanup done for ', contrast_name, sep=''))

<<<<<<< HEAD
        save.image(file = paste("EDGER_DEU",contrast_name,"SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)

=======
>>>>>>> DTU
    }, error=function(e){
        print(warnings)
        file.create(paste("DEU","EDGER",contrast_name,"results.tsv.gz", sep="_"))
        # file.create(paste("EDGER_DEU_",contrast_name,"_MD.png",sep=""))
        cat("WARNING :",conditionMessage(e), "\n")
    } )
}

save.image(file = "EDGER_DEU_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
