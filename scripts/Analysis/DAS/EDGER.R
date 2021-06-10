suppressPackageStartupMessages({
    require(dplyr)
    require(BiocParallel)
    require(edgeR)
    require(rtracklayer)
})

options(echo=TRUE)

## ARGS
args <- commandArgs(trailingOnly = TRUE)
anname          <- args[1]
countfile       <- args[2]
gtf             <- args[3]
outdir          <- args[4]
combi           <- args[5]
cmp             <- args[6]
availablecores  <- as.integer(args[7])
print(args)

### FUNCS
get_gene_name <- function(id, df){
    name_list <- df$gene_name[df['gene_id'] == id]
    if(length(unique(name_list)) == 1){
        return(name_list[1])
    }else{
        message(paste("WARNING: ambigous gene id: ", id))
        return (paste(unique(name_list), sep="|"))
    }
}


### SCRIPT
print(paste('Will run EdgeR DAS with ',availablecores,' cores',sep=''))

## set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)
gtf_gene <- droplevels(subset(gtf.df, type == "gene"))

## Annotation
sampleData_all <- as.data.frame(read.table(gzfile(anname), row.names=1))
colnames(sampleData_all) <- c("condition", "type", "batch")
sampleData_all$condition <- as.factor(sampleData_all$condition)
sampleData_all$batch <- as.factor(sampleData_all$batch)
sampleData_all$type <- as.factor(sampleData_all$type)
samples <- rownames(sampleData_all)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

## check combi
if (combi == "none"){
    combi <- ''
}

## readin counttable
countData_all <- read.table(countfile, header = TRUE, row.names=1)

#Check if names are consistent
if (!all(rownames(sampleData_all) %in% colnames(countData_all))){
    stop("Count file does not correspond to the annotation file")
}

comparison_objs <- list()

WD <- getwd()
setwd(outdir)

## Analyze according to comparison groups
for(compare in comparisons[[1]]){

    contrast_name <- strsplit(compare, ":")[[1]][1]
    contrast_groups <- strsplit(strsplit(compare, ":")[[1]][2], "-vs-")

    print(paste("Comparing ", contrast_name, sep=""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names=FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names=FALSE)

    #subset Datasets for pairwise comparison
    countData <- cbind(countData_all[ , grepl( paste(B, '_', sep='') , colnames( countData_all ) )], countData_all[ ,grepl(paste(A, '_', sep='') , colnames( countData_all )) ])
    sampleData <- droplevels(rbind(subset(sampleData_all, B == condition), subset(sampleData_all, A == condition)))
    sampleData$condition <- relevel(sampleData$condition, ref=B)
    samples <- rownames(sampleData)
    ## name types and levels for design
    bl <- sapply("batch", paste0, levels(sampleData$batch)[1:length(levels(sampleData$batch))-1])
    tl <- sapply("type", paste0, levels(sampleData$type)[1:length(levels(sampleData$type))-1])

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == condition)$type)) > 1 | length(unique(subset(sampleData, B == condition)$type)) > 1){
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1){
            des <- ~type+batch+condition
            design <- model.matrix(des, data=sampleData)
            # colnames(design) <- c(levels(sampleData$condition), tl, bl)
        } else{
            des <- ~type+condition
            design <- model.matrix(des, data=sampleData)
            # colnames(design) <- c(levels(condition), tl)
        }
    } else{
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1){
            des <- ~batch+condition
            design <- model.matrix(des, data=sampleData)
            # colnames(design) <- c(levels(sampleData$condition), bl)
        } else{
            des <- ~condition
            design <- model.matrix(des, data=sampleData)
            # colnames(design) <- levels(sampleData$condition)
        }
    }
    print(design)

    ## get genes names out
    genes <- rownames(countData)

    ## get genes and exon names out
    splitted <- strsplit(rownames(countData), ":")
    exons <- sapply(splitted, "[[", 2)
    genesrle <- sapply(splitted, "[[", 1)

    dge <- DGEList(counts=countData, group=sampleData$condition, samples=samples, genes=genesrle)

    ## Addd exons to dge
    dge$genes$exons <- exons

    ## filter low counts
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]

    #relevel to base condition B
    dge$samples$group <- relevel(dge$samples$group, ref = B[[1]])

    ## normalize with TMM
    dge <- calcNormFactors(dge, method = "TMM", BPPARAM=BPPARAM)

    ## create file normalized table
    tmm <- as.data.frame(cpm(dge))
    colnames(tmm) <- t(dge$samples$samples)
    tmm$ID <- dge$genes$genes
    tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]

    tmm.genes <- tmm
    tmm.genes$Gene <- lapply(tmm.genes$ID, function(x){get_gene_name(x,gtf.df)})
    tmm.genes <- as.data.frame(apply(tmm.genes,2,as.character))
    write.table(as.data.frame(tmm.genes), gzfile(paste("Tables/DAS","EDGER",combi, contrast_name, "DataSet","table","Normalized.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)

    ## create file MDS-plot with and without sumarized replicates
    out <- paste("Figures/DAS","EDGER",combi, contrast_name, "DataSet","figure","MDS.png", sep="_")
    png(out)
    plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1)
    dev.off()

    ## estimate Dispersion
    dge <- estimateDisp(dge, design, robust=TRUE)

    ## create file BCV-plot - visualizing estimated dispersions
    out <- paste("Figures/DAS","EDGER",combi, contrast_name, "DataSet","figure","BCV.png", sep="_")
    png(out)
    plotBCV(dge)
    dev.off()

    ## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
    fit <- glmQLFit(dge, design, robust=TRUE)

    ## create file quasi-likelihood-dispersion-plot
    out <- paste("Figures/DAS","EDGER",combi, contrast_name, "DataSet","figure","QLDisp.png", sep="_")
    png(out)
    plotQLDisp(fit)
    dev.off()

    tryCatch({

        # # determine contrast
        # A <- strsplit(contrast_groups[[1]][1], "\\+")
        # B <- strsplit(contrast_groups[[1]][2], "\\+")
        # minus <- 1/length(A[[1]])*(-1)
        # plus <- 1/length(B[[1]])
        # contrast <- cbind(integer(dim(design)[2]), colnames(design))
        # for(i in A[[1]]){
        #   contrast[which(contrast[,2]==i)]<- minus
        # }
        # for(i in B[[1]]){
        #   contrast[which(contrast[,2]==i)]<- plus
        # }
        # contrast <- as.numeric(contrast[,1])

        # create files topSpliced by gene, simes and exon method
        # sp <- diffSpliceDGE(fit, contrast=contrast, geneid="genes", exonid="exons", verbose=FALSE)
        sp <- diffSpliceDGE(fit, geneid="genes", exonid="exons", verbose=FALSE)

        # add comp object to list for image
        comparison_objs <- append(comparison_objs, sp)

        tops <- topSpliceDGE(sp, test="gene", n=length(fit$counts))
        tops$Gene  <- lapply(strsplit(rownames(tops), split = ":")[[1]][1] , function(x){get_gene_name(x,gtf.df)})
        tops$Gene_ID <- rownames(tops)
        tops <- tops[,c(7,6,3,4,5,2)]
        tops <- as.data.frame(apply(tops,2,as.character))
        write.table(as.data.frame(tops), gzfile(paste("Tables/DAS","EDGER",combi,contrast_name, "table","resultsGeneTest.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        tops <- topSpliceDGE(sp, test="simes", n=length(fit$counts))
        tops$Gene  <- lapply(strsplit(rownames(tops), split = ":")[[1]][1] , function(x){get_gene_name(x,gtf_gene)})
        tops$Gene_ID <- rownames(tops)
        tops <- tops[,c(6,5,2,3,4)]
        tops <- as.data.frame(apply(tops,2,as.character))
        write.table(as.data.frame(tops), gzfile(paste("Tables/DAS","EDGER",combi,contrast_name, "table","resultsSimesTest.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        tops <- topSpliceDGE(sp, test="exon", n=length(fit$counts))
        tops$Gene  <- lapply(strsplit(rownames(tops), split = ":")[[1]][1] , function(x){get_gene_name(x,gtf_gene)})
        tops$Gene_ID <- rownames(tops)
        tops <- tops[,c(8,7,4,5,6,2,3)]
        tops <- as.data.frame(apply(tops,2,as.character))
        write.table(as.data.frame(tops), gzfile(paste("Tables/DAS","EDGER",combi,contrast_name, "table","resultsDiffSpliceExonTest.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        # create files diffSplicePlots
        tops <- topSpliceDGE(sp, test="simes", n=length(fit$counts))
        figures <- data.frame("geneID" = character(), "exonID" = character(), "file"=character(),stringsAsFactors=FALSE)
        sigs <- which(tops$P.Value < 0.05)

        limit   <- 100
        counter <- 1
        message("create splice plots")
        for(gene in sigs){
            if(counter>limit){break}

            name1 <- paste("Figures/DAS","EDGER",combi,contrast_name,tops$genes[gene],"figure","transcripts.png",sep="_")
            png(name1)
            print(plotSpliceDGE(sp, geneid=tops$genes[gene], genecol="genes"))
            dev.off()
            figures[counter, ] <- c(rownames(tops)[gene], tops$genes[gene], paste(outdir,name1, sep="/"))

            counter <- counter+1
        }
    write.table(figures, paste("Figures/DAS","EDGER",combi,contrast_name, "list","topSpliceSimes.tsv", sep="_"), sep="\t", quote=F, row.names=FALSE, col.names=TRUE)

    # cleanup
    rm(lrt,tops,BPPARAM)
    print(paste('cleanup done for ', contrast_name, sep=''))

    })
}

## Create design-table considering different types (paired, unpaired) and batches
des <- ~0+condition
design <- model.matrix(des, data=sampleData_all)
colnames(design) <- levels(sampleData_all$condition)
print(design)
## get genes names out
genes <- rownames(countData_all)
samples <- rownames(sampleData_all)

## get genes and exon names out
splitted <- strsplit(rownames(countData), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)

dge <- DGEList(counts=countData_all, group=sampleData_all$condition, samples=samples, genes=genesrle)

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

tmm.genes <- tmm
tmm.genes$Gene <- lapply(tmm.genes$ID, function(x){get_gene_name(x,gtf_gene)})
tmm.genes <- as.data.frame(apply(tmm.genes,2,as.character))
write.table(as.data.frame(tmm.genes), gzfile(paste("Tables/DAS","EDGER",combi, "DataSet","table","AllConditionsNormalized.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)

## create file MDS-plot with and without sumarized replicates
out <- paste("Figures/DAS","EDGER",combi, "DataSet","figure","AllConditionsMDS.png", sep="_")
png(out)
plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1)
dev.off()

## estimate Dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- paste("Figures/DAS","EDGER",combi, "DataSet","figure","AllConditionsBCV.png", sep="_")
png(out)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust=TRUE)

## create file quasi-likelihood-dispersion-plot
out <- paste("Figures/DAS","EDGER",combi, "DataSet","figure","AllConditionsQLDisp.png", sep="_")
png(out)
plotQLDisp(fit)
dev.off()


save.image(file = paste("DAS_EDGER",combi, "SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
