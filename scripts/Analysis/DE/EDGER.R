suppressPackageStartupMessages({
    require(BiocParallel)
    require(dplyr)
    require(edgeR)
    require(rtracklayer)
    require(RUVSeq)
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
spike           <- args[8]

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

### SCRIPT
print(paste('Run EdgeR DE with ',availablecores,' cores',sep=''))

# set thread-usage
BPPARAM = MulticoreParam(workers=availablecores)

# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)

## Annotation
sampleData <- as.data.frame(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("group","type","batch")
sampleData <- as.data.frame(sampleData)
groups <- factor(sampleData$group)
samples <- rownames(sampleData)
types <- factor(sampleData$type)
batches <- factor(sampleData$batch)

## Combinations of conditions
comparisons <- strsplit(cmp, ",")

## check combi
if (combi == "none"){
    combi <- ''
}

## readin counttable
countData <- read.table(countfile,header = TRUE,row.names=1)

# Check if names are consistent
if (!all(rownames(sampleData) %in% colnames(countData))){
    stop("Count file does not correspond to the annotation file")
}

## name types and levels for design
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


## get genes names and spike-ins

genes <- rownames(countData)
if (spike != ''){
    ctrlgenes <- readLines(spike)
    countData <- countData %>% filter(row.names(countData) %!in% ctrlgenes)
    genes <- rownames(countData)
    counts_norm <-RUVg(countData, ctrlgenes, k=1)

    dge_norm <- DGEList(counts=counts(counts_norm), group=groups, samples=samples, genes=genes)

    ## filter low counts
    keep <- filterByExpr(dge_norm)
    dge_norm <- dge_norm[keep, , keep.lib.sizes=FALSE]

    ## normalize with TMM
    dge_norm <- calcNormFactors(dge_norm, method = "TMM", BPPARAM=BPPARAM)


    ## create file normalized table
    tmm_norm <- as.data.frame(cpm(dge_norm))
    colnames(tmm_norm) <- t(dge_norm$samples$samples)
    tmm_norm$ID <- dge_norm$genes$genes
    tmm_norm <- tmm_norm[c(ncol(tmm_norm),1:ncol(tmm_norm)-1)]

    setwd(outdir)
    write.table(as.data.frame(tmm_norm), gzfile(paste("Tables/    DE","EDGER",combi,"DataSet","table","AllConditionsNormalized_norm.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)

    ## create file MDS-plot with and without summarized replicates
    out <- paste("Figures/    DE","EDGER",combi,"DataSet","figure","AllConditionsMDS_norm.png", sep="_")
    png(out)
    plotMDS(dge_norm, col = as.numeric(dge_norm$samples$group), cex = 1)
    dev.off()

    ## estimate Dispersion
    dge_norm <- estimateDisp(dge_norm, design, robust=TRUE)

    ## create file BCV-plot - visualizing estimated dispersions
    out <- paste("Figures/    DE","EDGER",combi,"DataSet","figure","AllConditionsBCV_norm.png", sep="_")
    png(out)
    plotBCV(dge_norm)
    dev.off()

    ## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
    fit_norm <- glmQLFit(dge_norm, design, robust=TRUE)

    ## create file quasi-likelihood-dispersion-plot
    out <- paste("Figures/    DE","EDGER",combi,"DataSet","figure","AllConditionsQLDisp_norm.png", sep="_")
    png(out)
    plotQLDisp(fit_norm)
    dev.off()

 }

# Same without spike-in normalization
dge <- DGEList(counts=countData, group=groups, samples=samples, genes=genes)

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

setwd(outdir)
write.table(as.data.frame(tmm), gzfile(paste("Tables/DE","EDGER",combi,"DataSet","table","AllConditionsNormalized.tsv.gz",sep="_")), sep="\t", quote=F, row.names=FALSE)

## create file MDS-plot with and without summarized replicates
out <- paste("Figures/DE","EDGER",combi,"DataSet","figure","AllConditionsMDS.png", sep="_")
png(out)
plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1)
dev.off()

## estimate Dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- paste("Figures/DE","EDGER",combi,"DataSet","figure","AllConditionsBCV.png", sep="_")
png(out)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust=TRUE)

## create file quasi-likelihood-dispersion-plot
out <- paste("Figures/DE","EDGER",combi,"DataSet","figure","AllConditionsQLDisp.png", sep="_")
png(out)
plotQLDisp(fit)
dev.off()

comparison_objs <- list()

## Analyze according to comparison groups
for(compare in comparisons[[1]]){

    contrast_name <- strsplit(compare,":")[[1]][1]
    contrast_groups <- strsplit(strsplit(compare,":")[[1]][2], "-vs-")

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

        # Testing
        qlf <- glmQLFTest(fit, contrast=contrast) ## glm quasi-likelihood-F-Test

        # add comp object to list for image
        comparison_objs <- append(comparison_objs, qlf)

        # # Add gene names  (check how gene_id col is named )
        qlf$table$Gene  <- lapply(rownames(qlf$table) , function(x){get_gene_name(x,gtf.df)})
        qlf$table$Gene_ID <- rownames(qlf$table)
        res <- qlf$table[,c(6,5,1,2,3,4)]
        res <- as.data.frame(apply(res,2,as.character))

        # create results table
        write.table(as.data.frame(res), gzfile(paste("Tables/DE","EDGER",combi,contrast_name,"table","results.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        # create sorted results Tables
        tops <- topTags(qlf, n=nrow(qlf$table), sort.by="logFC")
        tops$table$Gene  <- lapply(rownames(tops) , function(x){get_gene_name(x,gtf.df)})
        tops$table$Gene_ID <- rownames(tops$table)
        tops <- tops$table[,c(8,7,2,3,4,5,6)]
        tops <- as.data.frame(apply(tops,2,as.character))
        write.table(tops, file=paste("Tables/DE","EDGER",combi,contrast_name,"table","resultsLogFCsorted.tsv.gz", sep="_"), sep="\t", quote=F, row.names=FALSE)
        tops <- topTags(qlf, n=nrow(qlf$table), sort.by="PValue")
        tops$table$Gene  <- lapply(rownames(tops) , function(x){get_gene_name(x,gtf.df)})
        tops$table$Gene_ID <- rownames(tops$table)
        tops <- tops$table[,c(8,7,2,3,4,5,6)]
        tops <- as.data.frame(apply(tops,2,as.character))
        write.table(tops, file=paste("Tables/DE","EDGER",combi,contrast_name,"table","resultsPValueSorted.tsv.gz", sep="_"), sep="\t", quote=F, row.names=FALSE)

        ## plot lFC vs CPM
        out <- paste("Figures/DE","EDGER",combi,contrast_name,"figure","MD.png",sep="_")
        png(out)
        plotMD(qlf, main=contrast_name)
        abline(h=c(-1, 1), col="blue")
        dev.off()


        if (spike != ''){  # Same for spike-in normalized
            rm(qlf, tops, res)
            # Testing
            qlf <- glmQLFTest(fit_norm, contrast=contrast) ## glm quasi- likelihood-F-Test

            # add comp object to list for image
            comparison_objs <- append(comparison_objs, qlf)

            # # Add gene names  (check how gene_id col is named )
            qlf$table$Gene  <- lapply(rownames(qlf$table) ,     function(x){get_gene_name(x,gtf.df)})
            qlf$table$Gene_ID <- rownames(qlf$table)
            res <- qlf$table[,c(6,5,1,2,3,4)]
            res <- as.data.frame(apply(res,2,as.character))

            # create results table
            write.table(as.data.frame(res), gzfile(paste("Tables/    DE","EDGER",combi,contrast_name,"table","results_norm.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

            # create sorted results Tables
            tops <- topTags(qlf, n=nrow(qlf$table), sort.by="logFC")
            tops$table$Gene  <- lapply(rownames(tops) ,     function(x){get_gene_name(x,gtf.df)})
            tops$table$Gene_ID <- rownames(tops$table)
            tops <- tops$table[,c(8,7,2,3,4,5,6)]
            tops <- as.data.frame(apply(tops,2,as.character))
            write.table(tops, file=paste("Tables/    DE","EDGER",combi,contrast_name,"table","resultsLogFCsorted_norm.tsv.gz", sep="_"), sep="\t", quote=F, row.names=FALSE)
            tops <- topTags(qlf, n=nrow(qlf$table), sort.by="PValue")
            tops$table$Gene  <- lapply(rownames(tops) ,     function(x){get_gene_name(x,gtf.df)})
            tops$table$Gene_ID <- rownames(tops$table)
            tops <- tops$table[,c(8,7,2,3,4,5,6)]
            tops <- as.data.frame(apply(tops,2,as.character))
            write.table(tops, file=paste("Tables/    DE","EDGER",combi,contrast_name,"table","resultsPValueSorted_norm.tsv.gz", sep="_"), sep="\t", quote=F, row.names=FALSE)

            ## plot lFC vs CPM
            out <- paste("Figures/    DE","EDGER",combi,contrast_name,"figure","MD_norm.png",sep="_")
            png(out)
            plotMD(qlf, main=contrast_name)
            abline(h=c(-1, 1), col="blue")
            dev.off()

        }

        # cleanup
        rm(qlf, res, tops, BPPARAM)
        print(paste('cleanup done for ', contrast_name, sep=''))
    })
}

save.image(file = paste("DE_EDGER",combi,"SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
