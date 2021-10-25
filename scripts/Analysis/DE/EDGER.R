suppressPackageStartupMessages({
    require(utils)
    require(BiocParallel)
    require(edgeR)
    require(rtracklayer)
    require(RUVSeq)
    require(dplyr)
    require(GenomeInfoDb)
})

options(echo=TRUE)

## ARGS
args <- commandArgs(trailingOnly = TRUE)
argsLen         <- length(args);
anname          <- args[1]
countfile       <- args[2]
gtf             <- args[3]
outdir          <- args[4]
cmp             <- args[5]
combi           <- args[6]
availablecores  <- as.integer(args[7])
spike           <- if (argsLen > 7) args[8] else ''

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
print(paste('Run EdgeR DE with ', availablecores, ' cores', sep=''))

# set thread-usage
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
comparison <- strsplit(cmp, ",")

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
for(contrast in comparison[[1]]){

    contrast_name <- strsplit(contrast, ":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast, ":")[[1]][2], "-vs-")

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
            #colnames(design) <- c(levels(sampleData$condition), tl, bl)
        } else{
            des <- ~type+condition
            design <- model.matrix(des, data=sampleData)
            #colnames(design) <- c(levels(condition), tl)
        }
    } else{
        if (length(unique(subset(sampleData, A == condition)$batch)) > 1 | length(unique(subset(sampleData, B == condition)$batch)) > 1){
            des <- ~batch+condition
            design <- model.matrix(des, data=sampleData)
            #colnames(design) <- c(levels(sampleData$condition), bl)
        } else{
            des <- ~condition
            design <- model.matrix(des, data=sampleData)
            #colnames(design) <- levels(sampleData$condition)
        }
    }
    print(design)

    ## check genes and spike-ins
    if (spike != ''){
        print("Spike-in used, data will be normalized to spike in separately")
        spike = strsplit(spike, "=")[[1]][2]
        setwd(WD)
        ctrlgenes <- readLines(spike)
        setwd(outdir)
        counts_norm <-RUVg(newSeqExpressionSet(as.matrix(countData)), ctrlgenes, k=1)
        genes <- rownames(countData)
        countData <- countData %>% subset(!row.names(countData) %in% ctrlgenes) # removing spike-ins for standard analysis

        sampleData_norm <- cbind(sampleData, pData(counts_norm))
        design_norm <- model.matrix(as.formula(paste(deparse(des), colnames(pData(counts_norm))[1], sep=" + ")), data=sampleData_norm)
        #colnames(design_norm) <- c(colnames(design),"W_1")

        dge_norm <- DGEList(counts=counts(counts_norm), group=sampleData$condition, samples=samples, genes=genes)

        ## filter low counts
        keep <- filterByExpr(dge_norm)
        dge_norm <- dge_norm[keep, , keep.lib.sizes=FALSE]

        #relevel to base condition B
        dge_norm$samples$group <- relevel(dge_norm$samples$group, ref = B[[1]])

        ## normalize with TMM
        dge_norm <- calcNormFactors(dge_norm, method = "TMM", BPPARAM=BPPARAM)

        ## create file normalized table
        tmm_norm <- as.data.frame(cpm(dge_norm))
        colnames(tmm_norm) <- t(dge_norm$samples$samples)
        tmm_norm$ID <- dge_norm$genes$genes
        tmm_norm <- tmm_norm[c(ncol(tmm_norm), 1:ncol(tmm_norm)-1)]

        write.table(as.data.frame(tmm_norm), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "DataSet", "table", "Normalized_norm.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        ## create file MDS-plot with and without summarized replicates
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "MDS_norm.png", sep="_")
        png(out)
        plotMDS(dge_norm, col = as.numeric(dge_norm$samples$group), cex = 1)
        dev.off()

        ## estimate Dispersion
        dge_norm <- estimateDisp(dge_norm, design_norm, robust=TRUE)

        ## create file BCV-plot - visualizing estimated dispersions
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "BCV_norm.png", sep="_")
        png(out)
        plotBCV(dge_norm)
        dev.off()

        ## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
        fit_norm <- glmQLFit(dge_norm, design_norm, robust=TRUE)

        ## create file quasi-likelihood-dispersion-plot
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "QLDisp_norm.png", sep="_")
        png(out)
        plotQLDisp(fit_norm)
        dev.off()

    }

    # Same without spike-in normalization
    genes <- rownames(countData)
    dge <- DGEList(counts=countData, group=sampleData$condition, samples=samples, genes=genes)

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
    tmm <- tmm[c(ncol(tmm), 1:ncol(tmm)-1)]

    write.table(as.data.frame(tmm), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "DataSet", "table", "Normalized.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

    ## create file MDS-plot with and without summarized replicates
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "MDS.png", sep="_")
    png(out)
    plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1)
    dev.off()

    ## estimate Dispersion
    dge <- estimateDisp(dge, design, robust=TRUE)

    ## create file BCV-plot - visualizing estimated dispersions
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "BCV.png", sep="_")
    png(out)
    plotBCV(dge)
    dev.off()

    ## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
    fit <- glmQLFit(dge, design, robust=TRUE)

    ## create file quasi-likelihood-dispersion-plot
    out <- paste("Figures/DE", "EDGER", combi, contrast_name, "DataSet", "figure", "QLDisp.png", sep="_")
    png(out)
    plotQLDisp(fit)
    dev.off()

    tryCatch({

        # determine contrast, only for complex cases, not needed for our pairwise comparisons now
        #A <- strsplit(contrast_groups[[1]][1], "\\+")
        #B <- strsplit(contrast_groups[[1]][2], "\\+")
        #minus <- 1/length(A[[1]])*(-1)
        #plus <- 1/length(B[[1]])
        #contrast <- cbind(integer(dim(design)[2]), colnames(design))
        #for(i in A[[1]]){
        #    contrast[which(contrast[,2]==i)]<- minus
        #}
        #for(i in B[[1]]){
        #    contrast[which(contrast[,2]==i)]<- plus
        #}
        #contrast <- as.numeric(contrast[,1])

        ## Testing
        #qlf <- glmQLFTest(fit, contrast=contrast) ## glm quasi-likelihood-F-Test
        AvsB <- makeContrasts(TreatvsUntreat = paste("condition", A, sep=""), levels=design)
        qlf <- glmQLFTest(fit, contrast=AvsB) ## glm quasi-likelihood-F-Test
        # add comp object to list for image
        comparison_objs[[contrast_name]] <- qlf

        # # Add gene names  (check how gene_id col is named )
        qlf$table$Gene  <- lapply(rownames(qlf$table) , function(x){get_gene_name(x, gtf_gene)})
        qlf$table$Gene_ID <- rownames(qlf$table)
        res <- qlf$table[, c(6,5,1,2,3,4)]
        res <- as.data.frame(apply(res,2, as.character))

        # create results table
        write.table(as.data.frame(res), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "results.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        # create sorted results Tables
        tops <- topTags(qlf, n=nrow(qlf$table), sort.by="logFC")
        tops <- tops$table[, c(7,6,4,2,3,5,8)]
        tops <- as.data.frame(apply(tops,2, as.character))
        write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsLogFCsorted.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        tops <- topTags(qlf, n=nrow(qlf$table), sort.by="PValue")
        tops <- tops$table[, c(7,6,4,2,3,5,8)]
        tops <- as.data.frame(apply(tops,2, as.character))
        write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsPValueSorted.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

        ## plot lFC vs CPM
        out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "MD.png", sep="_")
        png(out)
        plotMD(qlf, main=contrast_name)
        abline(h=c(-1, 1), col="blue")
        dev.off()

        if (spike != ''){  # Same for spike-in normalized
            rm(qlf, tops, res)

            # determine contrast, only for complex cases, not needed for our pairwise comparisons now
            #A <- strsplit(contrast_groups[[1]][1], "\\+")
            #B <- strsplit(contrast_groups[[1]][2], "\\+")
            #minus <- 1/length(A[[1]])*(-1)
            #plus <- 1/length(B[[1]])
            #contrast <- cbind(integer(dim(design)[2]), colnames(design))
            #for(i in A[[1]]){
            #    contrast[which(contrast[,2]==i)]<- minus
            #}
            #for(i in B[[1]]){
            #    contrast[which(contrast[,2]==i)]<- plus
            #}
            #contrast <- as.numeric(contrast[,1])

            ## Testing
            #qlf <- glmQLFTest(fit, contrast=contrast) ## glm quasi-likelihood-F-Test
            AvsB <- makeContrasts(TreatvsUntreat = paste("condition", A, sep=""), levels=design)
            qlf <- glmQLFTest(fit, contrast=AvsB) ## glm quasi-likelihood-F-Test
            # add comp object to list for image
            comparison_objs <- append(comparison_objs, qlf)

            # # Add gene names  (check how gene_id col is named )
            qlf$table$Gene  <- lapply(rownames(qlf$table), function(x){get_gene_name(x, gtf_gene)})
            qlf$table$Gene_ID <- rownames(qlf$table)
            res <- qlf$table[, c(6,5,1,2,3,4)]
            res <- as.data.frame(apply(res,2, as.character))

            # create results table
            write.table(as.data.frame(res), gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "results_norm.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

            # create sorted results Tables
            tops <- topTags(qlf, n=nrow(qlf$table), sort.by="logFC")
            tops <- tops$table[, c(7,6,4,2,3,5,8)]
            tops <- as.data.frame(apply(tops,2, as.character))
            write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsLogFCsorted_norm.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

            tops <- topTags(qlf, n=nrow(qlf$table), sort.by="PValue")
            tops <- tops$table[, c(7,6,4,2,3,5,8)]
            tops <- as.data.frame(apply(tops,2, as.character))
            write.table(tops, gzfile(paste("Tables/DE", "EDGER", combi, contrast_name, "table", "resultsPValueSorted_norm.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

            ## plot lFC vs CPM
            out <- paste("Figures/DE", "EDGER", combi, contrast_name, "figure", "MD_norm.png", sep="_")
            png(out)
            plotMD(qlf, main=contrast_name)
            abline(h=c(-1, 1), col="blue")
            dev.off()
        }

        # cleanup
        rm(qlf, res, tops)
        print(paste('cleanup done for ', contrast_name, sep=''))
    })
}

# Same for all samples without design specific normalization
## name types and levels for design

## Create design-table considering different types (paired, unpaired) and batches
des <- ~0+condition
design <- model.matrix(des, data=sampleData_all)
colnames(design) <- levels(sampleData_all$condition)
print(design)

genes <- rownames(countData_all)
samples <- rownames(sampleData_all)

dge <- DGEList(counts=countData_all, group=sampleData_all$condition, samples=samples, genes=genes)

## filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

## normalize with TMM
dge <- calcNormFactors(dge, method = "TMM", BPPARAM=BPPARAM)

## create file normalized table
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm), 1:ncol(tmm)-1)]

write.table(as.data.frame(tmm), gzfile(paste("Tables/DE", "EDGER", combi, "DataSet", "table", "AllConditionsNormalized.tsv.gz", sep="_")), sep="\t", quote=F, row.names=FALSE)

## create file MDS-plot with and without summarized replicates
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsMDS.png", sep="_")
png(out)
plotMDS(dge, col = as.numeric(dge$samples$group), cex = 1)
dev.off()

## estimate Dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

## create file BCV-plot - visualizing estimated dispersions
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsBCV.png", sep="_")
png(out)
plotBCV(dge)
dev.off()

## fitting a quasi-likelihood negative binomial generalized log-linear model to counts
fit <- glmQLFit(dge, design, robust=TRUE)

## create file quasi-likelihood-dispersion-plot
out <- paste("Figures/DE", "EDGER", combi, "DataSet", "figure", "AllConditionsQLDisp.png", sep="_")
png(out)
plotQLDisp(fit)
dev.off()

save.image(file = paste("DE_EDGER", combi, "SESSION.gz", sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
