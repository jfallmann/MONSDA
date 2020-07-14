#https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
suppressPackageStartupMessages({
  require(utils)
  require(BiocParallel)
  require(DESeq2)
})

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

anname <- args[1]
inname <- args[2]
outdir <- args[3]
cmp    <- args[4]
availablecores <- as.integer(args[5])

### MAIN ###
############

anno <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(anno) <- c("condition","type")
anno <- as.data.frame(anno)
#head(anno)
comparison<-strsplit(cmp, ",")
countData <- as.matrix(read.table(gzfile(inname),header=T,row.names=1))
#head(countData)

setwd(outdir)

#Check if names are consistent
if (!all(rownames(anno) %in% colnames(countData))){
    stop("Count file does not correspond to the annotation file")
}

if (length(levels(anno$types))>1){
    design <- ~0 + condition + type
} else {
    design <- ~0 + condition
}

#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = anno,
                              design= design)

#filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#print(head(dds))

#run for each pair of conditions
BPPARAM = MulticoreParam(workers=availablecores)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=BPPARAM)#, betaPrior=TRUE)

                                        #Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

rld<- rlogTransformation(dds, blind=FALSE)
vsd<-varianceStabilizingTransformation(dds, blind=FALSE)

pdf(paste("DESeq2","PCA.pdf",sep="_"))
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

                                        #We also write the normalized counts to file
write.table(as.data.frame(assay(rld)), gzfile("DESeq2_rld.txt.gz"), sep="\t", col.names=NA)
write.table(as.data.frame(assay(vsd)), gzfile("DESeq2_vsd.txt.gz"), sep="\t", col.names=NA)


for(contrast in comparison[[1]]){

    contrast_name <- strsplit(contrast,":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast,":")[[1]][2], "-vs-")

    message(paste("Comparing ",contrast_name, sep=""))

    tryCatch({

                                        # determine contrast
        A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"),use.names=FALSE)
        B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"),use.names=FALSE)

        tempa <- droplevels(anno[anno$condition %in% A,])
        tempb <- droplevels(anno[anno$condition %in% B,])

        plus <- 1/length(A)
        minus <- 1/length(B)*-1

        BPPARAM = MulticoreParam(workers=availablecores)

                                        #initialize empty objects
        res=""
        resOrdered=""

        res <- results(dds,contrast=list(paste('condition',levels(tempa$condition),sep=''),paste('condition',levels(tempb$condition),sep='')), listValues=c(plus,minus), parallel=TRUE, BPPARAM=BPPARAM)

                                        #sort and output
        resOrdered <- res[order(res$log2FoldChange),]

                                        #write the table to a csv file
        write.table(as.data.frame(resOrdered), gzfile(paste(contrast_name,'_DESeq2.csv.gz',sep="")), sep="\t")

                                        #plotMA
        pdf(paste(contrast_name,"DESeq2_MA.pdf",sep="_"))
        plotMA(res, ylim=c(-3,3))
        dev.off()

        rm(res,resOrdered, BPPARAM)


        print(paste('cleanup done for ', contrast_name, sep=''))
    }, error=function(e){
        rm(res,resOrdered)
        file.create(paste(contrast_name,'_DESeq2.csv.gz',sep=""))
        print(warnings)
        cat("WARNING :",conditionMessage(e), "\n")
    } )
}

#Here we choose blind so that the initial conditions setting does not influence the outcome, ie we want to see if the conditions cluster based purely on the individual datasets, in an unbiased way. According to the documentation, the rlogTransformation method that converts counts to log2 values is apparently better than the old varienceStabilisation method when the data size factors vary by large amounts.
par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'black')
pdf(paste("DESeq2_VST","and_log2.pdf",sep="_"))
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type='l', lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

##############################
library('RColorBrewer')
library('gplots')
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
pdf(paste("DESeq2","heatmap1.pdf",sep="_"))
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10,6))
dev.off()
pdf(paste("DESeq2","heatmap2.pdf",sep="_"))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 6))
dev.off()
pdf(paste("DESeq2","heatmap3.pdf",sep="_"))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(10, 6))
dev.off()

#The above shows heatmaps for 30 most highly expressed genes (not necessarily the biggest fold change). The data is of raw counts (left), regularized log transformation (center) and from variance stabilizing transformation (right) and you can clearly see the effect of the transformation has by shrinking the variance so that we don’t get the squish effect shown in the left hand graph.
##############################
#Now we calculate sample to sample distances so we can make a dendrogram to look at the clustering of samples.
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
#updated in latest vignette (See comment by Michael Love)
#this line was incorrect
#heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(16, 16))
#From the Apr 2015 vignette
hc <- hclust(distsRL)
pdf(paste("DESeq2","heatmap_samplebysample.pdf",sep="_"))
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))
dev.off()

##############################
rm(rld, vsd)

save.image(file = "DESeq2_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
