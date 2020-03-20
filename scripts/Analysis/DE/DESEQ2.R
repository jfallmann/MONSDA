#https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
suppressPackageStartupMessages({
  require(DESeq2)
  require(utils)
  require("BiocParallel")
})

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

anname <- args[1]
inname <- args[2]
outdir <- args[3]
cmp    <- args[4]
availablecores <- as.integer(args[5])

#register(MulticoreParam(availablecores))

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

#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = anno,
                              design= ~ condition+type)

#filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#print(head(dds))

#run for each pair of conditions
BPPARAM = MulticoreParam(workers=availablecores)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=BPPARAM)

                                        #Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

rld<- rlogTransformation(dds, blind=FALSE)
vsd<-varianceStabilizingTransformation(dds, blind=FALSE)

pdf(paste("DESeq2","PCA.pdf",sep="_"))
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

                                        #We also write the normalized counts to file
write.table(as.data.frame(assay(rld)), gzfile("DESeq2_rld.txt.gz"), sep="\t", col.names=NA)
write.table(as.data.frame(assay(vsd)), gzfile("DESeq2_vsd.txt.gz"), sep="\t", col.names=NA)


for(pair in comparison[[1]]){

    cname=""
    comp <- strsplit(pair,"-vs-")
    cname=pair
    print(cname)

                                        #initialize empty objects
    res <- NULL
    resOrdered <- NULL

    tryCatch({
        res <- results(dds,contrast=c("condition",as.character(comp[[1]][1]),as.character(comp[[1]][2])), parallel=TRUE, BPPARAM=BPPARAM)
                                        #sort and output
        resOrdered <- res[order(res$log2FoldChange),]
                                        #write the table to a csv file
        write.table(as.data.frame(resOrdered), gzfile(paste(cname,'_DESEQ2.csv.gz',sep="")), sep="\t")

                                        #plotMA
        pdf(paste(cname,"DESeq2_MA.pdf",sep="_"))
        plotMA(res, ylim=c(-3,3))
        dev.off()

        rm(res,resOrdered)


        print(paste('cleanup done for ', cname, sep=''))
    }, error=function(e){
        rm(res,resOrdered)
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
