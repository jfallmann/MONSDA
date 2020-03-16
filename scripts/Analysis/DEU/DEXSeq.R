## From https://github.com/vivekbhr/Subread_to_DEXSeq/blob/master/load_SubreadOutput.R
## Load Fcount output from : DEXSeq_after_Fcount.R into DEXSeq
## Copyright 2015 Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de). Licence: GPLv3.

suppressPackageStartupMessages({
    require(DEXSeq)
    require(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

anname    <- args[1]
countfile <- args[2]
flatanno  <- args[3]
outdir    <- args[4]
cmp       <- args[5]
availablecores <- as.integer(args[6])

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type")
sampleData <- as.data.frame(sampleData)
#head(sampleData)
comparison<-strsplit(cmp, ",")
## Combinations of conditions
#condcomb<-as.data.frame(combn(unique(sampleData$condition),2))[1:2,]
##countfile <- as.matrix(read.table(gzfile(inname),header=T,row.names=1))
##head(countData)

## Read Fcount output and convert to dxd
DEXSeqDataSetFromFeatureCounts <- function (countfile, sampleData,
                                            design = ~sample + exon + type:exon + condition:exon, flattenedfile = NULL)

{
    ##  Take a fcount file and convert it to dcounts for dexseq
    message("Reading and adding Exon IDs for DEXSeq")

    read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
    colnames(dcounts) <- c("GeneID", rownames(sampleData) )
    id <- as.character(dcounts[,1])
    n <- id
    split(n,id) <- lapply(split(n ,id), seq_along )
    rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
    dcounts <- dcounts[,2:ncol(dcounts)]

    dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name
    ##filter low counts
    keep <- rowSums(dcounts) >= 10
    dcounts <- dcounts[keep,]

    ## get genes and exon names out
    splitted <- strsplit(rownames(dcounts), ":")
    exons <- sapply(splitted, "[[", 2)
    genesrle <- sapply(splitted, "[[", 1)

    ## parse the flattened file
    if (!is.null(flattenedfile)) {
        aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE,
                                 header = FALSE)
        colnames(aggregates) <- c("chr", "source", "class", "start",
                                  "end", "ex", "strand", "ex2", "attr")
        aggregates$strand <- gsub("\\.", "*", aggregates$strand)
        aggregates <- aggregates[which(aggregates$class == "exonic_part"), # exonic_part
                                 ]
        aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
        aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1",
                                  aggregates$attr)
        # trim the gene_ids to 255 chars in order to match with featurecounts
        longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
        warning(paste0(longIDs,
                       " aggregate geneIDs were found truncated in featureCounts output"),
                call. = FALSE)
        aggregates$gene_id <- substr(aggregates$gene_id,1,255)

        transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1",
                            aggregates$attr)
        transcripts <- strsplit(transcripts, "\\+")
        exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", # exonic_part_number
                        aggregates$attr)
        exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start,
                                                                  end = aggregates$end), strand = aggregates$strand)
        names(exoninfo) <- paste(aggregates$gene_id, exonids,
                                 sep = ":E")

        names(transcripts) <- names(exoninfo)
        if (!all(rownames(dcounts) %in% names(exoninfo))) {
            stop("Count files do not correspond to the flattened annotation file")
        }
        matching <- match(rownames(dcounts), names(exoninfo))
        stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
        stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons,
                             genesrle, exoninfo[matching], transcripts[matching])
        return(dxd)
    }
    else {
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons,
                             genesrle)
        return(dxd)
    }

}

### MAIN ###
#read in count table and normalize

dxd = DEXSeqDataSetFromFeatureCounts(countfile, sampleData, design = ~sample + exon + condition:exon, flattenedfile = flatanno)

setwd(outdir)

print(paste('Will run DEXSeq with ',availablecores,' cores',sep=''))

BPPARAM = MulticoreParam(workers=availablecores)

for (pair in comparison[[1]]){#n in 1:ncol(condcomb)){

    cname=""
    #cname=paste(condcomb[,n],collapse='_vs_')
    comp <- strsplit(pair,"-vs-")
    cname=pair
    print(cname)

                                        #initialize empty objects
    dxdpair <- NULL
    dxr1 <- NULL

    tryCatch({

                                        #dxdpair = dxd[,which(dxd$condition == condcomb[1,n] | dxd$condition == condcomb[2,n])]#, drop=True]
        dxdpair = dxd[,which(dxd$condition == as.character(comp[[1]][1]) | dxd$condition == as.character(comp[[1]][2]))]#, drop=True]

        dxdpair = estimateDispersions( dxdpair, BPPARAM=BPPARAM)

        pdf(paste("DEXSeq",cname,"DispEsts.pdf",sep="_"))
        plotDispEsts( dxdpair )
        dev.off()

        dxdpair = testForDEU( dxdpair,BPPARAM=BPPARAM )

        dxdpair = estimateExonFoldChanges( dxdpair, fitExpToVar="condition", BPPARAM=BPPARAM)

        dxr1 = DEXSeqResults( dxdpair )

        rm(dxdpair)

        csvout <- paste(paste('DEXSeqResults',cname,sep='_'),'.tsv.gz', sep='')
        write.table(as.data.frame(dxr1), gzfile(csvout), sep="\t")

        htmlout <- paste(paste('DEXSeq',cname,sep='_'),'.html', sep='')
        pathout <- paste('DEXSeqReport',cname,sep='_')
        DEXSeqHTML( dxr1, FDR=0.1, color=c("#FF000080", "#0000FF80"), path=pathout, file=htmlout, BPPARAM=BPPARAM)

        rm(dxr1)
        print(paste('cleanup done for ', cname, sep=''))


    }, error=function(e){
        rm(dxdpair,dxr1)
        cat("WARNING :",conditionMessage(e), "\n")
    })
}
