## Modified from https://github.com/vivekbhr/Subread_to_DEXSeq/blob/master/load_SubreadOutput.R
## Load Fcount output from : DEXSeq_after_Fcount.R into DEXSeq
## Copyright 2015 Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de). Licence: GPLv3.

suppressPackageStartupMessages({
    require(dplyr)
    require(BiocParallel)
    require(DEXSeq)
    require(rtracklayer)
})

options(echo=TRUE)

## define notin
`%notin%` = Negate(`%in%`)

## ARGS
args <- commandArgs(trailingOnly = TRUE)

anname          <- args[1]
countfile       <- args[2]
gtf             <- args[3]
flatanno        <- args[4]
outdir          <- args[5]
combi           <- args[6]
cmp             <- args[7]
availablecores  <- as.integer(args[8])
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
# load gtf
gtf.rtl <- rtracklayer::import(gtf)
gtf.df <- as.data.frame(gtf.rtl)

BPPARAM = MulticoreParam(workers=availablecores)

## Annotation
sampleData <- as.data.frame(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type","batch")
#head(sampleData)
comparisons <- strsplit(cmp, ",")
print(paste("Will analyze conditions ",comparisons,sep=""))

## check combi
if (combi == "none"){
    combi <- ''
}

## Create design-table considering different types (paired, unpaired) and batches
if (length(levels(sampleData$type)) > 1){
    if (length(levels(sampleData$batch)) > 1){
        full = ~ sample + exon + type:exon + batch:exon + condition:exon
        reduced = ~ sample + exon + type:exon + batch:exon
    } else{
        full = ~ sample + exon + type:exon + condition:exon
        reduced = ~ sample + exon + type:exon
    }
} else{
    if (length(levels(sampleData$batch)) > 1){
        full = ~ sample + exon + batch:exon + condition:exon
        reduced = ~ sample + exon + batch:exon
    } else{
        full = ~ sample + exon + condition:exon
        reduced = ~ sample + exon
    }
}

print(paste('FITTING DESIGN: ',full, reduced, sep=""))

## Read Fcount output and convert to dxd
DEXSeqDataSetFromFeatureCounts <- function (countfile, sampleData, design = design, flattenedfile = NULL)

{
    ##  Take a fcount file and convert it to dcounts for dexseq
    print("Reading and adding Exon IDs for DEXSeq")

    read.table(countfile,skip = 2) %>% dplyr::arrange(V1,V3,V4) -> dcounts
    colnames(dcounts) <- c("GeneID", rownames(sampleData) )
    id <- as.character(dcounts[,1])
    n <- id
    split(n,id) <- lapply(split(n ,id), seq_along )
    rownames(dcounts) <- sprintf("%s%s%03.f",id, ":E",as.numeric(n))
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
        aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
        colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
        aggregates$strand <- gsub("\\.", "*", aggregates$strand)
        aggregates <- aggregates[which(aggregates$class == "exonic_part"),]
        aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
        aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
        # trim the gene_ids to 255 chars in order to match with featurecounts
        longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
        warning(paste0(longIDs, " aggregate geneIDs were found truncated in featureCounts output"), call. = FALSE)
        aggregates$gene_id <- substr(aggregates$gene_id,1,255)

        transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
        transcripts <- strsplit(transcripts, "\\+")
        exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
        exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
        names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")

        names(transcripts) <- names(exoninfo)
        if (!all(rownames(dcounts) %in% names(exoninfo))) {
            print(head(rownames(dcounts)))
            print(head(names(exoninfo)))
            stop("Count files do not correspond to the flattened annotation file")
        }
        matching <- match(rownames(dcounts), names(exoninfo))
        stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
        stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle, exoninfo[matching], transcripts[matching])
        return(dxd)
    }else{
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle)
        return(dxd)
    }
}

### MAIN ###
#read in count table and normalize

dxd = DEXSeqDataSetFromFeatureCounts(countfile, sampleData, design = full, flattenedfile = flatanno)

comparison_objs = list()

setwd(outdir)
print(paste('Will run DEXSeq with ',availablecores,' cores',sep=''))

for(contrast in comparisons[[1]]){

    contrast_name <- strsplit(contrast, ":")[[1]][1]
    contrast_groups <- strsplit(strsplit(contrast, ":")[[1]][2], "-vs-")

    print(paste("Comparing ",contrast_name, sep=""))

    # initialize empty objects
    dxdpair=""
    dxr1=""

    tryCatch({

        # determine contrast
        A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"),use.names=FALSE)
        B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"),use.names=FALSE)
        C <- c(A,B)

        dxdpair = dxd[, dxd$condition %in% C]

        dxdpair = estimateSizeFactors( dxdpair )
        dxdpair = estimateDispersions( dxdpair, BPPARAM=BPPARAM )

        png(paste("Figures/DEU","DEXSEQ",combi,contrast_name, "figure","DispEsts.png",sep="_"))
        plotDispEsts(dxdpair)
        dev.off()

        dxdpair = testForDEU( dxdpair, reducedModel=reduced, fullModel=full, BPPARAM=BPPARAM )

        dxdpair = estimateExonFoldChanges( dxdpair, fitExpToVar="condition", BPPARAM=BPPARAM)

        dxr1 = DEXSeqResults( dxdpair )

        comparison_objs[[contrast_name]] = dxr1

        png(paste("Figures/DEU","DEXSEQ",combi,contrast_name, "figure","plotMA.png",sep="_"))
        print(plotMA( dxr1, cex=0.8))
        dev.off()

        out <- paste('Tables/DEU','DEXSEQ',combi,contrast_name,'table','results.tsv.gz', sep='_')
        write.table(as.data.frame(dxr1), gzfile(out), sep="\t", row.names=FALSE, quote=F)

        htmlout <- paste('DEXSeq_',combi,'_',contrast_name,'.html', sep='')
        pathout <- paste('DEXSeqReport',combi,contrast_name,sep='_')
        DEXSeqHTML( dxr1, FDR=0.1, color=c("#FF000080", "#0000FF80"), path=pathout, file=htmlout, BPPARAM=BPPARAM)

        figures <- data.frame("geneID" = character(), "dxr1ID" = character(), "file"=character(),stringsAsFactors=FALSE)
        sigs <- as.data.frame(which(dxr1$padj < 0.01))
        colnames(sigs) <- c("id")
        sigs$padj <- dxr1$padj[sigs$id]
        sigs <- sigs$id[order(sigs$padj)]
        limit   <- 100
        counter <- 1
        message("create transcripts plots")
        for(gene in sigs){
            if(counter>limit){break}
            name1 <- paste("Figures/DEU","DEXSEQ",combi,contrast_name,counter, "figure","transcripts.png",sep="_")
            png(name1, height = 1000, width = 1000)
            print(plotDEXSeq( dxr1, dxr1$groupID[gene],legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ))
            dev.off()
            figures[counter, ] <- c(dxr1$groupID[gene], as.character(gene), paste(outdir,name1, sep="/"))
            # figures <- rbind(figures, c(as.character(gene), dxr1$groupID[gene], paste(outdir,name1, sep="/")))

            counter <- counter+1
        }
        write.table(figures, paste("Figures/DEU","DEXSEQ",combi,contrast_name, "list","sigGroupsFigures.tsv", sep="_"), sep="\t", quote=F, row.names=FALSE, col.names=TRUE)

        # cleanup
        rm(dxdpair, dxr1)
        print(paste('cleanup done for ', contrast_name, sep=''))
    })
}

save.image(file = paste("DEU_DEXSEQ",combi, "SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
