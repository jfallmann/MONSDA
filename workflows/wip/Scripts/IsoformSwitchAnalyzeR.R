#https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#quick-start
#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#  BiocManager::install()
#}

suppressPackageStartupMessages({
    require(IsoformSwitchAnalyzeR)
})

#define notin
`%notin%` = Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)

anname      <- args[1]
countfile   <- args[2]
annotation  <- args[3]
fasta       <- args[4]
outdir      <- args[5]
cmp         <- args[6]
availablecores <- as.integer(args[7])

BPPARAM = MulticoreParam(workers=availablecores)

### MAIN ###
############

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname)))
colnames(sampleData) <- c("sampleID","condition","type","batch")
sampleData <- as.data.frame(sampleData)
#head(sampleData)

## Create design-table considering different types (paired, unpaired) and batches
if (length(levels(sampleData$type)) > 1){
  if (length(levels(sampleData$batch)) > 1){
    design = sampleData
  } else{
    design = subset(sampleData, select = -c(batch))
  }
} else{
  if (length(levels(sampleData$batch)) > 1){
    design = subset(sampleData, select = -c(type))
  } else{
    design = subset(sampleData, select = -c(batch,type))
  }
}

#get comparables
comparisons <- strsplit(cmp, ",")
print(paste("Will analyze conditions ",comparisons,sep=""))

countfiles = unlist(strsplit(countfile,','))

for(contrast in comparisons[[1]]){
  
  contrast_name <- strsplit(contrast,":")[[1]][1]
  contrast_groups <- unlist(strsplit(strsplit(contrast,":")[[1]][2], "-vs-"))
  
  print(paste("Comparing ",contrast_name, sep=""))
  
  tryCatch({
    
    print(contrast_groups)
    cf <- grep(paste(contrast_groups,collapse="|"),countfiles,value=TRUE)
    des <- droplevels(design[design$condition %in% contrast_groups,])
  
    ### Import Salmon data
    salmonQuant <- importIsoformExpression(
        sampleVector = cf,
        addIsofomIdAsColumn = TRUE
    )
    
    ####fix colnames
    colnames(salmonQuant$counts) <- sub(".sf.gz", '', basename(colnames(salmonQuant$counts)))
    colnames(salmonQuant$abundance) <- sub(".sf.gz", '', basename(colnames(salmonQuant$abundance)))

    ### Create switchAnalyzeRlist
    SwitchList <- importRdata(
      isoformCountMatrix   = salmonQuant$counts,
      isoformRepExpression = salmonQuant$abundance,
      designMatrix         = des,
      isoformExonAnnoation = annotation,
      isoformNtFasta       = fasta,
      showProgress = TRUE
    )

    ###Part 1 of Analysis
    ####Filter
    SwitchListFiltered <- preFilter(
      switchAnalyzeRlist = SwitchList,
      geneExpressionCutoff = 10,
      isoformExpressionCutoff = 3,
      removeSingleIsoformGenes = TRUE
    )
    
    ###Analyze with DEXSeq
    SwitchListDEX <- isoformSwitchTestDEXSeq(
      switchAnalyzeRlist = SwitchListFiltered,
      alpha = 0.05,
      dIFcutoff = 0.1,
      correctForConfoundingFactors=TRUE,
      overwriteIFvalues=TRUE,
      reduceToSwitchingGenes = TRUE,
      reduceFurtherToGenesWithConsequencePotential = TRUE,
      onlySigIsoforms = FALSE,
      showProgress = TRUE,
      quiet = FALSE
    )
    
    extractSwitchSummary( SwitchListDEX )
    
    ###Analyze with DRIMSeq
    SwitchListDRIM <- isoformSwitchTestDRIMSeq(
      switchAnalyzeRlist = SwitchListFiltered,
      alpha = 0.05,
      dIFcutoff = 0.1,
      testIntegration = 'intersect',
      reduceToSwitchingGenes = TRUE,
      reduceFurtherToGenesWithConsequencePotential = TRUE,
      onlySigIsoforms = FALSE,
      dmFilterArgs=list(
        min_feature_expr = 4,
        min_samps_feature_expr = min(
          SwitchListFiltered$conditions$nrReplicates
        )
      ),
      showProgress = TRUE,
      quiet = FALSE
    )
    
    extractSwitchSummary( SwitchListDRIM )
  
    ###Analyze ORFs
    SwitchListDEX <- analyzeORF(
      SwitchListDEX,
      genomeObject = NULL,
      minORFlength=100,
      orfMethod = "longest",
      cds = NULL,
      PTCDistance = 50,
      startCodons="ATG",
      stopCodons=c("TAA", "TAG", "TGA"),
      showProgress=TRUE,
      quiet=FALSE
    )
  
    SwitchListDRIM <- analyzeORF(
      SwitchListDRIM,
      genomeObject = NULL,
      minORFlength=100,
      orfMethod = "longest",
      cds = NULL,
      PTCDistance = 50,
      startCodons="ATG",
      stopCodons=c("TAA", "TAG", "TGA"),
      showProgress=TRUE,
      quiet=FALSE
    )
    
    ###Extract Sequences
    SwitchListDEX <- extractSequence(
      SwitchListDEX, 
      onlySwitchingGenes = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1,
      extractNTseq = TRUE,
      extractAAseq = TRUE,
      removeShortAAseq = TRUE,
      removeLongAAseq  = FALSE,
      alsoSplitFastaFile = FALSE,
      removeORFwithStop=TRUE,
      addToSwitchAnalyzeRlist = TRUE,
      pathToOutput = outdir,
      outputPrefix = paste("DEXSwitch_",contrast_name,sep=""),
      writeToFile = TRUE
    )
    
    SwitchListDRIM <- extractSequence(
      SwitchListDRIM, 
      onlySwitchingGenes = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1,
      extractNTseq = TRUE,
      extractAAseq = TRUE,
      removeShortAAseq = TRUE,
      removeLongAAseq  = FALSE,
      alsoSplitFastaFile = FALSE,
      removeORFwithStop=TRUE,
      addToSwitchAnalyzeRlist = TRUE,
      pathToOutput = outdir,
      outputPrefix = paste("DRIMSwitch_",contrast_name,sep=""),
      writeToFile = TRUE
    )
    
    
    ###Alternative Splicing
    SwitchListDEX <- analyzeAlternativeSplicing(
      switchAnalyzeRlist = SwitchListDEX,
      onlySwitchingGenes=TRUE,
      alpha=0.05,
      dIFcutoff = 0.1,
      showProgress=TRUE
    )
    
    SwitchListDRIM <- analyzeAlternativeSplicing(
      switchAnalyzeRlist = SwitchListDRIM,
      onlySwitchingGenes=TRUE,
      alpha=0.05,
      dIFcutoff = 0.1,
      showProgress=TRUE
    )
    
    #### Global splicing analysis
    pdf(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingSummary.pdf",sep="_"))
    extractSplicingSummary( SwitchListDEX )
    dev.off()
    pdf(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingEnrichment.pdf",sep="_"))
    extractSplicingEnrichment( SwitchListDEX )
    dev.off()
    pdf(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingGenomeWide.pdf",sep="_"))
    extractSplicingGenomeWide( SwitchListDEX )
    dev.off()
    
    pdf(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingSummary.pdf",sep="_"))
    extractSplicingSummary( SwitchListDRIM )
    dev.off()
    pdf(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingEnrichment.pdf",sep="_"))
    extractSplicingEnrichment( SwitchListDRIM )
    dev.off()
    pdf(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingGenomeWide.pdf",sep="_"))
    extractSplicingGenomeWide( SwitchListDRIM )
    dev.off()
    
    
    ###IntronRetention
    SwitchListDEX <- analyzeIntronRetention(
      SwitchListDEX,
      onlySwitchingGenes = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1,
      showProgress = TRUE,
      quiet = FALSE
    )
    
    SwitchListDRIM <- analyzeIntronRetention(
      SwitchListDRIM,
      onlySwitchingGenes = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1,
      showProgress = TRUE,
      quiet = FALSE
    )
    
    ###Get biological meaning
    consequencesOfInterest <- c('intron_retention','NMD_status','ORF_seq_similarity','tss','tts','intron_structure')
    
    SwitchListDEX <- analyzeSwitchConsequences(
      SwitchListDEX,
      consequencesToAnalyze = consequencesOfInterest,
      alpha=0.05,
      dIFcutoff=0.1,
      onlySigIsoforms=TRUE,
      ntCutoff=50,
      ntFracCutoff=NULL,
      ntJCsimCutoff=0.8,
      AaCutoff=10,
      AaFracCutoff=0.5,
      AaJCsimCutoff=0.9,
      removeNonConseqSwitches=TRUE,
      showProgress=TRUE
    )
    
    SwitchListDRIM <- analyzeSwitchConsequences(
      SwitchListDRIM,
      consequencesToAnalyze = consequencesOfInterest,
      alpha=0.05,
      dIFcutoff=0.1,
      onlySigIsoforms=TRUE,
      ntCutoff=50,
      ntFracCutoff=NULL,
      ntJCsimCutoff=0.8,
      AaCutoff=10,
      AaFracCutoff=0.5,
      AaJCsimCutoff=0.9,
      removeNonConseqSwitches=TRUE,
      showProgress=TRUE
    )
  
    ###Summarize
    
    topdex <- extractTopSwitches(
      SwitchListDEX,
      filterForConsequences = TRUE,
      extractGenes=FALSE,
      alpha=0.05,
      dIFcutoff = 0.1,
      n = Inf,
      inEachComparison=TRUE,
      sortByQvals=TRUE 
      )
    csvout <- paste('IsoformSwitchAnalyzeR_DEX_',contrast_name,'.tsv.gz', sep='')
    write.table(as.data.frame(topdex), gzfile(csvout), sep="\t", row.names=FALSE, quote=F)
    
    topdrim <- extractTopSwitches(
      SwitchListDRIM,
      filterForConsequences = TRUE,
      extractGenes=FALSE,
      alpha=0.05,
      dIFcutoff = 0.1,
      n = Inf,
      inEachComparison=TRUE,
      sortByQvals=TRUE 
    )
    csvout <- paste('IsoformSwitchAnalyzeR_DRIM_',contrast_name,'.tsv.gz', sep='')
    write.table(as.data.frame(topdrim), gzfile(csvout), sep="\t", row.names=FALSE, quote=F)
    save.image(file = paste("IsoformSwitchAnalyzeR",contrast_name,"SESSION.gz",sep="_"), version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
    
  }, error=function(e){
    rm(res,resOrdered)
    file.create(paste("IsoformSwitchAnalyseR_",contrast_name,'.tsv.gz',sep=""))
    file.create(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingSummary.pdf",sep="_"))
    file.create(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingEnrichment.pdf",sep="_"))
    file.create(paste("IsoformSwitchAnalyzeR_DEX",contrast_name,"SplicingGenomeWide.pdf",sep="_"))
    file.create(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingSummary.pdf",sep="_"))
    file.create(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingEnrichment.pdf",sep="_"))
    file.create(paste("IsoformSwitchAnalyzeR_DRIM",contrast_name,"SplicingGenomeWide.pdf",sep="_"))
        
    print(warnings)
    cat("WARNING :",conditionMessage(e), "\n")
  } )
}

save.image(file = "IsoformSwitchAnalyzeR_SESSION.gz", version = NULL, ascii = FALSE, compress = "gzip", safe = TRUE)
