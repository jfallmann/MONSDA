#https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#quick-start

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
comparison <- strsplit(cmp, ",")
print(paste("Will analyze conditions ",comparison,sep=""))

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
    reduceToSwitchingGenes=FALSE
  )
  
  extractSwitchSummary( SwitchListDEX )
  
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
  
  ORFDEX <- analyzeORF(
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
  
  ORFDRIM <- analyzeORF(
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
  
  SwitchListDEXout <- extractSequence(
    SwitchListDEX, 
    pathToOutput = paste(output,paste("DEXSwitch_",contrast_name,sep=""),sep="/"),
    writeToFile = TRUE
  )
  
  SwitchListDRIMout <- extractSequence(
    SwitchListDRIM, 
    pathToOutput = paste(output,paste("DRIMSwitch_",contrast_name,sep=""),sep="/"),
    writeToFile = TRUE
  )
  
  }, error=function(e){
    rm(res,resOrdered)
    file.create(paste("IsoformSwitchAnalyseR_",contrast_name,'.tsv.gz',sep=""))
    file.create(paste("IsoformSwitchAnalyseR",contrast_name,"MA.pdf",sep="_"))
    print(warnings)
    cat("WARNING :",conditionMessage(e), "\n")
  } )
}
