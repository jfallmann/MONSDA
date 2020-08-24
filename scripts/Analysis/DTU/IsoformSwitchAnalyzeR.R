#https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#quick-start

 suppressPackageStartupMessages({
    require(IsoformSwitchAnalyzeR)
})

#define notin
`%notin%` = Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)

anname    <- args[1]
countfile <- args[2]
flatanno  <- args[3]
outdir    <- args[4]
cmp       <- args[5]
availablecores <- as.integer(args[6])

### MAIN ###
############

## Annotation
sampleData <- as.matrix(read.table(gzfile(anname),row.names=1))
colnames(sampleData) <- c("condition","type","batch")
sampleData <- as.data.frame(sampleData)
#head(sampleData)
comparisons <- strsplit(cmp, ",")
print(paste("Will analyze conditions ",comparisons,sep=""))

### Import Salmon data
salmonQuant <- importIsoformExpression(
    parentDir = ,
    addIsoformIdAsColumn = TRUE
)
