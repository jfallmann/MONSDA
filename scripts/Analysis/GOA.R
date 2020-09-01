#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("topGO")

suppressPackageStartupMessages({
    require(topGO)
    require(Rgraphviz)
})

#define notin
`%notin%` = Negate(`%in%`)

args <- commandArgs(TRUE)
background <- args[1]
test <- args[2]
GOs <- args[3]

#run GO
geneID2GO <- readMappings(GOs, sep = "\t", IDsep = ",")
expressedGenes <- read.table(background,sep="\t")
Genes <- read.table(test,sep="\t")
GenesOI <- Genes$V1
geneList <- factor(as.integer(expressedGenes$V1 %in% GenesOI))
names(geneList) <- expressedGenes$V1

GOdata <- new("topGOdata",
              ontology = "MF",
              allGenes = geneList,
              geneSel = GenesOI,
              annot = annFUN.gene2GO,  # the new annotation function
              gene2GO = geneID2GO)    ## the gene ID to GOs dataset

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
pvalFis <- score(resultFisher)
resultWeight <- getSigGroups(GOdata, test.stat)
pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
cor(pvalFis, pvalWeight)
geneData(resultWeight)
allRes <- GenTable(GOdata, classic = resultFisher, weight = resultWeight, orderBy = "weight", ranksOf = "classic", topNodes = 20)
write.table(allRes, file = paste("TopGO_MF",test,sep="_"), append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = paste("TopGO_MFGraph",test,sep="_"), useInfo = "all", pdfSW = TRUE)

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = GenesOI,
              annot = annFUN.gene2GO,  # the new annotation function
              gene2GO = geneID2GO)    ## the gene ID to GOs dataset

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
pvalFis <- score(resultFisher)
resultWeight <- getSigGroups(GOdata, test.stat)
pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
cor(pvalFis, pvalWeight)
geneData(resultWeight)
allRes <- GenTable(GOdata, classic = resultFisher, weight = resultWeight, orderBy = "weight", ranksOf = "classic", topNodes = 20)
write.table(allRes, file = paste("TopGO_BP",test,sep="_"), append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = paste("TopGO_BPGraph",test,sep="_"), useInfo = "all", pdfSW = TRUE)


