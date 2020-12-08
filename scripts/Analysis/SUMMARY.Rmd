```{r header, include=FALSE}
author    <- Sys.getenv("LOGNAME")
pres_date <- Sys.Date()
```

---
title:  "NextSnakes SUMMARY"
author: `r author`
date:   `r pres_date`
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(DRIMSeq)
  library(stageR)
  library(DEXSeq)
  library(DESeq2)
  library(edgeR)
})

# args <- commandArgs(trailingOnly = TRUE)
# input     <- args[1]
# formats   <- args[2]
# outdir    <- args[3]
# cutoffs    <- args[4]

input     <- "DTU/DRIMSEQ+DE/EDGER"
formats   <- "html+pdf"
outdir    <- "SUMMARY"
cutoffs   <- "DTU=pvalue:0.5+lfc:1.5"
condis    <- "mikelove-group1-out1+mikelove-group1-out2+mikelove-group1-out3+mikelove-group2-out1+mikelove-group2-out2+mikelove-group2-out3"

setwd("PROJECTS/DTU_dev")

analyses <- strsplit(input,'+', fixed=TRUE)[[1]]
works <- unique(sub("\\/.*", '', analyses))

```

***
# OVERVIEW
input: `r input`\
formats: `r formats`\
cutoffs: `r cutoffs`


```{r DE, echo=FALSE, results='asis', eval='DE' %in% works}
cat("
***
# DIFFERENTIAL GENE EXPRESSION
")
```

```{r DE_EDGER, echo=FALSE, results='asis', eval='DE/EDGER' %in% analyses}
cat("
## edgeR
results of DE EDGER analysis: \n
............\n
")

```

```{r DTU, echo=FALSE, results='asis', eval='DTU' %in% works}
cat("
***
# DIFFERENTIAL TRANSCRIPT USAGE
")
```

```{r DTU_DRIMSEQ, echo=FALSE, results='asis', eval='DTU/DRIMSEQ' %in% analyses}
cat("
## DRIMSeq
results of DTU DRIMSeq analysis:\n
............\n
")
list.files('DTU/DRIMSEQ')

```

```{r DTU_DEXSEQ, echo=FALSE, results='asis', eval='DTU/DEXSEQ' %in% analyses}
cat("
## DEXSeq

results of DTU DEXSeq analysis:\n
............\n
")
setwd('DTU/DEXSEQ')
list.files()
```
