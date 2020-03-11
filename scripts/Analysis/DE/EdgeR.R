require(edgeR)
require(statmod)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

tbl <- args[1]
smp <- args[2]
drct <- args[3]
cmp <- args[4]


f <- function(list){
  spektrum <- rainbow(max(list),alpha=1)
  cl <- c()
  for(i in list){
    cl <- c(cl,spektrum[i])
  }
  return(cl)
}


###READ IN
data <- read.table(tbl,skip=1,h=F,sep="\t",check.names=FALSE)
header <- strsplit(readLines(tbl, n=1),"\t")[[1]]
header[1] <- "gene"
comparison <- strsplit(cmp, ",")
#######FILTER COLOUMS WITH ONLY NA
data <- Filter(function(x) !all(is.na(x)), data)
colnames(data) <- header

samples <- read.table(smp,h=F,sep="\t")
colnames(samples) <- c("name","group")

###FIND GROUPS
groups_by_number <- as.numeric(samples$group)
groups_by_name <- samples$group

###GET GENES
l <- strsplit(as.character(data$V1), ":")
genes <- unlist(lapply(l,function(x) x[length(x)]))

###GET RID OF "GENE" COLUMN
counts <- data[,-c(1)]

###CONVERT other NA TO 0
counts[is.na(counts)] <- 0

summary(counts)

###CREATE DE-GENE-LIST
dge <- DGEList(counts=counts, group=groups_by_name, samples=samples$name, genes=genes)

###REMOVE LOW EXPRESSION

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]


###Normalize by TMM

dge <- calcNormFactors(dge, method = "TMM")
tmm <- as.data.frame(cpm(dge))
colnames(tmm) <- t(dge$samples$samples)
tmm$ID <- dge$genes$genes
tmm <- tmm[c(ncol(tmm),1:ncol(tmm)-1)]
write.table(tmm, file=paste(drct,"normalized_table.tsv",sep='/'), sep="\t", quote=F, row.names=FALSE)

###plot MDS
out <- paste(drct,"MDS.png",sep="/")
png(out, width = 1000, height = 800)
colors <- f(groups_by_number)
plotMDS(dge,col=colors)
dev.off()

#design <- model.matrix(~0+Conditions)
design <- model.matrix(~0+group, data=dge$samples)
colnames(design) <- levels(dge$samples$group)

#design
###Estimate dispersion
dge <- estimateDisp(dge, design, robust=TRUE)

###plot Dispersion

out <- paste(drct,"BCV.png",sep="/")
png(out, width = 350, height = 350)
plotBCV(dge)
dev.off()

###Estimate DE
####fit glms
fit <- glmFit(dge, design)

contrast_list = list()
for(pair in comparison[[1]]) {
  con <- integer(dim(design)[2])
  comp <- strsplit(pair,"-vs-")
  a <- comp[[1]][1]
  b <- comp[[1]][2]
  m <- match(c(a,b), colnames(design))
  con[m[1]] <- -1
  con[m[2]] <- 1
  ####likelihood ratiotest
  lrt <- glmLRT(fit, contrast = con)

  ####top DE
  #topTags(lrt)

  ###DE distribution
  #summary(decideTests(lrt))

  ###plot lFC vs CPM
  print(paste("compare ", pair[[1]]))

  out <- paste(drct,'/',pair[[1]],".png",sep="")
  png(out, width = 350, height = 350)
  plotMD(lrt, main=pair[[1]])
  abline(h=c(-1, 1), col="blue")
  dev.off()
}
