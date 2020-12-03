args <- commandArgs(trailingOnly = TRUE)

input     <- args[1]
formats   <- args[2]
outdir    <- args[3]
cutoff    <- args[4]

print(input)
print(formats)
print(outdir)
print(cutoff)

setwd(outdir)
test <- data.frame(2)
write.table(as.data.frame(test), paste("rmarkdown_summary.txt", sep=""), sep="\t", quote=F, row.names=FALSE)
# write.table(as.data.frame(test), gzfile(paste(outdir,"rmarkdown_summary.",formats,sep="_")), sep="\t", quote=F, row.names=FALSE)
