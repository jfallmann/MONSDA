##read from commandline
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
tab <- args[1]
sample <- paste(strsplit(tab, "_")[[1]][1],strsplit(tab, "_")[[1]][3],sep="_")

sample

library(ggplot2)
data <- read.table(tab, head=T, sep="\t")
colnames(data) <- sub("X","U",colnames(data))
library(reshape2)
library(plyr)
require(scales)

meltdata <- melt(data, id.vars = "Position")

meltdata$Distbin <- findInterval(meltdata$Position, c(seq(min(meltdata$Position),max(meltdata$Position))),all.inside=F)

p <- ggplot(meltdata, aes(x = factor(round_any(Distbin,1)), y = value, fill = variable)) + geom_boxplot() + theme_bw()
p <- p + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=1,show_guide = FALSE)
p <- p + facet_wrap(~ variable)
p <- p + theme(aspect.ratio=0.5)
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(sample)
p <- p + xlab("Position")
p <- p + ylab("Probability of being unpaired")
p
out<-paste("Accessibility_Profile_AllUs_",sample,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)
out<-paste("Accessibility_Profile_AllUs_",sample,".eps",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)
out<-paste("Accessibility_Profile_AllUs_",sample,".jpg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)

data$Distbin <- findInterval(data$Position, c(seq(min(data$Position),max(data$Position))),all.inside=F)
p <- ggplot(data, aes(x = factor(round_any(Distbin,1)), y = U9, fill = "firebrick")) + geom_boxplot() + theme_bw()
p <- p + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=3,show_guide = FALSE)
p <- p + scale_fill_manual(values = c("firebrick")) + guides(fill=FALSE)
p <- p + theme(aspect.ratio=0.5)
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(sample)
p <- p + xlab("Position")
p <- p + ylab("Probability of being unpaired")
p
out<-paste("Accessibility_Profile_U9_",sample,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)
out<-paste("Accessibility_Profile_U9_",sample,".eps",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)
out<-paste("Accessibility_Profile_U9_",sample,".jpg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)
q()
