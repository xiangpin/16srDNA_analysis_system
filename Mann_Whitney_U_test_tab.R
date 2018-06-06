#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to do the Mann-Whitney U test.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("Feature_tab", nargs=1, help="the Feature_tab file.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o","--output", default="wilcox_test_out.txt", help="the output file, default is wilcox_test_out.txt.")
args <- myparser$parse_args()

otu <- args$Feature_tab
groupfile <- args$sample
output <- args$output
library("reshape2")

da<-read.table(otu, header=T, row.names=1, check.names=FALSE, sep='\t', comment.char="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="", sep="\t")
positive <- as.character(levels(groupTable[,2])[1])
another <- as.character(levels(groupTable[,2])[2])
da$taxonomy <- NULL

positivesample <- as.vector(groupTable[groupTable$group==positive,]$sample)
anothersample <- as.vector(groupTable[groupTable$group==another,]$sample)

wilctest <- function(x){
	(wilcox.test(as.numeric(x[positivesample]), as.numeric(x[anothersample])))$p.value
	#print(wilcox.test(as.numeric(x[positivesample]), as.numeric(x[anothersample])))
	}
dat <- da
dat <- setNames(dat, as.vector(groupTable[,2]))
datmean <- sapply(split.default(dat, names(dat)), rowMeans)
colnames(datmean) <- sapply(colnames(datmean), paste, "mean", sep="_")
da$p.value <- apply(da, 1, wilctest)
da$FDR <- p.adjust(as.numeric(da$p.value), method = "BH")
da <- da[da$p.value<=0.05&!is.na(da$p.value),]
#da<-da[da$FDR<=0.05&!is.na(da$FDR),]

dat <- merge(datmean,da, by=0)
colnames(dat)[1] <- "Feature"
write.table(dat, output, sep="\t", row.names=F, col.names=T, quote=F)
