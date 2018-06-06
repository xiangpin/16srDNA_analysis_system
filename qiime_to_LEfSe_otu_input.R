#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to create the otu relative abundance table for LEfSe analysis.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab txt file of the qiime created.")
myparser$add_argument("samplefile", nargs=1, help="the sample files.")
myparser$add_argument("-o","--output", default="LEfSe_otu_input.txt", help="the output file, default is LEfSe_otu_input.txt.")
args <- myparser$parse_args()

otu <- args$otu_tab
sample <- args$samplefile
out<- args$output

a<-read.table(otu, header=T, skip=1,sep = '\t', row.names=1, check.names=F, comment.char="")
sampleTable <- read.table(sample, header=T, check.names=F)
a$taxonomy <- NULL
#a <- (a/colSums(a))*100
sample <- vector()
sample <- sampleTable[,1]
#head(groupTable)
sampleList <- names(a)
groupList <- vector()
for (i in 1:length(sampleList)) {
        groupList[i] <- as.character(sampleTable[sampleTable$sample==sampleList[i],2])
}
names(a) <- groupList

write.table(a, out, sep="\t", row.names=T, col.names=T, quote=F)
