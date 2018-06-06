#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to create the subsample infofile with the total sample infofile.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("sample", nargs=1, help="the infomation file of total sample.")
myparser$add_argument("groups", nargs=1, help="the groups wanted to extract.")
myparser$add_argument("-o","--output", default="subsample_info.txt", help="the output file, default is subsample_info.txt.")

args <- myparser$parse_args()

input <- args$sample
gr <- as.character(args$groups)
output <- args$output
groupids <- unlist(strsplit(gr, "[,]")[1])

sampleda <- read.table(input, header=T, sep="\t", comment.char="")
dd <- data.frame()
newgroupname <- vector()
for (i in 1:length(groupids)){
	newgroupname[i] <- groupids[i]
	tmp <- as.vector(sampleda[[groupids[i]]])
	if (i == 1){
		dd <- data.frame(cbind(as.vector((sampleda[,1])), tmp))
	}else{ 
		dd <- data.frame(cbind(dd, tmp))
	}
}
colnames(dd) <- c("sample", "group")
write.table(dd, output, sep="\t", row.names=F, col.names=T, quote=F)
