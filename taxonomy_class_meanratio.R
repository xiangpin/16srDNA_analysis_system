#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to caculate the different class level relative abundance.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
#myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-l", "--level", default=8, help="the tax levels to want to plot, default is 8, total.")
myparser$add_argument("-o", "--data", default="taxoutput_ratio.txt", help="the output data, default is taxoutput_ratio.txt.")
#myparser$add_argument("-o","--output", default="specie_bar_plot.svg", help="the output file, default is specie_bar_plot.svg.")

args <- myparser$parse_args()


library(ggplot2)
library(plyr)
library(reshape2)
otu <- args$otu_tab
#groupfile <- args$sample
#outfigure <- args$output
outputdata <- args$data

level <- as.integer(args$level)
da<-read.table(otu, header=T, skip=1, row.names=1,check.names=FALSE, sep='\t', comment.char="")
#groupTable <- read.table(groupfile, header=T, sep="\t", check.names=F, comment.char="")
tax <- as.vector(da$taxonomy)
species <- vector()
for (i in 1:length(tax)){
	if (level==8){
		tmp <- tax[i]
		species[i] <- tmp
	}
	else {
		tmp <- strsplit(tax[i], "; ")[[1]][level]
		species[i] <- tmp
	}
}
da$taxonomy <-NULL
da <- data.frame(cbind(species, da), check.names=F)
da <- data.frame(ddply(da, "species", numcolwise(sum)), check.names=F)
rownames(da) <- da$species
da$species <- NULL
da <- data.frame(prop.table(as.matrix(da), 2), check.names=F)
da <- da*100
da$mean <- apply(da, 1, mean)
da <- da[,"mean", drop=F]
write.table(da, outputdata, sep="\t", row.names=T, col.names=T, quote=F)
