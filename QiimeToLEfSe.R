#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to convert qiime to LEfSe.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-d","--tmpdata", default="LEfSe_sample_file.txt", help="the output file (header is sample ID), default is LEfSe_sample_file.txt.")
myparser$add_argument("-o","--output", default="LEfSe_input_file.txt", help="the output file (header is group ID), default is LEfSe_input_file.txt.")

args <- myparser$parse_args()

library(ggplot2)
library(plyr)
library(reshape2)
otu <- args$otu_tab
groupfile <- args$sample
tmpotu <- args$tmpdata
outfile <- args$output

da <- read.table(otu, sep="\t", skip=1, row.names=1, header=T, check.names=F, comment.char="")
sample <- read.table(groupfile, sep="\t", header=T, row.names=1, check.names=F, comment.char="")
tax <- as.vector(da$taxonomy)
da$taxonomy <- NULL
ktab <- vector()
ptab <- vector()
ctab <- vector()
otab <- vector()
ftab <- vector()
gtab <- vector()
stab <- vector()
for (i in 1:length(tax)){
	ktab[i] <- strsplit(tax[i], "; ")[[1]][1]
	k <- strsplit(tax[i], "; ")[[1]][1]
	p <- strsplit(tax[i], "; ")[[1]][2]
	c <- strsplit(tax[i], "; ")[[1]][3]
	o <- strsplit(tax[i], "; ")[[1]][4]
	f <- strsplit(tax[i], "; ")[[1]][5]
	g <- strsplit(tax[i], "; ")[[1]][6]
	s <- strsplit(tax[i], "; ")[[1]][7]
	ptab[i] <- paste(k, p, sep="|")
	ctab[i] <- paste(k, p, c, sep="|")
	otab[i] <- paste(k, p, c, o, sep="|")
	ftab[i] <- paste(k, p, c, o, f, sep="|")
	gtab[i] <- paste(k, p, c, o, f, g, sep="|")
	stab[i] <- paste(k, p, c, o, f, g, s, sep="|")
}
taxratio <- function(da, taxlist){
	data <- data.frame(cbind(taxlist, da), check.names=F)
	colnames(data)[1] <- "tax"
	data <- data.frame(ddply(data, "tax", numcolwise(sum)), check.names=F)
	rownames(data) <- data$tax
	data$tax <- NULL
	data <- data.frame(prop.table(as.matrix(data), 2), check.names=F)
	data <- data*100
	return (data)
}
kda <- taxratio(da, ktab)
pda <- taxratio(da, ptab)
cda <- taxratio(da, ctab)
oda <- taxratio(da, otab)
fda <- taxratio(da, ftab)
gda <- taxratio(da, gtab)
sda <- taxratio(da, stab)
data <- rbind(kda, pda, cda, oda, fda, gda, sda)
data <- data[(sort(rownames(data))),]
data <- data[rowSums(data)>0,]
#write.table(data, tmpotu, sep="\t", row.names=T, col.names=T, quote=F)
dat <- data.frame(t(data), check.names=F)

data <- merge(sample, dat, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL
data <- data.frame(t(data), check.names=F)
write.table(data, tmpotu, sep="\t", row.names=T, col.names=T, quote=F)
write.table(data, outfile, sep="\t", row.names=T,col.names=F, quote=F)
