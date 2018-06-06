#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to rare the otu by the depth.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.(have taxonomy column.)")
myparser$add_argument("-o","--output", default="otu_rare_norm.txt", help="the output file, default is otu_rare_norm.txt.")
myparser$add_argument('-F',"--flag", default="No", help="if the min depth of the total sample, set Yes, else is No, default is No.")
args <- myparser$parse_args()

library('vegan')
otu <- args$otu_tab
out<- args$output
tmpflag <- args$flag
da<-read.table(otu, header=T, skip=1, row.names=1, check.names=F, sep='\t', comment.char="", quote="")
samplelist <- as.vector(colnames(da))
tmptax <- data.frame(da$taxonomy)
colnames(tmptax) <- c("taxonomy")
rownames(tmptax) <- rownames(da)
if ('taxonomy' %in% samplelist ){
	da$taxonomy <- NULL
}else{
	da <- da
}
if (tmpflag == "No" ){
	data <- data.frame(t(decostand(data.frame(t(da),check.names=F), "hell")),check.names=F)
}else{
	data <- as.data.frame(t(da))
	set.seed(1000)
	newdata <- rrarefy(data, min(rowSums(data)))
	data <- data.frame(t(decostand(newdata, "hell")),check.names=F)
}
dat <- merge(data, tmptax, by=0)
#rownames(dat) <- dat$Row.names
colnames(dat)[1] <- "#OTU ID"
tmp <- data.frame("# Constructed from biom file")
write.table(tmp, out, sep="\t", col.names=F, row.names=F, quote=F)
write.table(dat, out, sep="\t", col.names=T, row.names=F, quote=F, append=T)
