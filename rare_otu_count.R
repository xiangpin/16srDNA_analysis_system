#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to rare the otu by the depth.\\n Desinger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.(have taxonomy column.)")
myparser$add_argument("-o","--output", default="new_otu_table.txt", help="the output file, default is new_otu_table.txt.")
args <- myparser$parse_args()

library('vegan')
otu <- args$otu_tab
out<- args$output
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
data <- as.data.frame(t(da))
set.seed(1000)
print (min(rowSums(data)))
newdata <- rrarefy(data, min(rowSums(data)))
#newdata <- decostand(newdata, "hell")
#head(newdata)
dat <- data.frame(t(newdata),check.names=F)
dat <- merge(dat, tmptax, by=0)
#rownames(dat) <- dat$Row.names
#dat$Row.names <- NULL

colnames(dat)[1] <- "#OTU ID"
tmpt <- data.frame("# Constructed from biom file")
write.table(tmpt, out, sep="\t", col.names=F, row.names=F, quote=F)
write.table(dat, out, sep="\t", col.names=T, row.names=F, quote=F, append=T)
