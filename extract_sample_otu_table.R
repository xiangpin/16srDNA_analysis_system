#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to extract the otu table of sample which is in the sample file.\\n Desinger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.(have taxonomy column.)")
myparser$add_argument("sample", nargs=1, help="the sample file of pca input.")
myparser$add_argument("-o","--output", default="new_otu_table.txt", help="the output file, default is new_otu_table.txt.")
myparser$add_argument('-t',"--taxonomy", default="T",help="Whether contain the taxonomy, T or N, default is T.")
args <- myparser$parse_args()

otu <- args$otu_tab
groupfile <- args$sample
out<- args$output
flagtax <- args$taxonomy

da<-read.table(otu, header=T, skip=1, row.names=1, check.names=F, sep='\t', comment.char="", quote="")
flagtaxtmp <- colnames(da)[length(colnames(da))]
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
sample <- colnames(da)
#head(sample)
if (flagtax == "T"){
	samplelist<- as.vector(groupTable$sample)
	samplelist <- c(samplelist, flagtaxtmp)
}else{
	samplelist <- groupTable$sample
}
#print(names(da))
#print (samplelist)
#names <- names(da)[!(names(da) %in% samplelist)]
da <- da[,colnames(da) %in% samplelist]
#for(i in 1:length(names)){
#	remove <- as.character(names[i])
#	#print (remove)
#	sample <- setdiff(sample, remove)
#	da[[remove]] <- NULL
#}
head(da)
if (flagtax == "T"){
	taxtab <- as.data.frame(da[[flagtaxtmp]])
	rownames(taxtab) <- rownames(da)
	colnames(taxtab) <- flagtaxtmp
	da[[flagtaxtmp]] <- NULL 
}
da = da[rowSums(da)!=0,]
if (flagtax == "T"){
	newda <- merge(da, taxtab, by=0)
	rownames(newda) <- newda$Row.names
	#newda$Row.names <- NULL
	da <- newda
}
#head(da)
colnames(da)[1] <- "#OTU ID"
tmpt <- data.frame("# Constructed from biom file")
write.table(tmpt, out, sep="\t", col.names=F, row.names=F, quote=F)
write.table(da, out, sep="\t", col.names=T, row.names=F, quote=F, append=T)
