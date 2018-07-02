#!/usr/bin/Rscript
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

sampleda <- read.table(input, header=T, sep="\t", comment.char="", check.names=F)
head(sampleda)
dd <- data.frame()
print (groupids)
for (i in 1:length(groupids)){
	tmp <- sampleda[sampleda[,2]==groupids[i],]
	#print(sampleda)
	dd <- rbind(dd, tmp)

}
#head(dd)
if (ncol(sampleda) == 3){
	dd[,2]<-NULL
}else{
	dd[,3] <- NULL
}
#head(dd)
colnames(dd) <- c("sample", "group")
write.table(dd, output, sep="\t", row.names=F, col.names=T, quote=F)
