#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to caculate the relative abundance of different level tax.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("rdptax", nargs=1, help="the rdptax file of rdp.")
myparser$add_argument("-l", "--level", default="Phylum", help="the level tax, default is Phylum.(Class Order Family Genus)")
myparser$add_argument("-o", "--output", default="tab.txt", help="the output file, default is tab.txt.")

args <- myparser$parse_args()

library(plyr)
otu <- args$otu_tab
rdptab <- args$rdptax
level <- as.character(args$level)
out <- as.character(args$output)

da <- read.table(otu, header=T, row.names=1, check.names=FALSE, sep='\t', comment.char="", quote="")
rdptax <- read.table(rdptab, sep="\t", row.names=1, header=F, check.names=F, comment.char="")

rdptax$V2 <- NULL
rdptax$V4 <- NULL
rdptax$V5 <- NULL
rdptax$V7 <- NULL
rdptax$V8 <- NULL
rdptax$V10 <- NULL
rdptax$V11 <- NULL
rdptax$V13 <- NULL
rdptax$V14 <- NULL
rdptax$V16 <- NULL
rdptax$V17 <- NULL
rdptax$V19 <- NULL
rdptax$V20 <- NULL
colnames(rdptax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
datax <- data.frame(rdptax[[level]])
rownames(datax) <- rownames(rdptax)
colnames(datax) <- level
dat <- merge(da, datax, by=0)
dat$Row.names <- NULL

dt <- ddply(dat, level, numcolwise(sum))
rownames(dt) <- as.vector(dt[[level]])
dt[[level]] <- NULL
df <- data.frame(prop.table(as.matrix(dt), 2), check.names=F)*100
output <- paste(level, out, sep="_")

write.table(df, output, sep="\t", row.names=T, col.names=T, quote=F)
