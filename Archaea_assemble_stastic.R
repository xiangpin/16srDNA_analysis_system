#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to merge the Raw and Clean reads tab.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("rawtab", nargs=1, help="the raw_nums_tab.")
myparser$add_argument("cleantab", nargs=1, help="the clean_nums_tab.")
myparser$add_argument("checktab", nargs=1, help="the Chimera otu tab.")
myparser$add_argument("otutab", nargs=1, help="the usearch global alignment otu tab.")
myparser$add_argument("Archaeatab", nargs=1, help="the Archaea otu tab.")
myparser$add_argument("-o","--output", default="Statistic_RawCleanNums.xls", help="the output file, default is Statistic_RawCleanNums.xls.")

args <- myparser$parse_args()

rawda <- args$rawtab
cleanda <- args$cleantab
chimera <- args$checktab
otutab <- args$otutab
outfile <- args$output
Archaeatab <- args$Archaeatab
rawda <- read.table(rawda, sep="\t", row.names=1, header=F)
cleanda <- read.table(cleanda, sep="\t", row.names=1, header=F)
chimerada <- read.table(chimera, sep="\t", header=T, row.names=1, check.names=F, comment.char="")
otutabda <- read.table(otutab, sep="\t", header=T, row.names=1, check.names=F, comment.char="")

Archaeatab <- read.table(Archaeatab, sep="\t", header=T, row.names=1, check.names=F, skip=1, comment.char="")
Archaeatab$taxonomy <- NULL

j = 1
otu <- data.frame()
for (i in (1:length(colnames(otutabda)))){
	d <- sum(otutabda[,i]>0)
	otu[j, 1] <- colnames(otutabda)[i]
	otu[j,2] <- d
	j <- j+1
}
colnames(otu) <- c("sample","OTUs")
rownames(otu) <- otu$sample
otu$sample <- NULL

z <- 1
ArchaeaOTUs <- data.frame()
for (i in (1:length(colnames(Archaeatab)))){
        d <- sum(Archaeatab[,i]>0)
        ArchaeaOTUs[z, 1] <- colnames(Archaeatab)[i]
        ArchaeaOTUs[z, 2] <- d
        z <- z+1
}
colnames(ArchaeaOTUs) <- c("sample", "ArchaeaOTUs")
rownames(ArchaeaOTUs) <- ArchaeaOTUs$sample
ArchaeaOTUs$sample <- NULL
ArchaeaOTUseqs <- data.frame(colSums(Archaeatab))
colnames(ArchaeaOTUseqs) <- "ArchaeaOTUseqs"
head(ArchaeaOTUseqs)
dt <- merge(rawda, cleanda, by=0)
rownames(dt) <- dt$Row.names
dt$Row.names <- NULL
colnames(dt) <- c("Raw", "Clean")
chimeranum <- data.frame(colSums(chimerada))
colnames(chimeranum) <- "Chimerta_check"
dt <- merge(dt, chimeranum, by=0)
rownames(dt) <- dt$Row.names
dt$Row.names <- NULL
OTUseqs <- data.frame(colSums(otutabda))
colnames(OTUseqs) <- "OTUseqs"
dt <- merge(dt, OTUseqs, by=0)
rownames(dt) <- dt$Row.names
dt$Row.names <- NULL
dt <- merge(dt, ArchaeaOTUseqs, by=0)
rownames(dt) <- dt$Row.names
dt$Row.names <- NULL
#dt <- merge(dt, ArchaeaOTUs, by=0)
#rownames(dt) <- dt$Row.names
#dt$Row.names <- NULL
dt <- merge(dt, otu, by=0)
rownames(dt) <- dt$Row.names
dt$Row.names <- NULL
dt <- merge(dt, ArchaeaOTUs, by=0)
#rownames(dt) <- dt$Row.names

colnames(dt)[1] <- "sample"
write.table(dt, outfile, sep="\t", row.names=F, col.names=T, quote=F)

