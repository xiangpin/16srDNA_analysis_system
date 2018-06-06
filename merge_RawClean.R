#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to merge the Raw and Clean reads tab.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("rawtab", nargs=1, help="the raw_nums_tab.")
myparser$add_argument("cleantab", nargs=1, help="the clean_nums_tab.")
#myparser$add_argument("-s", "--pointsize", default=2, help="the point size, default is 2.")
myparser$add_argument("-o","--output", default="merge_RawCleanNums.xls", help="the output file, default is merge_RawCleanNums.xls.")

args <- myparser$parse_args()

rawda <- args$rawtab
cleanda <- args$cleantab
outfile <- args$output
rawda <- read.table(rawda, sep="\t", row.names=1, header=F)
cleanda <- read.table(cleanda, sep="\t", row.names=1, header=F)

dt <- merge(rawda, cleanda, by=0)
colnames(dt) <- c("sample", "RawReadsNums", "CleanReadsNums")
write.table(dt, outfile, sep="\t", row.names=F, col.names=T, quote=F)
