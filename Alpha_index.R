#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the boxplot of alpha.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
#myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
#myparser$add_argument("-s","--size", default=10000, help="the size of the rare, default is 10000.")
#myparser$add_argument("-m","--method", default="Observed", help="the rare eveness, default is Observed, Choose (Chao1, se.chao1, Shannon, InvSimpson).")
#myparser$add_argument("-o","--output", default="alpha_box_plot.svg", help="the output file, default is alpha_box_plot.svg.")
myparser$add_argument("-d","--alphatab", default="alpha_tab.txt", help="the output file, default is alpha_tab.txt.")

args <- myparser$parse_args()

library(vegan)
library(ggplot2)
library(reshape2)
library(agricolae)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm),
          Q95  = quantile(xx[[col]], 0.95, na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

otu <- args$otu_tab
groupfile <- args$sample
outfigure <- args$output
alphatab <- args$alphatab
da <- read.table(otu, header=T, skip=1,row.names=1, check.names=FALSE, sep='\t', comment.char="")
da$taxonomy <- NULL
da <- data.frame(t(da), check.names=F)
#da <- da[rowSums(da) >= 30000,]
set.seed(1000)
print(min(rowSums(da)))
data <- rrarefy(da, min(rowSums(da)))
#head(da)
#data <- da
Chao <- estimateR(data)
Shannon <- diversity(data)
Simpson <- diversity(data, index="simpson")
J <- Shannon/log(specnumber(data))

alpha <- data.frame(Observe=Chao[1,], Chao1=Chao[2,], ACE=Chao[4,], Shannon, Simpson, J, check.names=F)
write.table(alpha, alphatab, sep="\t", col.names=T, row.names=T, quote=F)
