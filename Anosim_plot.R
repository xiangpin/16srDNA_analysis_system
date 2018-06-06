#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the Anosim.\\n Desinger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o","--outanosim", default="anosim_plot.svg", help="the output file of anosim, default is anosim_plot.svg.")
#myparser$add_argument("-o2","--outadonis", default="adonis_plot.svg", help="the output file of adonis, default is adonis_plot.svg.")


args <- myparser$parse_args()

library(vegan)
library(ggplot2)
mycolors <- c('#00AED7', '#FD9347', '#C1E168','#319F8C', "#984EA3", "#FF4040", "#228B22", "#FF7F00", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
#mycolors <- c('#00AED7', '#FD9347', '#C1E168','#319F8C', "#984EA3")
otu <- args$otu_tab
sampledat <- args$sample
out_anosim <- args$outanosim
#out_adonis <- args$outadonis

otutab <- read.table(otu, sep="\t", skip=1, header=T, row.names=1, check.names=F, comment.char = "")
#head(otutab)
otutab$taxonomy <- NULL

sampledat <- read.table(sampledat, header=T, sep="\t", check.names=F, comment.char = "")

otutab <- data.frame(t(otutab))
decorana(otutab)
otus_dist <- as.matrix(vegdist(otutab))
#attach(sampledat)
otuanosim <- anosim(otus_dist, sampledat[,2])
summary(otuanosim)
svg(out_anosim)
plot(otuanosim,col=mycolors,main="ANOSIM", ylab="Rank Dissimilarity")
dev.off()

#adonis_location <- adonis2(formula = otutab ~ group, data = sampledat)
#summary(adonis_location)
#svg(out_adonis)
#plot(adonis_location, col=mycolors, main="PermANOVA", ylab="Rank Dissimilarity")
#dev.off()

