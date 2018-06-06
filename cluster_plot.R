#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the cluster.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-m", "--hcmethod", default="average", help="the method of cluster, default is average(UPGMA), optional complete,single,median,centroid,mcquitty,ward.D2.")
myparser$add_argument("-d", "--dsmethod", default="euclidean", help="the method of distance, default is euclidean, default is euclidean,optional binary,manhattan,maximum.")
myparser$add_argument("-l", "--level", default=8, help="the tax levels to want to plot, default is 8, total, optional 2-8.")
myparser$add_argument("-s", "--labelsize", default=1, help="the size of label cluster, default is 1.")
myparser$add_argument("-o","--output", default="cluster_plot.svg", help="the output file, default is cluster_plot.svg.")

args <- myparser$parse_args()
suppressPackageStartupMessages(library(dendextend))
library(colorspace)
library(plyr)
otu <- args$otu_tab
groupfile <- args$sample
dmethod <- args$dsmethod
hcmethod <- args$hcmethod
size <- as.numeric(args$labelsize)
level <- as.numeric(args$level)
outfigure <- args$output

da <- read.table(otu, header=T, skip=1, row.names=1, check.names=FALSE, sep='\t', comment.char="", quote="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
tax <- as.vector(da$taxonomy)
species <- vector()
for (i in 1:length(tax)){
        if (level==8){
                tmp <- tax[i]
                species[i] <- tmp
        }
        else {
                tmp <- strsplit(tax[i], "; ")[[1]][level]
                species[i] <- tmp
        }
}
da$taxonomy <- NULL
da <- data.frame(cbind(species, da), check.names=F)
da <- data.frame(ddply(da, "species", numcolwise(sum)), check.names=F)
rownames(da) <- da$species
da$species <- NULL
da <- data.frame(prop.table(as.matrix(da), 2), check.names=F)
da <- da*100
tmpspecies <- row.names(da)
sample <- as.vector(groupTable[,1])
group <- as.vector(groupTable[,2]) 
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	b <- da[[tmp]]
	c<- data.frame(cbind(c,b))#, check.names=F)
	
}

a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- as.vector(groupTable[,1])
rownames(a) <- tmpspecies
#head(a)
dd <- data.frame(t(a), check.names=F)

dg <- cbind(dd, group)
tmpgroup <- dg[,colnames(dg)=="group", drop=F]
sample_d <- dist(dd, method = dmethod)
hc_sample <- hclust(sample_d, method = hcmethod)
group <- unique(as.vector(dg$group))
#print(group)
colors <- rainbow_hcl(length(group))
tt <- data.frame(cbind(group, colors))
#print(tt)
#group <- rev(levels(dg$group))
#print(colors)
dend <- as.dendrogram(hc_sample)
dend <- rotate(dend, 1:length(sample))
#print(labels_colors(dend))
tmppgroup <- as.vector(dg[order.dendrogram(dend),]$group)
labels_colors(dend) <- vector()
for (i in 1:length(tmppgroup)){
	labels_colors(dend)[i] <- as.character(tt[tt[,1]==tmppgroup[i],2])
}
#print(labels_colors(dend))
#print(colors)
#print(sort_levels_values(as.numeric(dg$group)))
#print(order.dendrogram(dend))
#labels_colors(dend) <- colors[sort_levels_values(as.numeric(dg$group)[order.dendrogram(dend)])]
print(labels_colors(dend))
tmpcolor <- data.frame(labels_colors(dend), check.names=F)
#print(tmpcolor)
colnames(tmpcolor) <- "colors"
tmpgroupcolor <- merge(tmpgroup, tmpcolor, by=0)
#print(tmpgroupcolor)
tmpcolor <- data.frame(t(tmpcolor), check.names=F)
#print (tmpcolor)
groupcolors <- vector()
for (i in 1:length(group)){
	tmp <- as.character(groupTable[groupTable[,2]==group[i], 1][1])
	#print (tmp)
	groupcolors[i] <- as.character(tmpcolor[[tmp]])
	
}
#print (groupcolors)
dend <- hang.dendrogram(dend, hang = -1)
dend <- set(dend, "labels_cex", size)
svg(outfigure)
par(mar = c(2,3,0.2,7))
plot(dend, horiz =  TRUE, nodePar = list(cex = .007))
legend("topleft", legend = group, fill = groupcolors, bty = "n")
dev.off()
