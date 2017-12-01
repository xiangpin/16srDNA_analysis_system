suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the cluster.\\n Desinger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-m", "--hcmethod", default="average", help="the method of cluster, default is average(UPGMA).")
myparser$add_argument("-d", "--dsmethod", default="binary", help="the method of distance, default is binary.")
myparser$add_argument("-s", "--labelsize", default=1, help="the size of label cluster, default is 1.")
myparser$add_argument("-o","--output", default="cluster_plot.svg", help="the output file, default is cluster_plot.svg.")

args <- myparser$parse_args()
#library(dendextend)
suppressPackageStartupMessages(library(dendextend))
library(colorspace)
otu <- args$otu_tab
groupfile <- args$sample
dmethod <- args$dsmethod
hcmethod <- args$hcmethod
size <- args$labelsize
outfigure <- args$output
da <- read.table(otu, header=T, row.names=1, check.names=FALSE, sep='\t', comment.char="", quote="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
da$taxonomy <- NULL
tmpspecies <- row.names(da)
sample <- as.vector(groupTable[,1])
group <- as.vector(groupTable[,2]) 
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	b <- da[[tmp]]
	c<- data.frame(cbind(c,b), check.names=F)
}

a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- as.vector(as.character(groupTable[,1]))
rownames(a) <- tmpspecies
dd <- data.frame(t(a), check.names=F)
dg <- cbind(dd, group)

sample_d <- dist(dd, method = dmethod)
hc_sample <- hclust(sample_d, method = hcmethod)
colors <- rainbow_hcl(length(levels(dg$group)))
#attributes(hc_sample)
#hc_sample$order
#hc_sample$labels
group <- rev(levels(dg$group))
dend <- as.dendrogram(hc_sample)
dend <- rotate(dend, 1:length(sample))
#dend <- color_branches(dend, k = length(group), groupLabels=group)
labels_colors(dend) <- colors[sort_levels_values(as.numeric(dg$group)[order.dendrogram(dend)])]
tmpcolor <- data.frame(labels_colors(dend))
colnames(tmpcolor) <- "colors"
tmpcolor <- data.frame(t(tmpcolor), check.names=T)
groupcolors <- vector()
for (i in 1:length(group)){
	tmp <- as.character(groupTable[groupTable[,2]==group[i], 1][1])
	groupcolors[i] <- as.vector(tmpcolor[[tmp]])
	
}
dend <- hang.dendrogram(dend, hang = -1)
dend <- set(dend, "labels_cex", size)
svg(outfigure)
plot(dend, horiz =  TRUE, nodePar = list(cex = .007))
legend("topleft", legend = group, fill = groupcolors, border="white", bty = "n")
dev.off()
