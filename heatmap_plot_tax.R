#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the heatmap of genus by group.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com."
parser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
	parser$add_argument("data", nargs=1, help="the species file.")
	parser$add_argument("sample", nargs=1, help="the sample file.")
	parser$add_argument("-t", "--top", default=30, help="the top number.")
	parser$add_argument("-w", "--width", default=12, help="the cell width of the heatmap, default is 12.")	
	parser$add_argument("-l", "--level", default=8, help="the tax levels to want to plot, default is 8, total.")
	parser$add_argument("-H", "--height", default=16, help="the cell height of the heatmap, default is 16.")
#	parser$add_argument("-H", "--height", default=6, help="the cell height of the heatmap, default is 6.")	
	parser$add_argument("-T", "--total", default="T", help="if use the total circRNA or others, default is T (yes).")
	parser$add_argument("-F", "--cluster", default="T", help="if plot the cluster sample. default is T(yes).")
	parser$add_argument("-o", "--output", default="heatmap_plot.svg", help="the Heatmap plot, default is heatmap_plot.svg.")
	args <- parser$parse_args()
library(ggplot2)
library(pheatmap)
library(dplyr)
library(recommenderlab)
library(Heatplus)
library(RColorBrewer)
library(gplots)
library(plyr)
dataFile <- args$data
groupFile <- args$sample
top <- as.integer(args$top)
total <- as.character(args$total)
cluster <- as.character(args$cluster)
level <- as.numeric(args$level)
output <- args$output
cell_width <- as.numeric(args$width)
cell_height <- as.numeric(args$height)


colors <- c("blue","white","red")
da <- read.table(dataFile, header=T, skip=1, sep = '\t', row.names=1, check.names=F, comment.char="")
grouptable <- read.table(groupFile, header=T, row.names=1, sep="\t", check.names=F)
groupTable <- data.frame(grouptable[,1])
rownames(groupTable) <- rownames(grouptable)
colnames(groupTable) <- "group"
groupnum <- nrow(data.frame(table(groupTable[,1])))
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
sample <- vector()
sample<- rownames(groupTable)#$sample
for (i in 1:length(sample)){
        tmp <- as.character(sample[i])
        b <- da[[tmp]]
        c<- data.frame(cbind(c,b))
}

a<-c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- sample
rownames(a) <- tmpspecies
a <- a[rowSums(a) > 0, ]
a$sum <- apply(a, 1, sum)
if (total == "T"){
	d<-a[rev(order(a$sum)),][1:length(a$sum),]
	}else{
		if ((length(a$sum)) < top){
        		d<-a[rev(order(a$sum)),][1:length(a$sum),]
		} else
		{
        		d<-a[rev(order(a$sum)),][1:top,]
		}
}
#anno_color <- c('#00AED7','#C1E168','#FD9347','#319F8C',"#F8766E","#7CAD00", "#00BEC3", "#C67CFF","#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
anno_color <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
group <- anno_color[1:groupnum]
names(group) <- levels(groupTable[,1])

ann_colors = list(group = group)
d<-d[colnames(d)!="sum"]
d <- data.frame(d, check.names=F)
#d
cellheight <- cell_height
cellwidth <- cell_width
if (cluster=="T"){
	svg(output, width=8, height=8)
	pheatmap(log2(d+1), scale="row", col=colorRampPalette(colors,space="rgb")(5000), border_color='#FFFAEA', cluster_rows=T, cluster_cols=T, show_rownames=T, show_colnames=F, cellheight=cellheight, cellwidth=cellwidth, annotation_colors=ann_colors, annotation_col=groupTable, annotation_names_col=F)
	dev.off()
}else{
	svg(output, width=8, height=8)
	pheatmap(log2(d+1), scale="row", col=colorRampPalette(colors,space="rgb")(5000), border_color='#FFFAEA', cluster_rows=T, cluster_cols=F, show_rownames=T, show_colnames=F, cellheight=cellheight, cellwidth=cellwidth, annotation_colors=ann_colors, annotation_col=groupTable, annotation_names_col=F)
	dev.off()
}

