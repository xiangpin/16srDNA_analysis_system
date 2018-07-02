#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the network with ggplot2 and igraph.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("nodestab", nargs=1, help="the nodes table file, the format such like\\n Feature\ttype\tmean\\n Feature1\tBacteria\t20.3")
myparser$add_argument("edgestab", nargs=1, help="the edges table file, the format such like\\n from\tto\tcorr\tpval\\n Feature1\tFeature3\t0.7\t0.02")
myparser$add_argument("-d", "--degree", default="network_degree_file.txt", help="the network degrees output, default is network_degree_file.txt.")
myparser$add_argument("-l", "--linestyle", default="curve", help="the style of network line, default is curve, you can chose curve or segment")
myparser$add_argument("-s", "--style", default=1, help="the style of network figure, default is 1, you can chose the 1 to 14")
myparser$add_argument("-g", "--gml", default="network_file_plot.gml", help="the gml file of igraph output, default is network_file_plot.gml.")
myparser$add_argument("-w", "--width", default=9, help="the width of output figures, default is 9.")
myparser$add_argument("-H", "--height", default=5, help="the height of output figures, default is 5.")
myparser$add_argument("-o", "--figure", default="network_figure_plot.svg", help="the output figures, default is network_figure_plot.svg.")

args <- myparser$parse_args()
suppressPackageStartupMessages(library(igraph))
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
set.seed(500)
nodes <- args$nodestab
link <- args$edgestab
degreefile <- args$degree
gmlfile <- args$gml
linestyle <- args$linestyle
style <- as.numeric(args$style)
wd <- as.numeric(args$width)
hd <- as.numeric(args$height)
figurefile <- args$figure
mycolors <- c(brewer.pal(n=9, name="Set1"), brewer.pal(n=8, name="Set2"))

nodes <- read.table(nodes, header=T, sep="\t", check.names=F)
links <- read.table(link, header=T, sep="\t", check.names=F)
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
netdeg <- degree(g, mode='all')
netdeg <- data.frame(netdeg)
netdeg <- cbind(rownames(netdeg), netdeg)
colnames(netdeg) <- c("Feature", "degree")
tmpnodes <- nodes
rownames(tmpnodes) <- tmpnodes$Feature
tmpnodes$Feature <- NULL
netdeg <- merge(netdeg, tmpnodes, by=0)
#colnames(netdeg)[1] <- "Feature"
netdeg$Row.names <- NULL
write.table(netdeg, degreefile, sep="\t", row.names=F, col.names=T, quote=F)
write.graph(g, gmlfile, format="gml")
if (style == 1){
	plotcord <- data.frame(layout.fruchterman.reingold(g))
}
if (style == 2){
	plotcord <- data.frame(layout.circle(g))
}
if (style == 3){
	plotcord <- data.frame(layout.kamada.kawai(g))
}
if (style == 4){
	plotcord <- data.frame(layout.davidson.harel(g))
}
if (style == 5){
	plotcord <- data.frame(layout.drl(g))
}
if (style == 6){
	plotcord <- data.frame(layout.gem(g))
}
if (style == 7){
	plotcord <- data.frame(layout.grid(g))
}
if (style == 8){
        plotcord <- data.frame(layout.lgl(g))
}
if (style == 9){
        plotcord <- data.frame(layout.mds(g))
}
if (style == 10){
        plotcord <- data.frame(layout.reingold.tilford(g))
}
if (style == 11){
        plotcord <- data.frame(layout.sphere(g))
}
if (style == 12){
        plotcord <- data.frame(layout.star(g))
}
if (style == 13){
        plotcord <- data.frame(layout.auto(g))
}
if (style == 14){
        plotcord <- data.frame(layout.graphopt(g))
}

rownames(plotcord) <- nodes$Feature
edgelist <- get.edgelist(g)
edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
colnames(edges) <- c("X1", "Y1", "X2", "Y2")
edges <- cbind(edges, links)
flagcor <- as.vector(edges$corr)
Flagcor <- vector()
for (i in 1:length(flagcor)){
	if (as.numeric(flagcor[i])>0){Flagcor[i] <- "positive"}
	if (as.numeric(flagcor[i])<0){Flagcor[i] <- "negative"}
}
edges <- cbind(edges,Flagcor)
rownames(nodes) <- nodes$Feature
plotcord <- merge(plotcord, nodes, by=0)
if (linestyle == "curve"){
	p <- ggplot()+
		geom_curve(aes(x=X1, y=Y1, xend = X2, yend =Y2, color=Flagcor, alpha=abs(corr)),
			size=0.6,
			data=edges)
}
if (linestyle == "segment"){
	p <- ggplot()+
		geom_segment(aes(x=X1, y=Y1, xend = X2, yend =Y2, color=Flagcor, alpha=abs(corr)),
			size=0.6,
			data=edges)
}
p <- p +
geom_point(data=plotcord, 
	aes(X1, X2, size=mean, shape=Kindom, fill=Phylum), 	
	color = "grey60",
	stroke = 0.8
	)+
geom_text_repel(data=plotcord, 
	aes(X1, X2, label=Feature), 
	family = "SimSun",
	fontface="bold",
	size=2,
	segment.size = 0,
	show.legend=F)+
scale_fill_manual(values=mycolors)+
#scale_fill_brewer(palette=mycolors)+
scale_shape_manual(values=c(21, 22, 24, 23, 25))+
theme_bw()+
labs(alpha="correlation size",
	color = "correlation type", 
	size="Abundance"
	#fill="Phylum"
	)+
guides(fill = guide_legend(override.aes = list(shape=21), order=3),
        size = guide_legend(override.aes = list(linetype = 0), order=4),
        alpha = guide_legend(override.aes = list(stroke = 0, shape=0), order=1),
	shape = guide_legend(order=5),
        color = guide_legend(override.aes = list(stroke = 0, shape=0, fill=0), order=2))+

theme(
	legend.text = element_text(size=7),
        legend.title = element_text(size=7),
	plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
	legend.key.height=unit(0.3,"line"),
        axis.ticks = element_blank()
              )
#guides(fill = guide_legend(override.aes =list(size=2)))


svg(figurefile, width=wd, height=hd)
p
dev.off()

