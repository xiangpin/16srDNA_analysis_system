#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the Venn image.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of pca input.")
myparser$add_argument("-o","--output", default="venn_plot.svg", help="the output file, default is venn_plot.svg.")
args <- myparser$parse_args()

library(VennDiagram)
#mycolor <- c("#009ACD", "#66CD00", "#EEAD0E", "darkorchid1","#8B0000")
mycolor <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
otu <- args$otu_tab
sampledat <- args$sample
output <- args$output
da<-read.table(otu, header=T, skip=1, row.names=1, check.names=F, sep='\t', comment.char="")
#head(da, n=3)
sampledat <- read.table(sampledat, header=T, sep="\t", check.names=F, comment.char="")
#head(sampledat)
tmpspecie <- row.names(da)
samplelist <- as.vector(colnames(da))
if ('taxonomy' %in% samplelist ){
	da$taxonomy <- NULL
}else{
	da <- da
}

group <- as.vector(sampledat[, 2])
sample <- as.vector(sampledat[, 1])
#ample
c <- data.frame()
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	b <- da[[tmp]]
	if (i == 1){
		c <- data.frame(b, check.names=F)
	}
		c <- data.frame(cbind(c,b), check.names=F)
}
head(c)
a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- group
rownames(a) <- tmpspecie
a <- as.data.frame(t(apply(a, 1, function(x) tapply(x, colnames(a), mean))))
#head(a)
vennlist <- list()
for (i in 1:length(group)){
	tmp <- rownames(a[a[[group[i]]]>0,])
	tmpname <- group[i]
	vennlist[[tmpname]] <- tmp
}
#head(vennlist)

if (ncol(a) == 2){
	venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:2], alpha = 0.80, cat.col = mycolor[1:2], fontfamily = "serif", fontface = "bold", cex = 1.5, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = 0.05, margin = 0.1, lty ='dotted', lwd = 3, imagetype = "svg")
	#venn.diagram(vennlist, height=5, width=5,filename = output, fill = mycolor[1:2], alpha = 0.80, cat.col = mycolor[1:2], fontfamily = "serif", fontface = "bold", cex = 1.8, cat.cex = 1.8, cat.default.pos = "outer", cat.dist = 0.05, margin = 0.1, lty ='dotted', lwd = 3, imagetype = "png")
}



if (ncol(a) == 3){
        venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:3], alpha = 0.60, cat.col = mycolor[1:3], fontfamily = "serif", fontface = "bold", cex = 1.1, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = 0.08, margin = 0.1, lwd = 2, lty ='dotted',imagetype = "svg")
	#venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:3], alpha = 0.60, cat.col = mycolor[1:3], fontfamily = "serif", fontface = "bold", cex = 1.1, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = 0.08, margin = 0.1, lwd = 2, lty ='dotted',imagetype = "png")
}
if (ncol(a) == 4){
        venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:4], alpha = 0.85, cat.col = mycolor[1:4], fontfamily = "serif", fontface = "bold",cex = 1.2, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = c(0.22,0.22,0.12,0.12), margin = 0.1, lwd = 3, lty ='dotted', imagetype = "svg")
	#venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:4], alpha = 0.85, cat.col = mycolor[1:4], fontfamily = "serif", fontface = "bold",cex = 1.2, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = c(0.22,0.22,0.12,0.12), margin = 0.1, lwd = 3, lty ='dotted', imagetype = "png")
}
if (ncol(a) == 5){
        venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor[1:5], alpha = 0.50, cat.col = mycolor[1:5], fontfamily = "serif", fontface = "bold", cex = 0.5, cat.cex = 0.8, cat.default.pos = "outer", margin = 0.1, lwd = 0.6, lty ='dotted',imagetype = "svg")
}
if (ncol(a) == 6){
        venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor, cat.col = mycolor, alpha = 0.50, fontfamily = "serif", fontface = "bold", cex = 1.5, cat.cex = 1.8, cat.default.pos = "outer", margin = 0.1, lwd = 0.6, lty ='dotted', imagetype = "svg")
	#venn.diagram(vennlist, height=5, width=5, filename = output, fill = mycolor, cat.col = mycolor, alpha = 0.50, fontfamily = "serif", fontface = "bold", cex = 0.2, cat.cex = 0.3, cat.default.pos = "outer", cat.dist = 0.2, margin = 0.2, lwd = 0.4, lty ='dotted',imagetype = "png")
}

