suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the cluster and bar of tax.\\n Developger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-m", "--hcmethod", default="average", help="the method of cluster, default is average(UPGMA), optional complete,single,median,centroid,mcquitty,ward.D2.")

myparser$add_argument("-d", "--dsmethod", default="euclidean", help="the method of distance, default is euclidean,optional binary,manhattan,maximum.")
myparser$add_argument("-n", "--topnum", default=30, help="the top number tax, default is 30.")
myparser$add_argument("-a", "--scalehight", default=0.01, help="the scale height of the cluster tree, default is 0.01")
myparser$add_argument("-z", "--magin", default=4, help="the blank size of the two plot, default is 4.")
myparser$add_argument("-s", "--labelsize", default=7, help="the size of label cluster, default is 7.")
myparser$add_argument("-o","--output", default="cluster_barplot.svg", help="the output file, default is cluster_barplot.svg.")

args <- myparser$parse_args()
library(ggdendro)
library(ggplot2)
library(reshape2)
#library(egg)
library(gridExtra)
library(colorspace)
mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C',"#F8766E","#7CAD00","#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")

otu <- args$otu_tab
groupfile <- args$sample
dmethod <- args$dsmethod
hcmethod <- args$hcmethod
topnum <- as.numeric(args$topnum)
scalehight <- as.numeric(args$scalehight)
size <- args$labelsize
outfigure <- args$output
magin <- as.numeric(args$magin)
da <- read.table(otu, header=T, row.names=1, check.names=FALSE, sep='\t', comment.char="", quote="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
da$taxonomy <- NULL
tmpspecies <- row.names(da)
sample <- as.vector(groupTable[,1])
group <- unique(as.vector(groupTable[,2]))
groupcolors <- rainbow_hcl(length(group))
group2colors <- data.frame(cbind(group, groupcolors))
dd <- data.frame(t(da), check.names=F)
#head(dd)
model <- hclust(dist(dd, dmethod), hcmethod)
dhc <- as.dendrogram(model)
ddata <- dendro_data(dhc, type = "rectangle")
sample <- ddata$labels$label
samplecolors <- vector()
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	grouptmp <- as.character(groupTable[groupTable[,1]==tmp, 2])
	samplecolors[i] <- as.character(group2colors[group2colors[,1]==grouptmp, 2][1])
	b <- da[[tmp]]
	c <- data.frame(cbind(c,b))
}

a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- sample
rownames(a) <- tmpspecies

a$sum <- apply(a, 1, sum)
if ((length(a$sum)) < topnum){
        d<-a[rev(order(a$sum)),][1:length(a$sum),]
} else{
        d<-a[rev(order(a$sum)),][1:topnum,]
}

d<-d[colnames(d)!="sum"]
factorLevel <- rownames(d)
d$species <- rownames(d)
b <- melt(d, id="species")
b$species <- factor(b$species, levels=factorLevel)

blank_theme <- theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid = element_blank(), panel.border = element_blank())

p_cl <- ggplot(segment(ddata))+
geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
scale_y_reverse(expand = c(0, 0))+
scale_x_continuous(expand = c(scalehight,0)) +
theme_bw()+
xlab("") +
ylab("") +
blank_theme +
theme(legend.position = "", axis.text.y = element_blank(), plot.margin=unit(c(3,0,0,0), "mm")) +
coord_flip()

p_bar <- ggplot(data=b) +
geom_bar(aes(y=value, x=variable, fill=species), position = "fill", width=0.45, stat="identity") +
theme_bw() +
xlab("") + 
blank_theme +
ylab("") +
#scale_x_continuous(expand = c(0.01, 0)) +
scale_fill_manual(values=mycolors)+
guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8, ncol=1)) +
theme(axis.text.y = element_text(colour = samplecolors, size=size), plot.margin=unit(c(3, 0, 0, -magin), "mm")) +
theme(legend.position = "right", legend.text = element_text(size = 7.5))+
labs(fill="") +
scale_y_continuous(expand=c(0,0)) +
coord_flip()

svg(outfigure)
grid.arrange(p_cl, p_bar, widths=c(4,13))
dev.off()

