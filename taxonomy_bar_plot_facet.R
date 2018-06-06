#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the bar of species.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("OTUtab", nargs=1, help="the otu tax table file.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-l", "--level", default=2, help="the tax levels to want to plot, default is 2, Phylum.")
myparser$add_argument("-o","--output", default="specie_bar_plot.svg", help="the output file, default is specie_bar_plot.svg.")

args <- myparser$parse_args()

mycolors <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")

library(ggplot2)
library(plyr)
library(reshape2)
otu <- args$OTUtab
groupfile <- args$sample
outfigure <- args$output
#outputdata <- args$data

level <- as.integer(args$level)
da<-read.table(otu, header=T, skip=1, row.names=1,check.names=FALSE, sep='\t', comment.char="", quote="")
groupTable <- read.table(groupfile, header=T, sep="\t", check.names=F, comment.char="")
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
da$taxonomy <-NULL
da <- data.frame(cbind(species, da), check.names=F)
da <- data.frame(ddply(da, "species", numcolwise(sum)), check.names=F)
rownames(da) <- da$species
da$species <- NULL
da <- data.frame(prop.table(as.matrix(da), 2), check.names=F)
da <- da*100
#write.table(da, outputdata, sep="\t", row.names=T, col.names=T, quote=F)
tmpspecies <- row.names(da)
sample<- as.vector(groupTable[,1])
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	b <- da[[tmp]]
	c<- data.frame(cbind(c,b), check.names=F)
}
a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- as.vector(groupTable[,1])
rownames(a) <- tmpspecies
a$sum <- apply(a, 1, sum)
if ((length(a$sum)) < 29){
        d<- data.frame(a[rev(order(a$sum)),][1:length(a$sum),],check.names=F)
} else
{
        d <- data.frame(a[rev(order(a$sum)),][1:29,], check.names=F)
}

d<-data.frame(d[colnames(d)!="sum"], check.names=F)
e<-data.frame(matrix(apply(d,2,sum), nrow=1))
colnames(e) <- colnames(d)
rownames(e) <- c("Others")
e <- 100 - e
d <- data.frame(rbind(d,e), check.names=F)
factorLevel <- rownames(d)
d$species <- rownames(d)
b<-melt(d, id="species", variable.name="sample")
b$species <- factor(b$species, levels=factorLevel)
grouplist <- vector()
samplelist <- as.vector(b$sample)
for (i in 1:length(samplelist)){
	grouplist[i] <- as.character(groupTable[groupTable[,1]==samplelist[i],2])
}
b$group <- grouplist
p <- ggplot(data=b) +
geom_bar(aes(y=value, x=sample, fill=species), 
position = "stack", 
width=0.7, 
stat="identity")+
#, fun.y = "mean") + #, stat="identity") +
#scale_x_discrete("group", labels = as.vector(groupTable$group)) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), 
#axis.text.x = element_blank(),
panel.grid = element_blank(), 
panel.border = element_blank(),
axis.line = element_line(size=0.5, colour = "black"), 
#axis.line.y = element_line(size=0.5, colour = "black")
#axis.ticks.x = element_blank()
) +
facet_wrap( ~ group, scale="free_x", nrow=1)+
xlab("") +
ylab("relative abundance (%)") +
scale_fill_manual(values=mycolors) +
guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
theme(legend.position = "bottom", 
legend.box = "horizontal", 
legend.text = element_text(size = 8), 
legend.title=element_blank(),
panel.spacing = unit(0.2, "mm"),
strip.text.x = element_text(size=12, face="bold"),
strip.background = element_rect(colour="white", fill="grey")) +
labs(fill="") + 
scale_y_continuous(expand=c(0,0))
svg(outfigure, height=5, width=8)
p
dev.off()

