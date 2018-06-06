#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to run the indicator.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-l", "--level", default=6, help="the tax levels to want to plot, default is 6, Genus.")
myparser$add_argument("-f", "--foldfdr", default=0.05, help="the fold value of p.value, default is 0.05.")
myparser$add_argument("-s", "--xsize", default=8.5, help="the size of the bar axis.text.y, default is 8.5.")
myparser$add_argument("-b", "--barsize", default=0.7, help="the bar width, default is 0.7.")
myparser$add_argument("-w", "--width", default=12, help="the outfigure width, default is 12.")
myparser$add_argument("-H", "--height", default=14, help="the outfigure height, default is 14.")
myparser$add_argument("-m", "--methodpadj", default="qvalue", help="the adjust method of p.value, default is p.adjust, another is qvalue.")
myparser$add_argument("-p", "--pfilter", default="No", help="whether filter the pvalue or FDR, default is No(FDR), Yes (pvalue).")
myparser$add_argument("-T", "--filterabundance", default=F, help="whether filter the low abundance taxon, default is F.")
myparser$add_argument("-a", "--foldabundance", default=0.01, help="the fold value of low abundance taxon, default is 0.01.")
myparser$add_argument("-d", "--data", default="indicator_output_data.txt", help="the output data, default is indicator_output_data.txt.")
myparser$add_argument("-o","--output", default="indicator_bar_plot.svg", help="the output figure file, default is indicator_bar_plot.svg.")

args <- myparser$parse_args()

mycolors <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")

library(ggplot2)
library(plyr)
library(reshape2)
library(indicspecies)
library(qvalue)
library(egg)
otu <- args$otu_tab
groupfile <- args$sample
foldFDR <- as.numeric(args$foldfdr)
barsize <- as.numeric(args$barsize)
xsize <- as.numeric(args$xsize)
width <- as.numeric(args$width)
height <- as.numeric(args$height)
filterabundance <- args$filterabundance
foldabundance <- as.numeric(args$foldabundance)
method <- as.character(args$methodpadj)
outfigure <- args$output
outputdata <- args$data
pfilter <- as.character(args$pfilter)
set.seed(1000)
level <- as.numeric(args$level)

da<-read.table(otu, header=T, skip=1, row.names=1,check.names=FALSE, sep='\t', comment.char="")
groupTable <- read.table(groupfile, header=T, sep="\t", check.names=F, comment.char="")
tax <- as.vector(da$taxonomy)
species <- vector()
phylum <- vector()
for (i in 1:length(tax)){
	if (level==8){
		tmp <- tax[i]
		species[i] <- tmp
	}
	else {
		tmp <- strsplit(tax[i], "; ")[[1]][level]
		phylum[i] <- strsplit(tax[i], "; ")[[1]][2]
		species[i] <- tmp
		
	}
}

taxdat <- data.frame(cbind(phylum, species))
uniqphylum <- vector()
for (i in 1:length(unique(taxdat$species))){
	tmp <- as.character(unique(taxdat$species)[i])
	uniqphylum[i] <- as.character(taxdat[taxdat[,2]==tmp,1][1])
}
newtaxdat <- data.frame(uniqphylum, unique(taxdat$species))
colnames(newtaxdat) <- c("phylum", "species")
rownames(newtaxdat) <- newtaxdat$species
da$taxonomy <-NULL
da <- data.frame(cbind(species, da), check.names=F)
da <- data.frame(ddply(da, "species", numcolwise(sum)), check.names=F)
rownames(da) <- da$species
da$species <- NULL
da <- data.frame(prop.table(as.matrix(da), 2), check.names=F)
da <- da*100
tmpspecies <- row.names(da)
sample<- as.vector(groupTable[,1])
for (i in 1:length(sample)){
	tmp <- as.character(sample[i])
	b <- da[[tmp]]
	c<- data.frame(cbind(c,b), check.names=F)
}
a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- as.vector(as.character(groupTable[,1]))
rownames(a) <- tmpspecies
da <- a
mean_abundance <- apply(da, 1, mean)
tmpab <- data.frame(mean_abundance)
rownames(tmpab) <- rownames(da)
dat <- data.frame(t(da), check.names=F)
group <- as.vector(groupTable[,2])
#phi <- multipatt(da, group, func = "r.g", control= how(nperm=999))
indval <- multipatt(dat, group, func = "indval.g", control= how(nperm=999))
strdata <- data.frame(indval$str, check.names=F)
strdata$p.value <- indval$sign$p.value
if (method=="p.adjust"){
	strdata$FDR <- p.adjust(strdata$p.value, method = "BH")
}else{
	strdata$FDR <- qvalue(strdata$p.value)$qvalues
}
if (pfilter=="Yes") {
	strdat <- strdata[strdata$p.value<=foldFDR&!is.na(strdata$p.value), ]
}else{
	strdat <- strdata[strdata$FDR<=foldFDR&!is.na(strdata$FDR), ]
}
#strdat <- strdata[strdata$FDR<=foldFDR&!is.na(strdata$FDR), ]
strdatab <- merge(tmpab, strdat, by=0)

rownames(strdatab) <- strdatab$Row.names
strdatab$Row.names <- NULL

newdata <- merge(newtaxdat, strdatab, by=0)
newdata$Row.names <- NULL

tmpphylum <- unique(newdata$phylum)
dat <- data.frame()
for (i in 1:length(tmpphylum)){
	tmpstr <- as.character(tmpphylum[i])
	tmp <-newdata[newdata$phylum==tmpstr,]
	dat <- rbind(dat, tmp)
}
if (filterabundance=="T"){
        dat <- dat[dat$mean_abundance >= foldabundance, ]
}else{
        dat <- dat
}

factorlevel <- dat$species
dat$species <- factor(dat$species, levels=factorlevel)
pointdat <- dat[,colnames(dat)%in%unique(group)]
pointdat$species <- as.vector(dat$species)
pointdat <- melt(pointdat, id="species")
pointdat$value <- cut(pointdat$value, breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0.01, 0.3, 0.4, 0.7), include.lowest=T)
colnames(pointdat) <- c("species", "groups", "indicator_value")
pointdat$species <- factor(pointdat$species, levels=factorlevel)
write.table(dat, outputdata, row.names=F, col.names=T, sep="\t", quote=F)
#if (filterabundance=="T"){
#	dat <- dat[dat$mean_abundance >= foldabundance, ]
#}else{
#	dat <- dat
#}

p1 <- ggplot(dat) +
geom_bar(aes(x=species,y=mean_abundance, fill=phylum),position="dodge", width=barsize,stat = "identity")+
theme_bw() +
theme(axis.text.y = element_text(size=xsize, face="bold"), legend.position = "left") +
guides(fill = guide_legend(ncol=1))+
xlab("") +
ylab("mean relative abundance (%)") +
coord_flip()

p2 <- ggplot(pointdat) +
geom_point(aes(x=species, y=groups, size=indicator_value, color=groups))+
theme_bw() +
xlab("") +
ylab("") +
scale_color_manual(values=mycolors) +
scale_size_manual(values=c(0.5,2,4,6), labels=c("0~0.25","0.25~0.5","0.5~0.75","0.75~1")) +
#theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank()) +
theme(axis.text.x=element_text(angle = 45, hjust = 1),axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
coord_flip()

svg(outfigure, height=height, width=width)
ggarrange(p1, p2, ncol=2, nrow=1, widths=c(8,5), byrow=F)
dev.off()

