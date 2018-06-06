#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the boxplot of alpha.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
#myparser$add_argument("-s","--size", default=10000, help="the size of the rare, default is 10000.")
#myparser$add_argument("-m","--method", default="Observed", help="the rare eveness, default is Observed, Choose (Chao1, se.chao1, Shannon, InvSimpson).")
myparser$add_argument("-o","--output", default="alpha_box_plot.svg", help="the output file, default is alpha_box_plot.svg.")
myparser$add_argument("-d","--alphatab", default="alpha_tab.txt", help="the output file, default is alpha_tab.txt.")

args <- myparser$parse_args()

library(vegan)
library(ggplot2)
library(reshape2)
library(agricolae)
otu <- args$otu_tab
groupfile <- args$sample
outfigure <- args$output
alphatab <- args$alphatab
da <- read.table(otu, header=T, skip=1, row.names=1, check.names=FALSE, sep='\t', comment.char="")
groupTable <- read.table(groupfile, header=T, check.names=FALSE, sep='\t')
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
sample <- rownames(alpha)
group <- vector()
for(i in 1:length(sample)){
	if (sample[i] %in% groupTable$sample){
		group[i]<-as.character(groupTable[groupTable$sample==sample[i],2])
	} else {
		cat("sample", sample[i], " not found in sample info.\n")
		remove <- sample[i]
		sample <- setdiff(sample, remove)
		alpha[remove,] <- NULL
	}
}

alpha <- cbind(alpha, group)
legendtmp<-unique(group)
if (length(legendtmp)>2){
	mycolors <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
	#mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C')

}else{

	mycolors <- c('#00AED7','#FD9347')
}
alpha <- melt(alpha, id.var = "group")

#svg(outfigure, width=10, height=4)
svg(outfigure, width=8, height=3.5)
#alpha$group<-factor(alpha$group,levels=c("T1","T5","T12"))
#alpha$group<-factor(alpha$group,levels=c("IFNa_T1","IFNa_T5","IFNa_T12","Mu_T1","Mu_T5","Mu_T12","WT_T1","WT_T5","WT_T12"))
ggplot(alpha) +
#geom_boxplot(aes(x=group, y=value, colour=group))+ 
geom_boxplot(aes(x=group, y=value, color=group))+#, fill=group))+
facet_wrap(~variable, scale="free", nrow=1) +
scale_color_manual(values= mycolors)+
theme_bw() +
ylab("Alpha Diversity Measure") +
theme(axis.title.x=element_blank(), axis.text.x = element_blank()) 
dev.off()

#png(paste(outfigure,"png",sep="."), width=14, height=6)
#ggplot(alpha) +
#geom_boxplot(aes(x=group, y=value, colour=group))+ 
#geom_boxplot(aes(x=group, y=value, color=group,fill=group))+
#facet_wrap(~variable, scale="free", nrow=1) +
#scale_color_manual(values= mycolors)+
#theme_bw() +
#ylab("Alpha Diversity Measure") +
#theme(axis.title.x=element_blank(), axis.text.x = element_blank())
#dev.off()
