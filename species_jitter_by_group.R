#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to calculate the abundant of species of group.\\n Designer: Shuangbin Xu\\n Email: xusbin@anjiemed.com."
parser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
	parser$add_argument("data", nargs=1, help="the species file.")
	parser$add_argument("sample", nargs=1, help="the sample file.")
	#parser$add_argument("-f","--outdata", default="outdata.txt", help="the output file.")
	parser$add_argument("-o", "--outimage", default="outimage.svg", help="the output image.")
	parser$add_argument("-t", "--tax", default="Phylum", help="the level of classify.")
	parser$add_argument("-H", "--height", default=5, help="the height of image, default is 5.")
	parser$add_argument("-w", "--width", default=6, help="the width of image, default is 6.")
	parser$add_argument("-z", "--size", default=9, help="the size of legend words, default is 9.")
	parser$add_argument("-n", "--topnum", default=20, help="the top number of classify, default is 20.")
	args <- parser$parse_args()
#myPalette <- c("#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
myPalette <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")

library("ggplot2")
library("reshape2")
library('scales')
data <- args$data
sampledat <- args$sample
#outdat <- args$outdata
outimage <- args$outimage
height <- as.numeric(args$height)
width <- as.numeric(args$width)
size <- as.numeric(args$size)
tax <- args$tax
topnum <- as.numeric(args$topnum)

specie <- read.table(data, header=T, sep="\t", row.names=1, check.names=F,comment.char="")
sampledat <- read.table(sampledat, header=T, sep="\t", check.names=F)
tmpspecie <- row.names(specie)
specie$sum <- apply(specie, 1, sum)
if ((length(specie$sum)) < 29){
        specie<- data.frame(specie[rev(order(specie$sum)),][1:length(specie$sum),],check.names=F)
} else
{
        specie <- data.frame(specie[rev(order(specie$sum)),][1:29,], check.names=F)
}
specie$sum <- NULL
specie$species <- rownames(specie)
factorLevel <- rownames(specie)
b<-melt(specie, id="species")
#b$species <- factor(b$species, levels=factorLevel)
sample <- as.vector(b$variable)
groups <- vector()
for (i in 1:length(sample)){
	tmp = as.character(sample[i])
	groups[i] <- as.character(sampledat[sampledat$sample==tmp,2])
}
b <- data.frame(cbind(b, groups))
#head(b)
b$species <- factor(b$species, levels=factorLevel)
svg(outimage, height=height, width=width)
ggplot(data=b) + 
#geom_point()+
geom_jitter(aes(groups, value, colour=species), width = 0.2)+
geom_smooth(aes(as.integer(groups), value, colour=species),method = 'lm', formula= y ~ x, alpha=0.3, show.legend=F) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face="bold",lineheight=25,hjust=0.5)) +
xlab("Group") +
ylab("Relative Abundance (%)") +
scale_color_manual(values=myPalette)+
guides(colour = guide_legend(keywidth = 0.85, keyheight = 0.85)) +
theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = size),legend.title=element_blank()) + 
labs(title=tax,fill="") 
dev.off()

