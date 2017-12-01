suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the PCA.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o","--output", default="pca_plot.svg", help="the output file, default is pca_plot.svg.")

args <- myparser$parse_args()

library(vegan)
library(ggplot2)
library(factoextra)
library(gridExtra)
library(ggExtra)
library(egg)
otu <- args$otu_tab
groupfile <- args$sample
outfigure <- args$output
da <- read.table(otu, header=T, skip=1, row.names=1, check.names=FALSE, sep='\t', comment.char="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
da$taxonomy <- NULL
sample <- colnames(da)
group <- vector()
#da$taxonomy <- NULL
for(i in 1:length(sample)){
	if (sample[i] %in% groupTable$sample){
		group[i]<-as.character(groupTable[groupTable$sample==sample[i],]$group)
	} else {
		cat("sample", sample[i], " not found in sample info.\n")
		remove <- sample[i]
		sample <- setdiff(sample, remove)
		da[[remove]] <- NULL
	}
}
legendtmp<-unique(group)
if (length(legendtmp)>2){
	mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C',"#F8766E","#7CAD00","#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
	#mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C')
	#mycolors <- c("#F8766E","#7CAD00", "#00BEC3", "#C67CFF")

}else{

	mycolors <- c('#00AED7','#C1E168')
}

pcadat <- t(da)
pcadat2 <- as.data.frame(cbind(pcadat, group))
pca <- prcomp(pcadat)

x<-data.frame(pca$x[, 1:2])
x<-data.frame(cbind(x, group))
ev <- pca$sdev^2
vp <- ev*100/sum(ev)
xlab_text <- paste("PC1 (", round(vp[1],2), "%)")
ylab_text <- paste("PC2 (", round(vp[2],2), "%)")

gg <- merge(x,aggregate(cbind(mean.x=PC1,mean.y=PC2)~group,x,mean),by="group")
p <- ggplot(data=x, aes(PC1, PC2, group)) + geom_point(aes(color=group,shape=group),size=3) +scale_color_manual(values=mycolors)+
theme_bw() +
geom_segment(aes(x=gg$mean.x, y=gg$mean.y, xend=gg$PC1, yend=gg$PC2,color=gg$group))+
xlab(xlab_text)+
ylab(ylab_text)+
geom_vline(xintercept = 0,linetype='dashed',size=1)+
geom_hline(yintercept = 0,linetype='dashed',size=1)+
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank())


p2 <- ggplot(data = x) + 
geom_boxplot(aes(x=group, y=PC1,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
labs(title="PCA - P1 vs P2")+
coord_flip()+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.y = element_text(angle = 50, hjust = 1, face="bold" ),legend.position="none",plot.title = element_text(face="bold",lineheight=25,hjust=0.5))

p3 <- ggplot(data = x) + geom_boxplot(aes(x=group, y=PC2,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle=40,hjust = 1,face="bold"),legend.position="none")
empty <- ggplot() + 
              geom_point(aes(1,1), colour="white") +
              theme(                              
                plot.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(), 
                panel.background = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank()
              )
 
svg(outfigure)
#png(outfigure)
ggarrange(p2,p,empty,p3,ncol=2, nrow=2, heights=c(3,9),widths=c(9,3),byrow=F)
dev.off()
