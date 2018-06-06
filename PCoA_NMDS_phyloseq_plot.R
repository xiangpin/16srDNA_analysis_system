#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the PCoA and NMDS with OTU table.\\n Desinger: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created (biom).")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o","--output", default="otu_tree_tax_phyloseq.class", help="the output file, default is otu_tree_tax_phyloseq.class.")
myparser$add_argument("-m2","--pcoaimage", default="PCoA_plot.svg", help="the PCoA plot out, default is PCoA_plot.svg")
myparser$add_argument("-m3", "--nmdsimage", default="NMDS_plot.svg", help="The NMDS plot out, default is NMDS_plot.svg")
myparser$add_argument('-s', '--size', default=2, help="the point size of the PCoA and NMDS,default is 2.")


args <- myparser$parse_args()
otu_tab <- args$otu_tab
sample <- args$sample
size <- as.numeric(args$size)
outputfile <- args$output
image <- args$richnessimage
pcoa <- args$pcoaimage
nmds <- args$nmdsimage

groupTable <- read.table(sample, header=T, check.names=F,comment.char="", row.names=1)
group <- groupTable$group
legendtmp<-unique(group)
if (length(legendtmp) > 2){
	mycolors <- c('#00AED7','#FD9347','#C1E168','#319F8C', "#FF4040", "#228B22", "gold", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")

}else{
	mycolors <- c('#00AED7','#FD9347')
}
shapes <- c(16, 17, 15, 3, 7, 8, 13, 11, 18, 9)
library(phyloseq)
library(ggplot2)
require(egg)

psotu <- import_biom(otu_tab)
pssample <- import_qiime_sample_data(sample)
ps = merge_phyloseq(psotu, pssample)
save(ps,file=outputfile)
PCoA = ordinate(ps, method="PCoA")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")
x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x,groupTable, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

gg <- merge(x,aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~group,x,mean),by="group")
#head(gg)
p1<-plot_ordination(ps, PCoA, color="group", shape="group") + 
xlab(xlab_text)+
ylab(ylab_text)+
theme(aspect.ratio=1.5)+
theme_bw()+
scale_colour_manual(values=mycolors)+
scale_shape_manual(values=shapes) +
geom_point(size = size) +
geom_vline(xintercept=0, linetype=2)+
geom_hline(yintercept = 0, linetype=2)+
geom_segment(aes(x=gg$mean.x, y=gg$mean.y, xend=gg$Axis.1, yend=gg$Axis.2, color=gg$group), size=0.3)+
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank())

p2 <- ggplot(data = x) +
geom_boxplot(aes(x=group, y=Axis.1,fill=group), outlier.size=0, outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
labs(title=title)+
coord_flip()+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.y = element_text(angle = 50, hjust = 1, face="bold" ),legend.position="none",plot.title = element_text(face="bold",lineheight=25,hjust=0.5))

p3 <- ggplot(data = x) + 
geom_boxplot(aes(x=group, y=Axis.2,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 40, hjust = 1,face="bold"),legend.position="none")

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

svg(pcoa)
ggarrange(p2,p1,empty,p3,ncol=2, nrow=2, heights=c(3,9),widths=c(9,3),byrow=F)
dev.off()

NMDS = ordinate(ps, method="NMDS")
svg(nmds,width=6,height=5)
plot_ordination(ps, NMDS, color="group", shape="group") +
labs(title="NMDS")+
theme(aspect.ratio=1.5)+
theme_bw()+
scale_colour_manual(values=mycolors)+
scale_shape_manual(values=shapes) +
geom_point(size=size)+
geom_vline(xintercept=0, linetype=2)+
geom_hline(yintercept = 0, linetype=2)+
theme(panel.grid = element_blank(),plot.title = element_text(face="bold",lineheight=25,hjust=0.5),panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

