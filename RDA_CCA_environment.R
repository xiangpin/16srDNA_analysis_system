suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the RDA and CCA.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("environment", nargs=1, help="the environment factor file.")
myparser$add_argument("-r", "--RDAoutput", default="rda_plot.svg", help="the output file, default is pca_plot.svg.")
myparser$add_argument("-c", "--CCAoutput", default="cca_plot.svg", help="the output file, default is cca_plot.svg.")

args <- myparser$parse_args()

library(vegan)
library(ggplot2)
library(factoextra)
library(gridExtra)
library(ggExtra)
library(egg)

otu <- args$otu_tab
groupfile <- args$sample
envfile <- args$environment
RDAfigure <- args$RDAoutput
CCAfigure <- args$CCAoutput

otuTable<-read.table(otu, header=T, skip=1, row.names=1, check.names=FALSE, sep='\t', comment.char="")
otuTable$taxonomy <- NULL
groupTable <- read.table(groupfile, header=T, as.is=T, check.names=F, comment.char="", row.names=1, sep='\t')
envFactor <- data.frame(t(read.table(envfile, header=T, row.names=1, check.names=FALSE, sep='\t', comment.char="")),check.names=F)

sample <- colnames(otuTable)

treatment <- vector()
for (i in 1:length(sample)){
	treatment[i] <- as.character(groupTable[sample[i],])
}

treatment <- factor(treatment)

legendtmp<-unique(treatment)
if (length(legendtmp)>2){
        mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C',"#FF4040", "#228B22", "#FF7F00", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
	#mycolors <- c('#00AED7','#C1E168','#FD9347','#319F8C')
}else{
        mycolors <- c("#FF7F00", "darkgreen")
}

##FFFF33
###################### transver the input data, such as otutable, envFactor.################
tmp <- data.frame(t(otuTable))
colnames(tmp) <- rownames(otuTable)
rownames(tmp) <- colnames(otuTable)
otuTable <- tmp

tmp <- data.frame(t(envFactor))
colnames(tmp) <- rownames(envFactor)
rownames(tmp) <- colnames(envFactor)
envFactor <- tmp
###########################RDA analysis########################
rda_res <- rda(otuTable, envFactor)
xlab_text <- paste("RDA1 (", round(rda_res$CCA$eig[1]/sum(rda_res$CCA$eig)*100,2), "%)", sep="")
ylab_text <- paste("RDA2 (", round(rda_res$CCA$eig[2]/sum(rda_res$CCA$eig)*100,2), "%)", sep="")
fit <- envfit(rda_res, envFactor, perm=10000)

for(i in 1:nrow(fit$vector$arrow)){
        if(fit$vector$pval[i]<=0.01) {
                rownames(fit$vectors$arrows)[i] <- paste(rownames(fit$vectors$arrows)[i], "**", sep="")
        }
        else if(fit$vector$pval[i]<=0.05) {
                rownames(fit$vectors$arrows)[i] <- paste(rownames(fit$vectors$arrows)[i], "*", sep="")
        }
}

cat("######################Ploting the RDA figure\n")
envfactor<-data.frame(fit$vectors$arrows)
scrs<-scores(rda_res,display=c("sp","wa","lc","bp","cn"))
env <- scrs$biplot
rownames(env) <- rownames(envfactor)
#multiplier <- vegan:::ordiArrowMul(env)
env <- as.data.frame(env*5)

#attributes(scrs)
#attributes(rda_res)
#x<-data.frame(rda_res$CCA$u[, 1:2])
x<-data.frame(scrs$sites)
x <- merge(x,groupTable, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
gg <- merge(x,aggregate(cbind(xcenter=RDA1,ycenter=RDA2)~group,x,mean),by="group")

p <- ggplot(data=x, aes(RDA1, RDA2, group)) + geom_point(aes(color=group,shape=group),size=3) +scale_color_manual(values=mycolors)+
theme_bw() +
geom_segment(aes(x=gg$xcenter, y=gg$ycenter, xend=gg$RDA1, yend=gg$RDA2,color=gg$group))+
xlab(xlab_text)+
ylab(ylab_text)+
geom_segment(data=env, aes(x=0,y=0, xend=RDA1, yend=RDA2), arrow=arrow(length=unit(0.3, "cm")),size=0.8)+
geom_text(data=as.data.frame(env*1.2),aes(RDA1, RDA2, label = rownames(env)),size=2.7)+
geom_vline(xintercept = 0,linetype='dashed',size=0.4)+
geom_hline(yintercept = 0,linetype='dashed',size=0.4)+
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank())

p2 <- ggplot(data = x) +
geom_boxplot(aes(x=group, y=RDA1,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
labs(title="RDA - RDA1 vs RDA2")+
coord_flip()+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.y = element_text(angle = 50, hjust = 1, face="bold" ),legend.position="none",plot.title = element_text(face="bold",lineheight=25,hjust=0.5))

p3 <- ggplot(data = x) + geom_boxplot(aes(x=group, y=RDA2,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle=40, hjust=1,face="bold"),legend.position="none")

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

svg(RDAfigure)
#png(RDAfigure, width = 880, height = 880)
ggarrange(p2,p,empty,p3,ncol=2, nrow=2, heights=c(3,9),widths=c(9,3),byrow=F)
dev.off()

cat("####################The RDA figure has been plotted.\n")


##########CCA analysis
cca_res <- cca(otuTable ~ ., data=envFactor)
x_lab <- paste("CCA1 (", round(cca_res$CCA$eig[1]/sum(cca_res$CCA$eig)*100,2), "%)", sep="")
y_lab <- paste("CCA2 (", round(cca_res$CCA$eig[2]/sum(cca_res$CCA$eig)*100,2), "%)", sep="")

scrs<-scores(cca_res,display=c("sp","wa","lc","bp","cn"))
fit <- envfit(cca_res, envFactor, perm=10000)
for(i in 1:nrow(fit$vector$arrow)){
        if(fit$vector$pval[i]<=0.01) {
                rownames(fit$vectors$arrows)[i] <- paste(rownames(fit$vectors$arrows)[i], "**", sep="")
        }
        else if(fit$vector$pval[i]<=0.05) {
                rownames(fit$vectors$arrows)[i] <- paste(rownames(fit$vectors$arrows)[i], "*", sep="")
        }
}

envfactor<-data.frame(fit$vectors$arrows)
env <- scrs$biplot
rownames(env) <- rownames(envfactor)
multiplier <- vegan:::ordiArrowMul(env)
#env <- as.data.frame(env*3)
env <- as.data.frame(env*3)

cat("######################Ploting the RDA figure\n")
x<-data.frame(scrs$sites)
x <- merge(x,groupTable, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL
gg <- merge(x,aggregate(cbind(xcenter=CCA1,ycenter=CCA2)~group,x,mean),by="group")

cca_p <- ggplot(data=x, aes(CCA1, CCA2, group)) + geom_point(aes(color=group,shape=group),size=3) +scale_color_manual(values=mycolors)+
theme_bw() +
geom_segment(aes(x=gg$xcenter, y=gg$ycenter, xend=gg$CCA1, yend=gg$CCA2,color=gg$group))+
xlab(x_lab)+
ylab(y_lab)+
geom_segment(data=env, aes(x=0,y=0, xend=CCA1, yend=CCA2), arrow=arrow(length=unit(0.3, "cm")),size=0.7)+
geom_text(data=as.data.frame(env*1.1),aes(CCA1, CCA2, label = rownames(env)),size=3)+
geom_vline(xintercept = 0,linetype='dashed',size=0.4)+
geom_hline(yintercept = 0,linetype='dashed',size=0.4)+
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank())

cca_p2 <- ggplot(data = x) +
geom_boxplot(aes(x=group, y=CCA1,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
labs(title="CCA - CCA1 vs CCA2")+
coord_flip()+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.y = element_text(angle = 50, hjust = 1, face="bold" ),legend.position="none",plot.title = element_text(face="bold",lineheight=25,hjust=0.5))

cca_p3 <- ggplot(data = x) + geom_boxplot(aes(x=group, y=CCA2,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle=40, hjust=1,face="bold"),legend.position="none")


svg(CCAfigure)
#png(CCAfigure,width = 880, height = 880)
ggarrange(cca_p2, cca_p,empty,cca_p3,ncol=2, nrow=2, heights=c(3,9),widths=c(9,3),byrow=F)
dev.off()

cat("####################The RDA figure has been plotted.\n")
