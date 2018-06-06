#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(clusterSim))
docstring = "Description: The script is designed to run the enterotype stratification with the PAM(Partitioning around medoids).\\nDeveloper: Shuangbin Xu\\nEmail: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("tab", nargs=1, help="the Genus ralative abundance(sums is 1).")
myparser$add_argument("-p","--pseudocount", default=0.000001, help="the pseudo count of abundance distributions to avoid zero in the numerator and/or denominator of equation, default is 0.000001.")
myparser$add_argument("-n","--clusternum", default=3, help="the number of cluster center, default is 3.")
myparser$add_argument("-o","--output", default="cluster_sample_enterotypes.txt", help="the output file, default is cluster_sample_enterotypes.txt.")
myparser$add_argument("-k","--indexplot", default="enterotypes_cluster_nums_choose.svg", help="the to choose the number of cluster, default is enterotypes_cluster_nums_choose.svg.")
myparser$add_argument("-c","--pcaplot", default="enterotypes_pca_plot.svg", help="the visualization of the cluster, default is enterotypes_pca_plot.svg.")

args <- myparser$parse_args()

suppressPackageStartupMessages(library(adehabitat))
suppressPackageStartupMessages(library(egg))
input <- args$tab
pseudo <- as.numeric(args$pseudocount)
clusternum <- as.numeric(args$clusternum)
clustersample <- args$output
index_plot <- args$indexplot
pcafigure <- args$pcaplot
mycolors <- c('#00AED7','#FD9347','#C1E168','#319F8C', "#FF4040", "#228B22", "gold", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
shapes <- c(16, 17, 15, 3, 7, 8, 13, 11, 18, 9)

data <- read.table(input, sep="\t", header=T, check.names=F, row.names=1)

dist.JSD <- function(inMatrix, pseudocount=pseudo, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

	inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) {
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
				as.vector(inMatrix[,j]))
			}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix)
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
	require(cluster)
	cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
	return(cluster)
}
#noise.removal <- function(dataframe, percent=0.001, top=NULL){
#	dataframe->Matrix
#	percent=100*(10^(mean(log10(rowSums(Matrix)/sum(rowSums(Matrix)))) + 1.96*sd(log10(rowSums(Matrix)/sum(rowSums(Matrix))))))
#	bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
#	Matrix_1 <- Matrix[bigones,]
#	return(Matrix_1)
#}
noise.removal<-function(dataframe, percent=NULL, top=NULL, bysample=TRUE){
noise.removal.global <- function(dataframe, percent=NULL, top=NULL){
dataframe->Matrix
if (is.null(top)){
if (is.null(percent)) {
percent=100*(10^(mean(log10(rowSums(Matrix)/sum(rowSums(Matrix)))) + 1.96*sd(log10(rowSums(Matrix)/sum(rowSums(Matrix)))))) # this is percent
}
bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent  #this is percent ### noize filter
Matrix_1 <- Matrix[bigones,]
print(percent)
}
else {
bigones <- rev(order(apply(Matrix, 1, mean)))[1:top]
Matrix_1 <- Matrix[bigones,]
}
return(Matrix_1)
}
noise.removal.spec <- function(dataframe, percent=NULL, top=NULL){
dataframe->Matrix
bigones<-NULL
for(i in 1:dim(Matrix)[2]){
        if (is.null(top)){
                if (is.null(percent)){
                        percent<- 100*(10^(median(log10(Matrix[Matrix[,i]>0,i])) + 1.96*sd(log10(Matrix[Matrix[,i]>0,i]))))
                        cat(percent,"\n")
                }
                        tmp <- which(Matrix[,i] > percent/100)
                        bigones<-unique(c(bigones,tmp))
        }
        else {
                tmp <- rev(order(Matrix[,i]))[1:top]
                bigones<-unique(c(bigones,tmp))
        }
}
return(Matrix[bigones,])
print(percent)
}
if(bysample){noise.removal.spec(dataframe, percent=percent, top=top)}
else {noise.removal.global(dataframe, percent=percent, top=top)}
}

#data = noise.removal(data, 0.00001)
data = data
pca <- prcomp(t(data))
#head(pca)
x<-data.frame(pca$x[, 1:2])
data.dist=dist.JSD(data)

data.cluster=pam.clustering(data.dist, k=clusternum)
group <- vector()
for (i in 1:length(data.cluster)){
        group[i] <- paste("cluster",data.cluster[i],sep="")
}
sample <- colnames(data)
clgroup <- data.frame(cbind(sample, group))
write.table(clgroup, clustersample, sep="\t", row.names=F, col.names=T, quote=F)
rownames(clgroup) <- clgroup$sample
clgroup$sample <- NULL

x <- merge(x, clgroup, by=0)
rownames(x) <- x$Row.names
x$Row.names <- NULL
#print(clgroup)
#require(clusterSim)
nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL

for (k in 1:20) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
print(obs.silhouette)

svg(index_plot, width=5, height=5)
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")
dev.off()

#obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
#obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1)
#dev.new()
#svg("test2.svg", width=5, height=5)
#s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")
#dev.off()

#plot 2
#obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
#print(attributes(obs.pcoa))
#dev.new()
#svg("test3.svg", width=5, height=5)
#s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")
#dev.off()

ev <- pca$sdev^2
vp <- ev*100/sum(ev)
xlab_text <- paste("PC1 (", round(vp[1],2), "%)")
ylab_text <- paste("PC2 (", round(vp[2],2), "%)")

gg <- merge(x,aggregate(cbind(mean.x=PC1,mean.y=PC2)~group,x,mean),by="group")
p <- ggplot(data=x, aes(PC1, PC2, group)) +
geom_point(aes(color=group,shape=group),size=2) +
stat_ellipse(aes(color=group), linetype="dashed", level=0.95)+
scale_color_manual(values=mycolors) +
scale_shape_manual(values=shapes) +
theme_bw() +
geom_segment(aes(x=gg$mean.x, y=gg$mean.y, xend=gg$PC1, yend=gg$PC2,color=gg$group),size=0.3)+
xlab(xlab_text)+
ylab(ylab_text)+
geom_vline(xintercept = 0,linetype='dashed',size=0.6)+
geom_hline(yintercept = 0,linetype='dashed',size=0.6)+
#guides(fill = guide_legend(ncol=1))+
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12),legend.title=element_blank())


p2 <- ggplot(data = x) +
geom_boxplot(aes(x=group, y=PC1,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
labs(title="PCA - PC1 vs PC2")+
coord_flip()+
theme_bw() +
theme(panel.grid = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text.y = element_text(angle = 50, hjust = 1, face="bold" ),
legend.position="none",
plot.title = element_text(face="bold",lineheight=25,hjust=0.5))

p3 <- ggplot(data = x) + 
geom_boxplot(aes(x=group, y=PC2,fill=group),outlier.size=0,outlier.shape=NA)+
scale_fill_manual(values=mycolors)+
xlab("")+
ylab("")+
theme_bw() +
theme(panel.grid = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text.x = element_text(angle=40,hjust = 1,face="bold"),
legend.position="none")

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


svg(pcafigure)
ggarrange(p2,p,empty,p3,ncol=2, nrow=2, heights=c(3,9),widths=c(9,3),byrow=F)
dev.off()

