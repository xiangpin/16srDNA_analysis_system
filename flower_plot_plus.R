#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the Venn image.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of pca input.")
myparser$add_argument("-o","--output", default="flower_plot.svg", help="the output file, default is flower_plot.svg.")
args <- myparser$parse_args()
outfigure <- args$output
library(plotrix)
#mycolor <- c("#009ACD", "#66CD00", "#EEAD0E", "darkorchid1","#8B0000")
mycolors <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C',"#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
#ramp <- colorRamp(c("red","blue"))
#colorsflower <- c(rgb(135, 206, 235, 150, max = 255),rgb(25, 26, 25, 150, max = 255), rgb(15, 6, 35, 150, max = 255), rgb(135, 206, 235, 150, max = 255))
colorsflower <- c('#00AED796', '#FD934796', '#C1E16896', '#319F8C96',"#FF404096", "#228B2296", "#FFFF3396", "#0000FF96", "#984EA396")
flower_plot <- function(sample, value, start, a, b,
                              ellipse_col = colorsflower, 
                              circle_col = rgb(0, 162, 214, max = 255),
                              circle_text_cex = 1, labels=labels) {
par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(0.1,0.1,0.1,0.1))
plot(c(0,10),c(0,10),type="n")
n   <- length(sample)
deg <- 360 / n
res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t]
        )

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = circle_text_cex
            )

    } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = circle_text_cex
            )
    }
})
plotrix::draw.circle(x = 5, y = 5, r = 0.81, col = circle_col, border = circle_col)
text(x = 5, y = 5, labels=labels)
}

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
a <- c[,2:ncol(c)]
a <- data.frame(matrix(as.numeric(as.matrix(a[1:nrow(a),1:ncol(a)])), ncol=ncol(a)))
colnames(a) <- group
rownames(a) <- tmpspecie
a <- as.data.frame(t(apply(a, 1, function(x) tapply(x, colnames(a), mean))))
#head(a)
groups <- colnames(a)
otunum <- vector()
for (i in 1:length(groups)){
	tmp <- as.character(groups[i])
	otunum[i] <- sum(a[[tmp]] > 0)
}
data <- data.frame(cbind(groups,otunum))
#head(data)

svg(outfigure)#,width=5,height=5)
flower_plot(data$groups, data$otunum, 90, 0.5, 2, labels="OTU")
dev.off()
