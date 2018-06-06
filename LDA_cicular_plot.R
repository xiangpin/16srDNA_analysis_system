#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to plot the cicular LDA.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("LEfSe_res", nargs=1, help="the LEfSe output.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o","--output", default="LEfSe_LDA_plot.svg", help="the output file, default is LEfSe_LDA_plot.svg.")
myparser$add_argument("-l","--level", default="N", help="the level tax, default is N.")
myparser$add_argument("-s","--labelsize", default=4, help="the tax label size, default is 4.")


args <- myparser$parse_args()

library(ggplot2)
mycolors <- c('#00AED7', '#FD9347','#C1E168','#319F8C',"#F8766E","#7CAD00","#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")
otu <- args$LEfSe_res
groupfile <- args$sample
outfigure <- args$output
level <- as.character(args$level)
sizelabel <- as.numeric(args$labelsize)
da <- read.table(otu, header=F, check.names=FALSE, sep='\t', comment.char="", quote="")
groupTable <- read.table(groupfile, header=T, check.names=F, comment.char="")
backdata <- data.frame(c("k","p","c","o","f","g","s"), c(1,2,3,4,5,6,7))
colnames(backdata) <- c("tax", "level")
group <- as.vector(groupTable[,2])
lenuniq <- length(unique(group))
dat <- data.frame()
for (i in 1:lenuniq){
	tmp <- da[da$V3==unique(group)[i],]
	dat <- rbind(dat, tmp)
}
tax <- as.vector(dat$V1)
#head(dat)
taxonomy <- vector()
if (level == "T"){
	taxonomy <- tax
}else{
	for (i in 1:length(tax)){
		tmp <- strsplit(tax[i], "[.]")[[1]]
		tmpgroup <- as.character(dat[dat$V1==tax[i],]$V3)
		tmplen <- length(tmp)
		#tmptaxT <- as.character(backdata[backdata$level==tmplen,1])
		taxonomy[i] <- paste(tmp[tmplen],tmpgroup, sep="_")
	}
}
#print(taxonomy)
dat$V1 <- NULL
data <- data.frame(cbind(taxonomy, dat), check.names=F)
colnames(data) <- c("taxonomy", "nomeans", "group", "LDA", "FDR")
data$taxonomy <- factor(data$taxonomy, levels=data$taxonomy)
#head(data)
sequence_length = length(unique(data$taxonomy))
first_sequence = c(1:(sequence_length%/%2))
second_sequence = c((sequence_length%/%2+1):sequence_length)
first_angles = c(90 - 180/length(first_sequence) * first_sequence)
second_angles = c(-90 - 180/length(second_sequence) * second_sequence)

#angle = 360/(2*pi)*rev( pi/2 + seq( pi/14, 2*pi-pi/sequence_length, len=sequence_length))

p <- ggplot(data, aes(x=taxonomy, y=LDA, fill=group))+
geom_bar(width=0.45, stat='identity')+
theme_light()+
theme(axis.text.y=element_text(angle=0))+
geom_hline(yintercept = 1, linetype='dashed', size=0.5, color="grey")+
geom_hline(yintercept = 2, linetype='dashed', size=0.5, color="grey")+
geom_hline(yintercept = 3, linetype='dashed', size=0.5, color="grey")+
ylab("LDA Score (log 10)")+
scale_fill_manual(values=mycolors)
p <- p +
theme(axis.text.x = element_text(angle=c(first_angles, second_angles), size=sizelabel, vjust = 1, hjust=1), axis.line.x=element_blank())
p <- p + coord_polar("x") #+ scale_y_continuous(breaks=cumsum(data$LDA))
#cumsum
svg(outfigure, width=10, height=10)
p
dev.off()
