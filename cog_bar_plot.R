#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
        cat("##################################################################\n")
        cat("Function: plot bar plot and p-value of data matrix by group file\n")
        cat("Usage: Rscript group_bar_plot_V2.R  dataTable sampleFile\n")
        cat("Make sure your dataTable format:\n
        ID A1 A2 A3 B1 B2 B3 ...\n
        ID1 1.0 2.0 1.5 2.1 1.3 1.3 ...\n
        ID2 1.2 2.1 1.4 2.0 1.4 1.4 ...\n
        ......\n")
        cat("##################################################################\n")
        quit()
}


library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(RGraphics)
library(scales)
library(vegan)

mycolors <- c('#00AED7', '#FD9347', '#C1E168', '#319F8C', "#FF4040", "#228B22", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32")


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

#Input dataTable sampleFile 
dataTable = args[1]
sampleFile = args[2]
outputFig = args[3]

a <- read.table(dataTable, header=T, sep='\t', row.names = 1, check.names=F, quote="",comment.char = "")
a <- data.frame(t(a), check.names=F)
set.seed(1000)
a <- rrarefy(a, min(rowSums(a)))
a <- data.frame(t(a), check.names=F)


sampleTable<-read.table(sampleFile, header=T, sep='\t', check.names=F)
s<-table(sampleTable$group)
#a <- read.table(dataTable, header=T, sep='\t', row.names = 1, check.names=F, quote="")
if (sum(s<2) != 0) {
        cat("the sample number for some group is smaller than 2!\n")
        cat("increase sample number or merge samples\n")
        quit()
}
#head(sampleTable)
sampleList <- names(a)
print(sampleList)
print(sampleTable[sampleTable$sample=="ZJ1.2-3",])
groupList <- vector()
for (i in 1:length(sampleList)) {
	#print (sampleTable[sampleTable$sample==sampleList[i],])
	groupList[i] <- as.character(sampleTable[sampleTable$sample==sampleList[i],]$group)
}

names(a) <- groupList
#test
ttest<-function(x){oneway.test(as.numeric(x) ~ names(a))$p.value}
a$p.value<-apply(a, 1, ttest)

#select pvalue less than 0.05
a<-a[a$p.value<=0.05&!is.na(a$p.value),]

if (nrow(a)==0){
        text=paste('There no significant different classification!')
        svg(outputFig,width=13)
        print(ggplot() + annotate("text", x = 1, y = 1, size=8, label = text) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))
        dev.off()
}else{
#sort according to pvalue
a<-a[order(a$p.value),]
a$attri<-rownames(a)
#select most significant pvalue and recalculate it
if (length(a$p.value)<20){
b<-a}else{
b<-a[1:20,]}
b$p.value<--1*log(b$p.value)/log(10)

#reshape data structure
c<-b[names(b)!="p.value"]
names(c) <- sampleList
names(c)[length(names(c))] <- "attri"

################################################################################
d<-melt(c, id="attri")
names(d) <- c("attri", "group", "value")

d_group <- as.character(d$group)
for (i in 1:length(d_group)) {
        d_group[i] <- as.character(sampleTable[sampleTable$sample==d_group[i],]$group)
}
d$group <- d_group
###############################################################################

t <- as.character(b$attri)
#print (t)
#b$ko <- substr(t,1,6)
tmpko <- strsplit(t, "[@]")
for (i in 1:length(tmpko)){b$ID[i] <- tmpko[[i]][1]}
#b$ko <- unlist(strsplit(t,"[@]"))[1]
#sort the attri
b$ID<-factor(b$ID, level=b[order(b$p.value),]$ID)
d$ID<-factor(b$ID, level=b[order(b$p.value),]$ID)

#calculate the mean sd se and ci
d_error <- summarySE(d, measurevar="value", groupvars=c("ID","group"))

#ko <- substr(unique(t),1,6)
ID <- vector()
for (i in 1:length(tmpko)){ID[i] <- tmpko[[i]][1]}
#ko <- unlist(strsplit(t,"[@]"))[1]
#anno <- substr(unique(t),8,100)
description <- vector()
for (i in 1:length(tmpko)){description[i] <- unlist(strsplit(tmpko[[i]][2], "[;]"))[1]}

df <- data.frame(ID,description)

df$ID<-factor(b$ID, level=b[order(b$ID),]$ID)

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
p1 <- tableGrob(df, theme = tt, rows=NULL)

p2 <- ggplot(d_error, aes(x=ID, y=value, fill=group)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.2, position=position_dodge(.9)) + coord_flip() + theme_bw() + xlab("") + theme(legend.position="left") + scale_fill_manual(values=mycolors) + scale_y_continuous(expand=c(0.01, 0))

p3 <- ggplot(data=b)+geom_point(aes(p.value, ID, color=p.value))+theme_bw()+theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())+ylab("")+xlab("-1 * log10(p-value)")+scale_colour_gradientn(colours=rainbow(20))

svg(outputFig, width=10, height=13); grid.arrange(p2, p3, p1, ncol=2, nrow=2, layout_matrix = rbind(c(1,2), c(3,3)), widths=c(2, 1), heights=c(4, 3.5))

dev.off()
}
