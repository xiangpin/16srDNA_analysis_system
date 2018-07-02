#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to do the ANOVA-test and build a table.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("sample", nargs=1, help="the sample file of qiime input.")
myparser$add_argument("-o", "--output", default="significant_feature_table.xls", help="the output file, default is significant_feature_table.xls.")

args <- myparser$parse_args()

library(reshape2)
library(plyr)
#library(agricolae)
library(openxlsx)
otu_tab <- args$otu_tab
sample <- args$sample
outfile <- args$output

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

round_df <- function(df, digits) {
	nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
	df[,nums] <- round(df[,nums], digits = digits)
	df
}

da <- read.table(otu_tab, sep="\t", row.names=1, header=T, check.names=F, quote="", comment.char="")
da <- data.frame(t(da),check.names=F)
sampleinfo <- read.table(sample, sep="\t", row.names=1, header=T, check.names=F ,comment.char="")
#da <- data.frame(t(da),check.names=F)
dat <- merge(da, sampleinfo, by=0)
rownames(dat) <- dat$Row.names
dat$Row.names <- NULL
colnames(dat)[length(colnames(dat))] <- "group"
data <- melt(dat, id="group")
data_sd <- summarySE(data, measurevar="value", groupvars=c("variable", "group"))
tmpgroup <- as.vector(dat$group)

dat <- data.frame(t(dat),check.names=F)
names(dat) <- tmpgroup
dat <- dat[rownames(dat)!=c("group"),]
datt <- dat
ttest<-function(x){oneway.test(as.numeric(x) ~ names(dat))$p.value}
fvalue <- function(x){oneway.test(as.numeric(x) ~ names(dat))$statistic}
datt$p.value <- apply(dat, 1, ttest)
datt$F.value <- apply(dat, 1, fvalue)
datt$FDR <- p.adjust(datt$p.value, method = "BH")
tmpdata <- datt[, c("p.value", "F.value", "FDR")]
tmpdata <- round_df(tmpdata, 4)
datt$p.value <- NULL
datt$F.value <- NULL
datt$FDR <- NULL
data_sd <- round_df(data_sd,2)
data_sd <- within(data_sd, value<-paste(value, sd, sep="±"))
data <- dcast(data_sd, variable~group, value.var = "value")
rownames(data) <- data$variable
data$variable <- NULL
newdata <- merge(data, tmpdata, by=0)
newdata <- newdata[order(newdata$FDR),]
colnames(newdata)[1] <- "Features"
#write.xlsx(newdata, outfile)
write.table(newdata, outfile, sep="\t", col.names=T, row.names=F, quote=F)
