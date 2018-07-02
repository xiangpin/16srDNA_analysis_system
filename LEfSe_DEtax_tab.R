#!/usr/bin/Rscript
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to build the DE tax tab with the LEfSe_output and input.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("LEfSe_res", nargs=1, help="the LEfSe output.")
myparser$add_argument("LEfSe_input", nargs=1, help="the LEfSe input file.")
myparser$add_argument("-o","--output", default="tax_LEfSe_tab.xls", help="the output file, default is tax_LEfSe_tab.xls.")


args <- myparser$parse_args()
library(stringr)
otu <- args$LEfSe_res
input2 <- args$LEfSe_input
output <- args$output

lefseres <- read.table(otu, header=F, row.names=1, sep="\t")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "\\.", "|")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "S24_", "S24-")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "615J_", "615J-")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "4C0d_", "4C0d-")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "JG30_KF_CM45", "JG30-KF-CM45")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "TM7_3", "TM7-3")
rownames(lefseres) <- str_replace_all(rownames(lefseres), "Rs_045", "Rs-045")
lefseinput <- read.table(input2, header=T, row.names=1, sep="\t", check.names=F)
#groups <- unique(colnames(lefseinput))
#tmplefse <- lefseinput[-1, ]
names <- as.vector(as.matrix(lefseinput[1,]))
groups <- unique(as.vector(as.matrix(lefseinput[1,])))
lefseinput <- lefseinput[-1, ]
tmplefse <- data.frame(matrix(as.numeric(as.matrix(lefseinput[1:nrow(lefseinput), 1:ncol(lefseinput)])), ncol=ncol(lefseinput)))
rownames(tmplefse) <- rownames(lefseinput)
colnames(tmplefse) <- colnames(lefseinput) 
lefseinput <- tmplefse
names(lefseinput) <- names
print(groups)
da <- data.frame()
for (i in 1:length(groups)){
	tmp <- lefseres[lefseres$V3==groups[i],]
	da <- rbind(da, tmp)
}
da$V2 <- NULL
colnames(da) <- c("Sign_Group", "LDA" ,"FDR")
head(da, 3)
keeptax <- rownames(da)
dat <- lefseinput[rownames(lefseinput)%in%keeptax, , drop=F]
dattmp <- data.frame(tmplefse[rownames(tmplefse)%in%keeptax, , drop=F])
Totalmean <- apply(dattmp, 1, mean)
datt <- cbind(Totalmean, dattmp)
head(dat,4)
dtmean <- sapply(split.default(dat, names(dat)), rowMeans)
head(dtmean)
colnames(dtmean) <- sapply(colnames(dtmean), paste, "mean", sep="_")
newdata <- data.frame(merge(da, dtmean, by=0),check.names=F)
rownames(newdata) <- newdata$Row.names
newdata$Row.names <- NULL
#head(newdata, 3)
newdata <- data.frame(merge(newdata, datt, by=0),check.names=F)
#newdata <- cbind(da, dtmean, datt)
rownames(newdata) <- newdata$Row.names
newdata$Row.names <- NULL
tax <- rownames(newdata)
newdata <- cbind(tax, newdata)
write.table(newdata, output, row.names=F, col.names=T, quote=F, sep="\t")
