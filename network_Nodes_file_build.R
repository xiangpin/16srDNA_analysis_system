#!/usr/bin/Rscript
# -*- coding: utf-8 -*-
suppressPackageStartupMessages(library(argparse))
docstring = "Description:\\n The script is designed to build the nodes files of networks with the edges file of SparCC output.\\n Developer: Shuangbin Xu\\n Email: xusbin@anjiemed.com"
myparser <- ArgumentParser(description=docstring, formatter_class="argparse.RawTextHelpFormatter")
myparser$add_argument("otu_tab", nargs=1, help="the otu_tab file of the qiime created.")
myparser$add_argument("edgestab", nargs=1, help="the edges file of SparCC output.")
myparser$add_argument("-l","--level", default=6, help="the tax levels, default is 6 (as Genus)")
myparser$add_argument("-o","--output", default="Feature_node_files.txt", help="the output nodes file, default is Feature_node_files.txt.")

library(plyr)
args <- myparser$parse_args()

otutab <- args$otu_tab
edgesfile <- args$edgestab
level <- as.numeric(args$level)
output <- args$output

da <- read.table(otutab, skip=1, sep="\t", header=T, row.names=1, check.names=F, comment.char="")
dt <- read.table(edgesfile, header=T, sep="\t")
level = 6
tax <- as.vector(da$taxonomy)
species <- vector()
phylum <- vector()
for (i in 1:length(tax)){
        if (level==8){
                tmp <- tax[i]
                species[i] <- tmp
        }
        else {
                tmp <- strsplit(tax[i], "; ")[[1]][level]
                phylum[i] <- strsplit(tax[i], "; ")[[1]][2]
                species[i] <- tmp
        }
}
taxdat <- data.frame(cbind(species, phylum))
uniqphylum <- vector()
for (i in 1:length(unique(taxdat$species))){
        tmp <- as.character(unique(taxdat$species)[i])
        uniqphylum[i] <- as.character(taxdat[taxdat[,1]==tmp,2][1])
}
newtaxdat <- data.frame(unique(taxdat$species),uniqphylum)
colnames(newtaxdat) <- c("Feature", "type")
rownames(newtaxdat) <- newtaxdat$Feature
da$taxonomy <-NULL
da <- data.frame(cbind(species, da), check.names=F)
da <- data.frame(ddply(da, "species", numcolwise(sum)), check.names=F)
rownames(da) <- da$species
da$species <- NULL
da <- data.frame(prop.table(as.matrix(da), 2), check.names=F)
da <- da*100
da$mean <- apply(da, 1, mean)
nodelist <- unique(c(unique(as.vector(dt[,1])), unique(as.vector(dt[,1]))))
dt <- da[rownames(da)%in%nodelist, "mean", drop=F]
dt <- merge(newtaxdat, dt, by=0)
dt$Row.names <- NULL

write.table(dt, output, sep="\t", row.names=F, col.names=T, quote=F)
