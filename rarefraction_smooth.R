#!/usr/bin/Rscript
args<-commandArgs(T)

if (length(args)!=3){
	cat("#############################\n")
	cat("Function:the program is designed to plot the rarefaction smooth\n")
	cat("Usage:Rscript table of rarefaction output ylab(chao1,shannon) \n")
	cat("############################\n")
	quit()
}	

library(ggplot2)
myPalette <- c("#FF4040", "#228B22", "#FFFF33", "#0000FF", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#458B74", "#A52A2A", "#8470FF", "#53868B", "#8B4513", "#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32","#68838B","#7A8B8B","#CDBE70","#D3D3D3","#EEA2AD","#EE9572","#FFA500","#8B4789","#548B54","#2E8B57","#F4A460","#3A5FCD","#6959CD","#8B4726","#B0E0E6","#EECBAD","#CD69C9","#EE4000","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#8B3626","#00C5CD","#008B45","#D2B48C","#EED2EE","#CDCD00","#CD3278","#836FFF","#FF4500","#7D26CD","#7FFF00","#6495ED","#006400","#9A32CD","slateblue2","slategray3","tan1","violetred2","thistle2","#00868B","royalblue2","palegreen1","palevioletred2","powderblue","purple3","orangered2","olivedrab2","maroon2","orange3","lightskyblue1","lightgoldenrod2","firebrick2","deepskyblue","cyan1","bisque","chartreuse","khaki2","midnightblue","seashell2","plum2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","#CD853F","#8B7B8B","#8B3626","#00C5CD","#008B45","#D2B48C","#EED2EE","#CDCD00","#CD3278","#836FFF","#FF4500","#7D26CD","#7FFF00","#6495ED","#006400","#9A32CD","slateblue2","slategray3","tan1","violetred2","thistle2","#00868B","royalblue2","palegreen1","palevioletred2","powderblue","purple3","orangered2","olivedrab2","maroon2","orange3","lightskyblue1","lightgoldenrod2","firebrick2","deepskyblue","cyan1","bisque","chartreuse","khaki2","midnightblue","seashell2","plum2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","#F4A460","#3A5FCD","#6959CD","#8B4726","#B0E0E6","#EECBAD","#CD69C9","#EE4000","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#8B3626","#00C5CD","#008B45","#D2B48C","#EED2EE","#CDCD00","#CD3278","#836FFF","#FF4500","#7D26CD","#7FFF00","#6495ED","#006400","#9A32CD","slateblue2","slategray3","tan1","violetred2","thistle2","#00868B","royalblue2","palegreen1","palevioletred2","powderblue","purple3","orangered2","olivedrab2","maroon2","orange3","lightskyblue1","lightgoldenrod2","firebrick2","deepskyblue","cyan1","bisque","chartreuse","khaki2","midnightblue","seashell2","plum2","yellowgreen","tomato3","violetred4","magenta","chocolate2","darkorange3","#CD853F","#8B7B8B","#8B3626","#00C5CD","#008B45","#D2B48C","#EED2EE","#CDCD00","#CD3278","#836FFF","#6495ED", "#8B6508", "#556B2F", "#CD5B45", "#483D8B", "#EEC591", "#8B0A50", "#696969", "#8B6914", "#008B00", "#8B3A62", "#20B2AA", "#8B636C", "#473C8B", "#36648B", "#9ACD32","#68838B","#7A8B8B","#CDBE70","#D3D3D3","#EEA2AD","#EE9572","#FFA500","#8B4789","#548B54","#2E8B57","#F4A460","#3A5FCD","#6959CD","#8B4726","#B0E0E6","#EECBAD","#CD69C9","#EE4000","#436EEE","#8B8B00","#8B7E66","#CD853F","#8B7B8B","#8B3626","#00C5CD","#008B45","#D2B48C","#EED2EE","#CDCD00","#CD3278","#836FFF","#FF4500","#7D26CD","#7FFF00","#6495ED","#006400","#9A32CD","slateblue2","slategray3","tan1","violetred2")
table = args[1]
output = paste(args[2],'svg',sep=".")

data <- read.table(table, header=T, sep='\t', check.names=F)
tmpsample <- colnames(data)[4:ncol(data)]
sample <- vector()
species <- vector()
depth <- vector()
for (i in 4:ncol(data)){
	samplename <- tmpsample[i-3]
	na <- c("n/a")
	tmpy = data[[tmpsample[i-3]]]
	tmp <- as.vector(tmpy[! tmpy %in% na])
	species <- as.numeric(c(species, tmp))
	tm <- data[,2][1:length(tmp)]
	depth <- as.numeric(c(depth, tm))
	for (j in 1:length(tmp)){
		sample <- c(sample, samplename)
	}
}
dat <- cbind.data.frame(depth, species, sample)

head(dat)
ylabs <- args[3]
svg(output,width=10,height=7)
ggplot(dat, aes(depth, species, color = sample)) + 
geom_smooth(se=F, method = "lm",formula = y ~ log(x)) + 
scale_color_manual(values=myPalette) + 
theme_bw() + 
theme(panel.grid=element_blank(), plot.title = element_text(lineheight=2.5, face="bold",hjust=0.5)) + 
ylab(ylabs) +xlab("Number of Reads Sampled") + 
labs(title="Rarefaction Curves") + 
guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) + 
theme(legend.position = "right", legend.text=element_text(size=7))

dev.off()
