rm(list=ls()) 

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputFilePrefix = args[2]

options(bitmapType='cairo')

#load Rcpp package first because of the error with reshape2 package
library(Rcpp)
library(grid)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

posmelt<-read.delim(inputFile, header=T,as.is=T, stringsAsFactor=F)

featureNumber<-length(unique(posmelt$Dinucleotide))
sampleNumber<-length(unique(posmelt$Sample))

height=max(featureNumber*100,2000)
width=max(sampleNumber*2000,3000)

axisTextSize=12
stripTextSize=12

ncols<-ifelse(featureNumber>=5, 2, 1)
maxPos<-max(posmelt$Position)

for (y in c("Count", "NormalizedCount")){
  if (! (y %in% colnames(posmelt))) {
    next
  }

  png( paste0(outputFilePrefix,"_",y, ".png"), width=width,height=height,res=300)
  m <- ggplot(posmelt, aes_string(x="Position",y=y,fill="Dinucleotide")) +
      geom_bar(stat="identity") +
      theme_bw()+
      theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(featureNumber)) + 
      xlim(-5, maxPos+5) +
      theme(legend.key.size = unit(0.4, "cm"), legend.position="right") +
      guides(fill= guide_legend(ncol=ncols,keywidth=1, keyheight=1.5))

  if (sampleNumber > 1){
    m <- m + facet_wrap(.~Sample)
  }

  print(m)
  dev.off()
}