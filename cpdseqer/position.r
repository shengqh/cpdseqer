rm(list=ls()) 
outFile='cpd.vis'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parSampleFile4=''
parFile1='out_chr1_chr2.txt'
parFile2=''
parFile3=''

setwd("C:/Users/sheng/projects/guoyan/cpd")

options(bitmapType='cairo')

maxFeature=50

groupFileList=parSampleFile1
visLayoutFileList=parSampleFile2
positionFile = parFile1
totalCountFile<-parFile2

#load Rcpp package first because of the error with reshape2 package
library(Rcpp)
library(grid)
library(reshape2)
library(ggplot2)
library(RColorBrewer)


position<-read.delim(positionFile, header=T,as.is=T)

posmelt<-melt(position, id='X')
colnames(posmelt)<-c("Position", "Dinucleotide", "Count")


featureNumber<-length(unique(posmelt$Dinucleotide))

height=max(featureNumber*100,2000)
width=3000

axisTextSize=12
stripTextSize=12

ncols<-ifelse(featureNumber>=5, 2, 1)
maxPos<-max(posmelt$Position)
png(paste0(outFile,".allPositionBar.png"),width=width,height=height,res=300)
m <- ggplot(posmelt, aes(x=Position,y=Count,fill=Dinucleotide)) +
    geom_bar(stat="identity") +
    theme_bw()+
    theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
    scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(featureNumber)) + 
    xlim(-5, maxPos+5) +
    theme(legend.key.size = unit(0.4, "cm"), legend.position="right") +
    guides(fill= guide_legend(ncol=ncols,keywidth=1, keyheight=1.5))

print(m)
dev.off()

