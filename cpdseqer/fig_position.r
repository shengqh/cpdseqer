args = commandArgs(trailingOnly=TRUE)

inputFile=args[1]
outfilePrefix=args[2]

options(bitmapType='cairo')

#load Rcpp package first because of the error with reshape2 package
required.pkg <- c("Rcpp", "grid", "reshape2", "ggplot2", "RColorBrewer")
for (pkg in required.pkg) {
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    require(pkg, character.only = TRUE)
  }
}


posmelt<-read.delim(inputFile, header=T,as.is=T, stringsAsFactor=F)

DINU4 = c("TT", "TC", "CT", "CC")
posmelt<-posmelt[posmelt$Dinucleotide %in% DINU4,]

dinus<-unique(posmelt$Dinucleotide)
dinus<-dinus[order(dinus)]

featureNumber<-length(dinus)
sampleNumber<-length(unique(posmelt$Sample))

width=max(sampleNumber*8,10)

axisTextSize=12
stripTextSize=12

ncols<-ifelse(featureNumber>=5, 2, 1)
maxPos<-max(posmelt$Position)

colors<-c("red", "blue", "yellow", "brown")
names(colors)<-dinus

for (ctype in c("Site", "Read")){
  for (y in c("Count", "NormalizedCount")){
    columnName=paste0(ctype, y)
    if (! (columnName %in% colnames(posmelt))) {
      next
    }
    cat(columnName, "\n")

    pdf( paste0(outfilePrefix,"_",columnName, ".pdf"), width=width,height=10, onefile=T)
    m <- ggplot(posmelt, aes_string(x="Position",y=columnName,fill="Dinucleotide")) +
        geom_bar(stat="identity") +
        theme_bw()+
        theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
        scale_fill_manual(values=colors) + 
        xlim(-5, maxPos+5) +
        theme(legend.key.size = unit(0.4, "cm"), legend.position="right") +
        guides(fill= guide_legend(ncol=ncols,keywidth=1, keyheight=1.5))

    if (sampleNumber > 1){
      m <- m + facet_wrap(.~Sample)
    }
    print(m)

    for (dinu in dinus) {
      dinuposmelt<-posmelt[posmelt$Dinucleotide == dinu,]
      m <- ggplot(dinuposmelt, aes_string(x="Position",y=columnName)) +
          geom_bar(stat="identity", color=colors[dinu]) +
          theme_bw()+
          theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
          xlim(-5, maxPos+5) +
          theme(legend.key.size = unit(0.4, "cm"), legend.position="right")+
          labs(title=paste0("Dinucleotide ", dinu))

      if (sampleNumber > 1){
        m <- m + facet_wrap(.~Sample)
      }
      print(m)
    }

    dev.off()
  }
}