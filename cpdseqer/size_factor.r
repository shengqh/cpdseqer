count_file='uv_control.count'
output_file='uv_control.sizefactor'
outfilePrefix='uv_control'
setwd('/scratch/cqs/shengq2/guoyan/test/size_factor')

library("data.table")
library(edgeR)

bedCM <- fread(count_file,header=T,data.table=F) # bed Count Matrix
bedCM<-bedCM[,-1]
normFac <- calcNormFactors(bedCM,method='TMM')
librarysize = colSums(bedCM)
sf<-data.frame("Sample"=names(normFac), "SizeFactor"=normFac, "LibrarySize"=librarysize)
write.table(sf, file=output_file, sep="\t", row.names=F, col.names=T, quote=F)
