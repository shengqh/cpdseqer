inputFile='genome.config.txt'
chromInfoFile='chromInfo.txt'
block='100000'
outfilePrefix='genome'
setwd('/gpfs23/scratch/cqs/shengq2/guoyan/20191104_cpd_analysis/fig_genome')

required.pkg <- c("data.table", "ggplot2")
for (pkg in required.pkg) {
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    require(pkg, character.only = TRUE)
  }
}

samples <- read.delim(inputFile, sep = '\t', header = T, stringsAsFactors = F)

if(nrow(samples) == 0) {
  stop(paste0("No data in ", inputFile))
}

if(!all(samples$Group %in% c(0, 1))){
  stop("Group should be either 0 for control or 1 for case in dinucleotide_file.list")
}

controlSamples<-samples$Sample[samples$Group==0]
caseSamples<-samples$Sample[samples$Group==1]
doBinomialTest<-(length(controlSamples) > 0) & (length(controlSamples > 0))

level.chr <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')

chrom_levels<-rev(level.chr)

info <- read.delim("chromInfo.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

colnames(info)<-c(paste0("V", c(1:ncol(info))))
info <- info[info$V1 %in% level.chr, c("V1", "V2")]
info <- data.frame(chrom = factor(info$V1, levels = chrom_levels),
                   chromStart = rep(0, nrow(info)),
                   chromEnd = info$V2)

ctype = "Site"
for (ctype in c("Site", "Read")){
  columnName = paste0(ctype, "Count")
  curPrefix = paste0(outfilePrefix, "_", columnName)
  pdf(paste0(curPrefix, ".pdf"), width=9, height=6, onefile = T)
  ## cases plot
  i<-1
  cnt.all <- NULL
  for (i in 1:nrow(samples)) {
    sampleFile <- samples$DinuFile[i]
    countFile <- samples$CountFile[i]
    sampleName <- samples$Sample[i]
    sampleGroup <- samples$Group[i]
    
    case <- fread(sampleFile, sep = "\t", header = T, data.table = F, stringsAsFactors = F)
    case <- case[case$Chrom %in% level.chr,]
    case$Chrom<-factor(case$Chrom, levels=chrom_levels)
    
    seg <- ggplot(data=info,aes(y = chrom, x = chromStart))+
      geom_segment(aes(y = chrom, yend = chrom, x = chromStart, xend = chromEnd),
                   lineend = "round", color = "Gainsboro", size = 5)+
      geom_segment(data = case,
                   aes_string(y = "Chrom", yend = "Chrom", x = "Start", xend = "End", color=columnName),
                   lineend = "butt", size = 5)+
      scale_y_discrete(drop=FALSE)+
      scale_color_gradient(high="DimGray",low="Gainsboro", name = paste0("Count per ", as.numeric(block) / 1000, " KB")) +
      scale_x_continuous(limits = c(0, 250e6),
                         expand = c(0, 3e6), ## expand dist = (maxlimit-minlimit)*a + b
                         breaks = seq(0, 250e6, by = 50e6),
                         labels = paste0(seq(0, 250, by = 50), "Mb"),
                         position = "top")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(10, 1, 1, 5, "pt"))+
      labs(x = "Genomic positions", y = NULL, title=paste0("Sample: ", sampleName))
    print(seg)
  }
  dev.off()
}