# singleSample(): When only one (group of) sample is specified, performs Chisq test for overall 5 categories and one-vs-other test for four DINUCs.
# INPUT smp: either file name or read-in data.frame for single sample. (When reg is specified, smp must be preprocessed to contain counts for the corresponding region.)
# INPUT gn: string for genome: hg38,hg19,saccer3
# INPUT reg: specifying the single region.
### choice 1 can be a keyword {'intron','exon','promoter'} etc.
### choice 2 can be a underscore-separated vector of four big integers (or proportions), freshly counted from corresponding CUSTOM region from refGenome.
# INPUT cnt5file='gn_reg_cnt5.tsv': one data file pre-storing count vectors for common regions of a few genomes we accommodate.
# OUTPUT res$res: 5x6 test result table in character matrix.
# OUTPUT res$unit: inform on the unit of sample counts used in statistical tests & presentation.
singleSample <- function(smp,gn,reg=NULL,cnt5file='gn_reg_cnt5.tsv0',nDigit=4,unit2=c(1e6,1e3)) { #nDigit: number of significant digits.
	if (!require(reshape2)) stop('Please install library reshape2 for dcast function')
	gn <- tolower(gn)
	DINUC4 <- c('TT','TC','CC','CT')
	if (is.null(reg)) { # Simplest case: no specification of region
		p4gn <- switch(gn,
			hg38=c(0.0986,0.0595,0.0515,0.0700), #0.7204
			hg19=c(0.0980,0.0594,0.0521,0.0700), #0.7206
			saccer3=rep(1/16,4) # place-holder p vector
		) # Limited to considered few refGenomes.
	} else { # parse out p4 vector when reg is supplied.
		p4gn <- parseReg(reg,gn,cnt5file=cnt5file)[1:4]
	} # expected probablility vector for five categories.
  if (!is.null(p4gn)) { # Somehow we might fail retrieving expected p-vector
    p5gn <- c(p4gn,1-sum(p4gn))
    names(p5gn) <- c(DINUC4,'others')
		if (!grepl(',',smp)) {
			smpCnt <- collapse.cntFile(smp,chrs=NULL)
		} else {
			smpCnt <- collapseGrp2Smp(smp)
		}
		if (!is.null(unit2)) {
			res <- convert.cntUnit(smpCnt,unit2)
			smpCnt <- res$cnt
			unit <- res$unit
		} else {
			unit <- 1
		}
    ########### Overall test spanning 5 categories ###########
		all.chisq.res <- chisq.test(smpCnt,p=p5gn)
    ########### Four iterations of ith vs. others ############
		oneVSother <- matrix(NA,nr=4,nc=6,dimnames=list(DINUC4,c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj')))
		for (di in DINUC4) {
			x <- smpCnt[di]
			n <- x+smpCnt['others']
			oneVSother[di,'exp'] <- exp <- signif(p5gn[di]/(p5gn[di]+p5gn['others']),nDigit)
			oneVSother[di,'obs'] <- obs <- signif(x/n,nDigit)
			p.binom <- pbinom(x-1,n,exp,lower.tail=F)
			oneVSother[di,'p.exact'] <- signif(max(p.binom,1e-15),nDigit)
			p.chisq <- chisq.test(c(x,smpCnt['others']),p=c(exp,1-exp))$p.val
			oneVSother[di,'p.chisq'] <- signif(max(p.chisq,1e-15),nDigit)
			oneVSother[di,'p.exact.adj'] <- min(signif(p.binom*length(DINUC4),nDigit),1)
			oneVSother[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		}
		expCnt <- floor(all.chisq.res$exp)
		smpCnt <- floor(smpCnt)
		expCnt <- paste(expCnt,collapse=':')
		smpCnt <- paste(smpCnt,collapse=':')
		all.row <- c(expCnt,smpCnt,'---','---',signif(all.chisq.res$p.val,nDigit),'---')
		res <- rbind(oneVSother,Overall=all.row)
  } else {
		warning('Cannot retrieve reference dinuc distribution vector for ',gn,' ',reg,'!!!')
		res <- matrix(NA,nr=5,nc=5,dimnames=list(DINUC4,c('exp','obs','95%CI','p.binom','p.chisq')))
  }
	return(list(res=res,unit=unit))
}

# twoReg(): Compare dinuc counts between two region sets in a single sample while taking into account reference dinuc distribution of the two regions in refGenome.
# INPUT smp1: filename for region1-limited sample cnt file.
# INPUT smp2: filename for region2-limited sample cnt file.
# INPUT reg1 & reg2: to inform reference p5 vector. A region keyword or a count5 string (comma-separated).
# NOTE smp1 & smp2: when supplied as multi-samples, smp1 & smp2 should be two groups of equal lengths. Haven't enforced such a pre-check, but it should be kept in mind.
twoReg <- function(smp1,smp2,reg1,reg2,gn,cnt5file='gn_reg_cnt5.tsv',nDigit=4,unit2=c(1e6,1e3),DINUC4=c('TT','TC','CC','CT')) {
	if (!grepl(',',smp1)) {
		smpCnt1 <- collapse.cntFile(smp1,chrs=NULL)
	} else {
		smpCnt1 <- collapseGrp2Smp(smp1)
	}
	if (!grepl(',',smp2)) {
		smpCnt2 <- collapse.cntFile(smp2,chrs=NULL)
	} else {
		smpCnt2 <- collapseGrp2Smp(smp2)
	}
	if (!is.null(unit2)) {
		res <- convert.cntUnit(rbind(smpCnt1,smpCnt2),unit2)
		smpCnt1 <- res$cnt[1,]
		smpCnt2 <- res$cnt[2,]
		unit <- res$unit
	} else {
		unit <- 1
	}
	ref5.1 <- parseReg(reg1,gn,scaleNow=F,cnt5file=cnt5file)
	ref5.2 <- parseReg(reg2,gn,scaleNow=F,cnt5file=cnt5file)
	res <- matrix(NA,nr=5,nc=7,dimnames=list(c(DINUC4,'Overall'),c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj','1stVS2nd')))
	chisq.stats <- numeric(length(DINUC4))
	names(chisq.stats) <- DINUC4
	for (di in DINUC4) {
		nRef.di <- c(ref5.1[di],ref5.2[di])
		p.di <- nRef.di/sum(nRef.di)
		obs.di <- c(smpCnt1[di],smpCnt2[di])
		chisq.res <- chisq.test(x=obs.di,p=p.di)
		p.chisq <- chisq.res$p.val
		chisq.stats[di] <- chisq.res$stat
		res[di,'exp'] <- paste(floor(nRef.di/unit),collapse=':') # Showing ratio (a:b)
		res[di,'obs'] <- paste(floor(obs.di),collapse=':') # Showing ratio (a:b)
		p.binom <- pbinom(obs.di[1]-1,sum(obs.di),p.di[1],lower.tail=F)
		res[di,'p.chisq'] <- signif(max(p.chisq,1e-15),nDigit)
		res[di,'p.exact'] <- signif(max(p.binom,1e-15),nDigit)
		res[di,'p.exact.adj'] <- min(signif(p.binom*length(DINUC4),nDigit),1)
		res[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		res[di,'1stVS2nd'] <- ifelse((obs.di[1]-1)/sum(obs.di)>p.di[1],'>','<')
	}
	p.exact.all <- sumPbyFisher(as.numeric(res[DINUC4,'p.exact']))
	res['Overall','p.exact'] <- signif(max(p.exact.all,1e-15),nDigit)
	p.chi2.all <- sumPbyChi2(chisq.stats)
	res['Overall','p.chisq'] <- signif(max(p.chi2.all,1e-15),nDigit)
	res['Overall',setdiff(colnames(res),c('p.exact','p.chisq'))] <- '---'
	return(list(res=res,unit=unit))
}
# twoGrp(): Given two sample file names (.cnt), return 5*7 statistical test result table.
# INPUT smp1: either file name or read-in data.frame for sample 1.
# INPUT smp2: either file name or read-in data.frame for sample 2.
# OUTPUT res$res: 5x7 test result table.
### Five rows (DINUC4+Overall), 7 columns (exp,obs,p.exact,p.exact.adj,p.chisq,p.chisq.adj,1stVS2nd)
# NOTE only chisquared test is employed for both overall and one-vs-others.
# NOTE number of digits printed out in kable is adjustable via nDigit
twoGrp <- function(smp1,smp2,nDigit=4,DINUC4=c('TT','TC','CC','CT'),unit2=c(1e6,1e3)) {
	if (!grepl(',',smp1)) {
		smpCnt1 <- collapse.cntFile(smp1,chrs=NULL)
	} else {
		smpCnt1 <- collapseGrp2Smp(smp1)
	}
	if (!grepl(',',smp2)) {
		smpCnt2 <- collapse.cntFile(smp2,chrs=NULL)
	} else {
		smpCnt2 <- collapseGrp2Smp(smp2)
	}
	if (!is.null(unit2)) {
		res <- convert.cntUnit(rbind(smpCnt1,smpCnt2),unit2)
		smpCnt1 <- res$cnt[1,]
		smpCnt2 <- res$cnt[2,]
		unit <- res$unit
	} else {
		unit <- 1
	}
	res <- matrix(NA,nr=5,nc=7,dimnames=list(c(DINUC4,'Overall'),c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj','1stVS2nd')))
	for (di in DINUC4) {
		N2.1 <- c(smpCnt1[di],smpCnt1['others']) # two numbers of grp1
		N2.2 <- c(smpCnt2[di],smpCnt2['others'])# two numbers of grp2
		chisq.res <- chisq.test(rbind(N2.1,N2.2))
		exp1 <- floor(chisq.res$expected[1,1])
		exp2 <- floor(chisq.res$expected[2,1])
		res[di,'exp'] <- paste(exp1,exp2,sep=', ')
		res[di,'obs'] <- paste(floor(N2.1[1]),floor(N2.2[1]),sep=', ') # Observed di cnts in millions
		p.chisq <- max(chisq.res$p.val,1e-15) # When too little p (say 0) occurs, upgrade to 1e-15.
		res[di,'p.chisq'] <- signif(p.chisq,nDigit)
		res[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		m=smpCnt1[di]+smpCnt2[di]
		n=smpCnt1['others']+smpCnt2['others']
		k=sum(N2.1)
		x=smpCnt1[di]
		p.hyper1 <- phyper(x-1,m,n,k,lower.tail=F)
		p.hyper2 <- phyper(x,m,n,k)
		p.hyper <- max(min(p.hyper1,p.hyper2),1e-15)
		res[di,'p.exact'] <- signif(p.hyper,nDigit)
		res[di,'p.exact.adj'] <- min(signif(p.hyper*length(DINUC4),nDigit),1)
		res[di,'1stVS2nd'] <- ifelse(N2.1[1]>=chisq.res$expected[1,1],'>','<')
	}
	all.chisq.res <- chisq.test(rbind(smpCnt1,smpCnt2))
	obs1 <- paste(floor(smpCnt1),collapse=':')
	obs2 <- paste(floor(smpCnt2),collapse=':')
	exp.M <- floor(all.chisq.res$exp)
	exp1 <- paste(exp.M[1,],collapse=':')
	exp2 <- paste(exp.M[2,],collapse=':')
	res['Overall','exp'] <- paste(exp1,exp2,collapse='; ')
	res['Overall','obs'] <- paste(obs1,obs2,collapse='; ')
	res['Overall','p.exact'] <- res['Overall','p.exact.adj'] <- res['Overall','p.chisq.adj'] <- '---'
	res['Overall','p.chisq'] <- signif(max(all.chisq.res$p.val,1e-15),nDigit)
	res['Overall','1stVS2nd'] <- ifelse(smpCnt1[1]>=all.chisq.res$exp[1,1],'>','<')
	return(list(res=res,unit=unit))
}

# collapse.cntFile(): Given one cnt filename, derive smpCnt (five cnts) vector.
# INPUT: only one essential input - sample file name (.cnt)
collapse.cntFile <- function(smp,chrs=NULL,DINUC4=c('TT','TC','CC','CT')) {
	if (!require(reshape2))
		stop('Please install library reshape2 for dcast function')
	if (length(as.vector(smp))==1)
		smp <- read.delim(smp,as.is=T)
	if (!is.null(chrs))
		smp <- smp[smp[,'Chromosome']%in%chrs,]
	dinuc16 <- genPerm2()
	smp <- smp[smp[,'Dinucleotide']%in%dinuc16,]
	smp <- dcast(smp,Dinucleotide~Chromosome,value.var='Count')
	smpMat <- as.matrix(smp[,-1])
	rownames(smpMat) <- smp$Dinucleotide
	colnames(smpMat) <- colnames(smp)[-1]
	smpCnt0 <- apply(smpMat,1,function(x) sum(as.numeric(x),na.rm=T))
	smpCnt.4 <- smpCnt0[DINUC4]
	smpCnt.5th <- sum(as.numeric(smpCnt0[setdiff(names(smpCnt0),DINUC4)]))
	smpCnt <- c(smpCnt.4,others=smpCnt.5th)
}
# genPerm2(): # generate permutations of two-tuples.
genPerm2 <- function(vocab=c('A','T','G','C')) {
	comb.2col <- t(combn(vocab,2))
	comb.col2 <- comb.2col[,2:1]# reverse 1st col & 2nd col of comb.2col
	sameRepeat <- cbind(vocab,vocab)
	prototype <- data.frame(rbind(comb.2col,comb.col2,sameRepeat))
	perm2 <- do.call(paste0,prototype)
	perm2
}
# parseReg(): parse reg to retrieve p5 vector. Typically extraction of the first 4 values will immediately follow.
# INPUT reg: either one keyword or five underscore separated frequency values.
# OUTPUT p5gn: a 4-component frequency vector. Program can abort at unexpected input.
# INPUT scaleNow: logic variable instructing if an across-5-dinuc scaling should be performed at the current step. set to TRUE by default.
parseReg <- function(reg,gn='hg38',scaleNow=T,cnt5file='gn_reg_cnt5.tsv',DINUC4=c('TT','TC','CC','CT')) {
	reg <- separated <- unlist(strsplit(reg,','))
	if (length(separated)==5) { # Underscore separated into 5 numeric strings.
		separated <- as.numeric(separated)
		p5gn <- separated
		if (max(p5gn)>1 & scaleNow) p5gn <- separated/sum(as.numeric(separated))
	} else if (length(separated)==1) {
		p5mat <- read.delim(cnt5file,row.names=1)
		if (is.element(paste(gn,reg,sep='_'),rownames(p5mat)))
			p5gn <- as.numeric(p5mat[paste(gn,reg,sep='_'),])
			if (scaleNow) p5gn <- p5gn/sum(p5gn)
	} else {
		stop('Please supply one keyword or five underscore-separated frequency values for one region!!!')
	}
	names(p5gn) <- c(DINUC4,'others')
	p5gn
}

# INPUT chi2.dinucs: chisq test statistics across two regions for each dinuc.
sumPbyChi2 <- function(chi2.dinucs) {
	sumchi2 <- sum(chi2.dinucs)
	dof <- length(chi2.dinucs)
	p <- pchisq(sumchi2,dof,lower.tail=F)
	p
}
# INPUT pvals: individual test p-values across two regions for each dinuc.
sumPbyFisher <- function(pvals,min.th=1e-147) {
	pvals[pvals<min.th] <- min.th
	aggChi2 <- (-2*sum(log(pvals)))
	dof <- 2*length(pvals)
	p <- pchisq(aggChi2,dof,lower.tail=F)
	p
}

# collapseGrp2Smp(): when sample CNT file name implies multiple files are supplied, flattenGrp2Smp() is invoked to collapse a group to a sample.
# NOTE: collapseGrp2Smp() calls on collapse.cntFile() to yield sample-wise 5-tuple count vectors, and then sum up counts by dinuc type.
# INPUT grp: a comma separated string, implying multiple CNT files.
# OUTPUT CNT5: a 5-tuple of counts, named DINUC4 + others.
collapseGrp2Smp <- function(grp) {
	samples <- unlist(strsplit(grp,','))
	cnt.samples <- matrix(nr=length(samples),nc=5)
	for (i in 1:length(samples)) {
		sample <- samples[i]
		cnt5 <- collapse.cntFile(sample)
		cnt.samples[i,] <- cnt5
	}
	CNT5 <- apply(cnt.samples,2,function(x) sum(x,na.rm=T))
	names(CNT5) <- names(cnt5)
	CNT5
}
# convert.cntUnit(): Convert dinuc count to bigger unit, say a million.
# INPUT cnt0: raw counts, most likely as a 5-value vector.
# INPUT unit2: two candidate count units, one bigger (say 1m) one smaller (say 1k). 
# OUTPUT res: list of two components
###### $cnt - counts converted to bigger unit.
###### $unit - the actual unit used in converted counts.
convert.cntUnit <- function(cnt0,unit2=c(1e6,1e3)) {
	uH <- max(unit2)
	uL <- min(unit2)
	cnt <- floor(cnt0/uH)
	unit <- uH
	if (min(cnt)<1) {
		cnt <- floor(cnt0/uL)
		unit <- uL
	}
	res <- list(cnt=cnt,unit=unit)
	res
}




