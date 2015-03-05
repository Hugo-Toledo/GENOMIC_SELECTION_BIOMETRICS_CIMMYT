## Parameters ############################################################################
 inputFile<-'~/Dropbox/CIMMYT/CASE_STUDIES/MAIZE_INDIA/raw_data/mars12_pheno_geno.txt'
 hasHeader<-TRUE
 idCol<-1
 nColSkip<-2
 sep='\t' #tab in this case
 na.pheno='NA'
 na.geno='?'
 
 minMAF<-1/100
 maxFreqNA<-1/10
## end of parameters #####################################################################


# Getting the folder name and the filename
 folder<-dirname(inputFile)
 filename<-basename(inputFile)
 
# stores current working directory (cwd) and sets the working folder to 'folder'
 cwd<-getwd()
 setwd(folder)
 	DATA<-read.table(filename,header=hasHeader,sep=sep,na.strings=na.pheno,stringsAsFactors=FALSE)
	nMrk<-ncol(DATA)-nColSkip
	nInd<-nrow(DATA)		
	X<-matrix(nrow=nInd,ncol=nMrk,NA)
	colnames(X)<-colnames(DATA)[-(1:nColSkip)]
	pheno<-DATA[,1:nColSkip]
	MAP<-data.frame(mkr=colnames(X),alleleZero=NA,alleleOne=NA,freqAlleleOne=NA,maf=NA,freqNA=NA,stringsAsFactors=FALSE)
	rownames(MAP)<-colnames(X)
	DATA<-as.matrix(DATA[,-c(1:nColSkip)])
	
	## Recoding
	 for(i in 1:nMrk){
		geno<-matrix(ncol=2,byrow=TRUE,data=unlist(strsplit(DATA[,i],split='_')))
		tmp<-(geno[,1]==na.geno)|(geno[,2]==na.geno)
		geno[tmp,]<-NA
		alleles<-table(geno)
		alleles<-names(alleles)[order(alleles,decreasing=TRUE)]
		MAP$alleleZero[i]<-alleles[1]
		if(length(alleles)>1){
			MAP$alleleOne[i]<-alleles[2]
			if(length(alleles)>2){ cat('Warning, marker ',colnames(X)[i],' has more than two alleles \n') }
		}
		X[,i]<-2-(geno[,1]==alleles[1])-(geno[,2]==alleles[1])		
	 }
	 tmp<-colMeans(X,na.rm=TRUE)/2
	 MAP$freqAlleleOne<-tmp
	 MAP$maf<-ifelse(tmp<=.5,tmp,1-tmp)
	 MAP$freqNA<-colMeans(is.na(X))	
	# end of recoding
	
	# edits
	 tmp<-(MAP$maf>minMAF)&(MAP$freqNA<=maxFreqNA)
	 MAP<-MAP[tmp,]
	 X<-X[,tmp]
	#
	
	# creadting genData object
	DATA<-genData(pheno=pheno,geno=X,map=MAP)
	save(DATA,file='genData.RData')
	
 # restore the working directory
 setwd(cwd)
 