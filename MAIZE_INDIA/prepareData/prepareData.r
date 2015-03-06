## Parameters ####################################################################
 inputFile<-'~/test/mars12_pheno_geno.txt' # change the path as needed
 hasHeader<-TRUE # does the file has variable names in the 1st row?
 idCol<-1 # which column contains the id of the lines
 nColSkip<-2 # how many columns contain information other than genotypes
 sep='\t' # column-separator (tab in this case)
 na.pheno='NA' # NA-code for phenotypes
 na.geno='?'   # NA-code for genotypes
 minMAF<-1/100 # minimum minor-allele frequency (used to edit markers)
 maxFreqNA<-1/10 # maximum frequency of missing genotypes (used to edit markers)
## end of parameters #############################################################

# Getting the folder name and the filename
 folder<-dirname(inputFile)
 filename<-basename(inputFile)

# stores current working directory (cwd) and sets the working folder to 'folder'
 cwd<-getwd()
 setwd(folder)
   DATA<-read.table(filename,header=hasHeader,sep=sep,
                   na.strings=na.pheno,stringsAsFactors=FALSE)
   nMrk<-ncol(DATA)-nColSkip
   nInd<-nrow(DATA)     
   X<-matrix(nrow=nInd,ncol=nMrk,NA)
   colnames(X)<-colnames(DATA)[-(1:nColSkip)]
   pheno<-DATA[,1:nColSkip]
   MAP<-data.frame(mkr=colnames(X),alleleZero=NA,alleleOne=NA,freqAlleleOne=NA,
                   maf=NA,freqNA=NA,stringsAsFactors=FALSE)
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
       if(length(alleles)>2){ 
         cat('Warning, marker ',colnames(X)[i],' has more than two alleles \n')
       }       
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
   # saving the data for future analyses
   save(X,MAP,pheno,file='genData.RData')   

   # appending allele one to mrk id
    colnames(X)<-apply(X=cbind(colnames(X),MAP$alleleOne),FUN=paste,collapse='_',MARGIN=1)
    rownames(MAP)<-colnames(X)

    # summaries and plots
      dim(X) # number of individuals and number of markers after edits
      sum(!is.na(pheno$Index)) # number of lines with phenotypic record
      plot(MAP$freqNA) # proportion of NAs by marker
      hist(MAP$maf,30) # histogram of maf
      hist(pheno$Index,30)  # histogram of phenotype
   # restore the working directory
   setwd(cwd) 
