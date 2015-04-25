# A function to read and recode GBS data (assuming input file with format commonly used by J. Poland's lab).
# returns a BGData object
# contact: gdeloscampos@gmail.com
readGBS<-function( fileIn,nColSkip,recode,header=TRUE,rsCol,sep,
                   class='matrix',n=NULL,p=NULL,na.strings='NA',returnData=TRUE,
                   verbose=TRUE,nChunks=NULL,vmode='byte',dimorder=1:2,
                   folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep='') 
                  ){
    ########
    # Reads and recode GBS data
    # Value: a BGData object with slots @map and @geno
    ########
    #setwd('~/Dropbox/wheat_gbs/genos/')

    #fileIn="G80_NonRedund2.csv"
    #nColSkip=1
    #header=TRUE
    #class='matrix'
    #n=NULL
    #p=NULL
    #na.strings='NA'
    #returnData=TRUE
    #verbose=TRUE
    #nChunks=NULL
    #vmode='byte'
    #dimorder=1:2
    #folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep='') 
    #sep=','
    #recode=FALSE
    #rsCol=1
    
    if(is.null(n)){
        tmp<-scan(gzfile(fileIn),nlines=1,what=character(),quiet=T,sep=sep)
        n<-length(tmp)-nColSkip
    }
    
    if(is.null(p)){
       tmp<-gzfile(fileIn,open='r')
       p <- 0
       while(length(scan(tmp,what=character(),sep=sep,nlines=1,quiet=TRUE))>0){
            p<-p+1
            if(verbose){
              cat('Determining number of markers (reading line ',p, ')\n')
            }
       }
       if(header){  p<-p-1}   
       close(tmp)
     }    
     if(header){
        tmp<-scan(fileIn,nlines=1,what=character(),quiet=TRUE,sep=sep)
        GID<-tmp[-(1:nColSkip)]
    }else{ GID<-paste0('ID_',1:n)}
    if(class=='matrix'){
        geno<-matrix(nrow=n,ncol=p)
    }else{
        geno<-new(class,nrow=n,ncol=p,vmode=vmode,folderOut=folderOut,
                  nChunks=nChunks,dimorder=dimorder)
    }
    rownames(geno)<-GID

    fileIn<-gzfile(fileIn,open='r')
    map<-matrix(nrow=p,ncol=nColSkip,NA)
    if(header){ tmp<-scan(fileIn,nlines=1,quiet=TRUE,what=character())
                colnames(map)<-tmp[1:nColSkip] 
              }

    for(i in 1:p){
        time<-proc.time()        
        if(recode){
	        x<-scan(fileIn,nlines=1,what=character(),quiet=TRUE,sep=sep,na.strings=na.strings)
	        alleleOne<-unlist(strsplit(x[allelesCol],split='/'))[1]
    	    map[i,]<-x[(1:nColSkip)]
        	x<-x[-(1:nColSkip)]
        	z<-rep(0,n)
        	z[x==alleleOne]<-2
        	z[x=='H']<-1            
        }else{
           map[i,]<-scan(fileIn,n=nColSkip,what=character(),quiet=TRUE,na.strings=na.strings,sep=sep)
           z<-scan(fileIn,n=n,what=integer(),quiet=TRUE,na.strings=na.strings,sep=sep)
        }
        geno[,i]<-z
        
        if(verbose){
            cat('Marker ',i,' (of ',p, '),',round(proc.time()[3]-time[3],3),
                ' sec./marker.\n',sep='')
        }
    }
    close(fileIn)

    # Adding names
    colnames(geno)<-map[,rsCol]
  
    pheno<-data.frame(GID=GID,stringsAsFactors=FALSE)
    MAP<-data.frame(type.convert(map[,1],as.is=TRUE),stringsAsFactors=FALSE)
    if(ncol(map)>2){
    	for(i in 2:ncol(map)){
    		MAP<-cbind(MAP,type.convert(map[,i],as.is=TRUE))
    	}
    }
    colnames(MAP)<-colnames(map)
 
    BGData<-new('BGData',geno=geno,pheno=pheno,map=MAP)
    if(class!='matrix'){
        attr(BGData,'origFile')<-list(path=fileIn,dataType=typeof(dataType))
        attr(BGData,'dateCreated')<-date()
        save(BGData,file=paste(folderOut,'/BGData.RData',sep=''))
    }
    if(returnData){
        return(BGData)
    }
} 
    
   # Functions to get summaries. 
   summary_num<-function(x){
	out=c(range(x,na.rm=TRUE),mean(is.na(x)),mean(x,na.rm=TRUE))
	names(out)<-c('min','max','freq-NA','mean')
	return(out)
   }
   summarySNPs<-function(X){
	TMP<-t(apply(FUN=summary_num,MARGIN=2,X=X))
	return(TMP)
  }

###
