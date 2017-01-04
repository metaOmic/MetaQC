getIQC <-
function(DList=NULL,overlapGenes=TRUE,filterGenes=TRUE,cutRatioByMean=0.1, cutRatioByVar=0.1){
  # DList - a list of studies
  # overlapGenes - get the overlap genes among studies
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  
  # check DList 
  if(is.null(DList)){
    stop("Error: Dlist is a required argument.")
  }
  K=length(DList)
  if(K<3){
    stop("Error: Insufficient number of studies.")
  }
  
  # get overlap genes
  if(overlapGenes){
    DList=metaOverlap(DList)
  }
  
  # filter genes
  if(filterGenes==TRUE){
    DList=metaFilterData(DList,cutRatioByMean, cutRatioByVar)
  }
  
  # get Pearson correlation
  Pcorr=getPcorr(DList)
  
  # List form
  PcorrVec=list()
  for(k in 1:K){
    PcorrVec[[k]]=Pcorr[[k]][upper.tri(Pcorr[[k]],diag=T)]
  }
  
  # get dissimilarity between studies
  d=matrix(1,K,K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      d[i,j]=(1-cor(x=PcorrVec[[i]],y=PcorrVec[[j]],method="spearman"))/2
      d[j,i]=d[i,j]
    }
  }
  
  # wilcoxon test
  pValue=rep(0,times=K)
  for(k in 1:K){
    Dstar=d[-k,k]
    tmp=d[-k,-k]
    Dpound=tmp[upper.tri(tmp)]
    pValue[k]=wilcox.test(x=Dstar,y=Dpound,alternative="greater")$p.value
  }
  
  # tranform and IQC
  z95 = qnorm(0.95,0,1)
  gp=1-pnorm(qnorm(pValue,mean=z95,sd=1),mean=-z95,sd=1)
  IQC=-log(gp,base=10)
  names(IQC)=names(DList)
  
  IQC[IQC==Inf]=410
  return(IQC)
}
