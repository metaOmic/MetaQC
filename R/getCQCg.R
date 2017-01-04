getCQCg <-
function(DList=NULL,colLabel=NULL,overlapGenes=TRUE,filterGenes=F,cutRatioByMean=0.1, cutRatioByVar=0.1){
  
  # DList - a list of studies
  # overlapGenes - get the overlap genes among studies
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  # pvalAdjustGene - adjust p-value or not for DE gene
  
  # check DList 
  if(is.null(DList)){
    stop("Error: DList is a required argument.")
  }
  K=length(DList)
  if(K<3){
    stop("Error: Insufficient number of studies.")
  }
  
  # check colLabel
  if(is.null(colLabel)){
    stop("Error: colLabel is a required argument.")
  }
  if(length(colLabel)!=K){
    stop("Error: Unequal length of DList and colLabel")
  }
  for(k in 1:K){
    if(prod(colLabel[[k]]%in%c(0,1))==F){
      stop(paste("Error: colLabel for study ",k," has to be either 0 or 1 ."))
    }
  }
  
  # get overlap genes
  if(overlapGenes){
    DList=metaOverlap(DList)
  }
  
  # filter genes
  if(filterGenes){
    DList=metaFilterData(DList,cutRatioByMean, cutRatioByVar)
  }
  
  # calculate CQCg
  G=nrow(DList[[1]])
  CQCg=rep(0,times=K)
  for(k in 1:K){
    # single
    pSingle=DEsingle(DList[[k]],colLabel[[k]],pvalAdjust=FALSE,method="t")
    rankSingle=rank(pSingle)
    
    # meta
    pMeta=DEmeta(DList[-k],colLabel[-k],pvalAdjust=FALSE,method="t")
    rankMeta=rank(pMeta)
    
    # test
    rho=1-6*sum((rankSingle-rankMeta)^2)/G/(G^2-1)  ## --> problematic here, maybe minus value
    if(rho<=-1){
      print("Warining, rho<=-1, invalide testing")
      CQCp[k]=0
    }else if(rho>=1){
      print("Warining, rho>=1, invalide testing")
      CQCp[k]=0
    }else{
      tValue=rho*sqrt((G-2)/(1-rho^2))
      #pValue=1-pt(tValue,df=(G-2))
      pValue=pt(-tValue,df=(G-2))
      CQCg[k]=-log(pValue)
    }
    
  } # end k
  # correction
  CQCg[which(CQCg==Inf)]=410
  names(CQCg)=names(DList)
  return(CQCg) 
}
