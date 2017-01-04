getCQCp <-
function(DList=NULL,colLabel=NULL,GList=NULL,overlapGenes=TRUE,filterGenes=F,cutRatioByMean=0.1, cutRatioByVar=0.1,
                    filterPathway=TRUE,minNumGenes=5,maxNumGenes=200){
  
  # DList - a list of studies
  # overlapGenes - get the overlap genes among studies
  # GList - a list of pathway
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  # pvalCutGene - p value cutoff to define DE gene
  # pvalAdjustGene - adjust p-value or not for DE gene
  # pvalAdjustPath - adjust p-value or not for enriched pathway
  # filterPathway - wehther to filter pathway or not
  # minNumGenes - miminum number of genes inside a pathway
  # maxNumGenes - maxinum number of genes inside a pathway
  
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
  
  # check GList
  if(is.null(GList)){
    stop("Error: GList is a required argument.")
  }
  
  if(length(GList)<=5){
    stop("Error: Insufficient number of qualified pathways.")
  }
  
  # get overlap genes
  if(overlapGenes){
    DList=metaOverlap(DList)
  }
  
  # filter genes
  if(filterGenes){
    DList=metaFilterData(DList,cutRatioByMean, cutRatioByVar)
  }
  
  # pathway filtering : filter out small and large pathways by overlapping number 
  totalGene=rownames(DList[[1]])
  if(filterPathway){
    pathwayLen=sapply(GList,function(x) return(length(intersect(x,totalGene))))
    rmInd=which(pathwayLen<minNumGenes | pathwayLen>maxNumGenes)
    if(length(rmInd)==length(GList)){
      stop("Error: No qualified pathway left.")
    }
    if(length(rmInd)>0){
      GList=GList[-rmInd]
    }
  }
  
  # calculate CQCp
  G=nrow(DList[[1]])
  CQCp=rep(0,times=K)
  for(k in 1:K){
    print(k)
    
    # single
    pSingle=DEsingle(DList[[k]],colLabel[[k]],pvalAdjust=FALSE,method="t")
    #rankA=rank(pSingle,ties.method="first")
    pPathSingle=getPathPKS(totalGene,pSingle,GList,pvalAdjust=FALSE)
    rankSingle=rank(pPathSingle,ties.method="first")
    
    # meta
    pMeta=DEmeta(DList[-k],colLabel[-k],pvalAdjust=FALSE,method="t")
    #rankB=rank(pMeta,ties.method="first")
    pPathMeta=getPathPKS(totalGene,pMeta,GList,pvalAdjust=FALSE)
    rankMeta=rank(pPathMeta,ties.method="first")
    
    # test
    rho=1-6*sum((rankSingle-rankMeta)^2)/G/(G^2-1)  ## --> problematic here, maybe minus value
    if(rho<=-1){
      print("Warining, rho<=-1, invalide testing")
      CQCp[k]=0
    }else if(rho>=1){
      print("Warining, rho>=1, invalide testing")
      CQCp[k]=Inf
    }else{
      tValue=rho*sqrt((G-2)/(1-rho^2))
      #pValue=1-pt(tValue,df=(G-2))
      pValue=pt(-tValue,df=(G-2))
      CQCp[k]=-log(pValue)
    }
    
    #rhoAB=1-6*sum((rankA-rankB)^2)/G/(G^2-1)
    #tValueAB=rhoAB*sqrt((G-1)/(1-rhoAB^2))
    
  } # end k
  CQCp[which(CQCp==Inf)]=410
  names(CQCp)=names(DList)
  return(CQCp) 
}
