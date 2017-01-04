getAQCg <-
function(DList=NULL,colLabel=NULL,overlapGenes=TRUE,filterGenes=F,cutRatioByMean=0.1, cutRatioByVar=0.1,
                    pvalCutGene=.05, pvalAdjustGene=TRUE){
  
  # DList - a list of studies
  # overlapGenes - get the overlap genes among studies
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  # pvalCutGene - p value cutoff to define DE gene
  # pvalAdjustGene - adjust p-value or not
  
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
  
  # calculate AQCg
  AQCg=rep(0,times=K)
  for(k in 1:K){
    # single
    pSingle=DEsingle(DList[[k]],colLabel[[k]],pvalAdjustGene,method="t")
    indSingle=which(pSingle<=pvalCutGene)
    
    # meta
    pMeta=DEmeta(DList[-k],colLabel[-k],pvalAdjustGene,method="t")
    indMeta=which(pMeta<=pvalCutGene)
    
    num11=length(intersect(indSingle,indMeta)) # in single, in meta
    num12=length(indSingle)-num11  # in single , out meta
    num21=length(indMeta)-num11  # out single, in meta
    num22=nrow(DList[[k]])-length(indMeta)-num12  # out single, out meta
    
    pAQC=fisher.test(matrix(c(num11,num12,num21,num22),2,2),alternative="greater")$p.value
    AQCg[k]=-log(pAQC)
  } # end k
  AQCg[AQCg==Inf]=410
  names(AQCg)=names(DList)
  return(AQCg) 
}
