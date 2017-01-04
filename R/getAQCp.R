getAQCp <-
function(DList=NULL,colLabel=NULL,GList=NULL,overlapGenes=TRUE,filterGenes=F,cutRatioByMean=0.1, cutRatioByVar=0.1,
                    pvalAdjustGene=TRUE,pvalCutPath=0.05,pvalAdjustPath=TRUE,filterPathway=TRUE,minNumGenes=5,maxNumGenes=200){
  
  # DList - a list of studies
  # overlapGenes - get the overlap genes among studies
  # GList - a list of pathway
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  # topGene - use top number of DE genes instead of pvalue cutoff
  # topGeneNum - top number of genes if topGene==T
  # pvalCutGene - p value cutoff to define DE gene if topGene==F
  # pvalAdjustGene - adjust p-value or not for DE gene if topGene==F
  # pvalCutPath - p value cutoff to define enriched pathway
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
  
  # calculate AQCp
  AQCp=rep(0,times=K)
  for(k in 1:K){
    print(k)
    
    # single
    pSingle=DEsingle(DList[[k]],colLabel[[k]],pvalAdjustGene,method="t")
    pPathSingle=getPathPKS(totalGene,pSingle,GList,pvalAdjustPath)
    indSingle=which(pPathSingle<=pvalCutPath)
    
    # meta
    pMeta=DEmeta(DList[-k],colLabel[-k],pvalAdjustGene,method="t")
    pPathMeta=getPathPKS(totalGene,pMeta,GList,pvalAdjustPath)
    indMeta=which(pPathMeta<=pvalCutPath)
    
    # test
    num11=length(intersect(indSingle,indMeta)) # in single, in meta
    num12=length(indSingle)-num11  # in single , out meta
    num21=length(indMeta)-num11  # out single, in meta
    num22=nrow(DList[[k]])-length(indMeta)-num12  # out single, out meta
    
    pAQC=fisher.test(matrix(c(num11,num12,num21,num22),2,2),alternative="greater")$p.value
    AQCp[k]=-log(pAQC)
  } # end k
  AQCp[AQCp==Inf]=410
  names(AQCp)=names(DList)
  return(AQCp) 
}
