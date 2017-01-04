getEQC <-
function(DList,GList,overlapGenes=TRUE,filterGenes=TRUE,cutRatioByMean=0.1, cutRatioByVar=0.1, 
                   filterPathway=TRUE,minNumGenes=5,maxNumGenes=200, B=1000){
  # DList - a list of studies
  # GList - a list of pathway
  # overlapGenes - get the overlap genes among studies
  # filterGenes - meta filter gene or not
  # cutRatioByMean - meta mean filter ratio
  # cutRatioByVar - meta var filter ratio
  # filterPathway - whether to filter pathway or not
  # minNumGenes - miminum number of genes inside a pathway
  # maxNumGenes - maxinum number of genes inside a pathway
  # B - permutation times
  
  # check DList 
  if(is.null(DList)){
    stop("Error: DList is a required argument.")
  }
  K=length(DList)
  if(K<=4){
    stop("Error: Insufficient number of studies.")
  }
  
  # check GList
  if(is.null(GList)){
    stop("Error: GList is a required argument.")
  }
  
  # get overlap genes
  if(overlapGenes){
    DList=metaOverlap(DList)
  }
  
  # filter genes
  if(filterGenes){
    DList=metaFilterData(DList,cutRatioByMean, cutRatioByVar)
  }
  
  # get Pearson correlation
  Pcorr=getPcorr(DList)
  
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
  
  EQC=rep(0,times=K)   # recording EQC value
  for(k in 1:K){
    print(k)
    
    geneName=rownames(DList[[k]])
    pathway=lapply(GList,function(x) return(intersect(geneName,x)))  # intersect genes
    M=length(pathway)
    mLen=sapply(pathway,length)
    
    # get the denominator
    Gnum=nrow(DList[[k]])
    rho=Pcorr[[k]][upper.tri(Pcorr[[k]])]
    denominator=sqrt(sum(rho^2)/(Gnum *(Gnum-1)/2))
    
    # calculate tNull, based on pathway length
    tNull=list()
    mLenUniq=sort(unique(mLen))
    mLenUniq=mLenUniq[which(mLenUniq!=0)]
    for(m in mLenUniq){
      tNull[[m]]=rep(0,times=B)
      for(b in 1:B){
        ind=sample(1:Gnum,m)
        #print(c(b,ind))
        tmpCorr=Pcorr[[k]][ind,ind]
        rho=tmpCorr[upper.tri(tmpCorr)]
        numerator=sqrt(sum(rho^2)/(m*(m-1)/2))
        tNull[[m]][b]=numerator/denominator
      }
    }
    
    # calculate the p-value for tNull
    tNullp=list()
    for(m in mLenUniq){
      tNullp[[m]]=floor(rank(tNull[[m]]))/(B+1)
    }
    
    # calculate tValue
    tValue=rep(0,times=M)
    for(m in 1:M){
      # tValue
      ind=rownames(DList[[k]])%in%pathway[[m]]
      tmpCorr=Pcorr[[k]][ind,ind]
      rho=tmpCorr[upper.tri(tmpCorr)]
      numerator=sqrt(sum(rho^2)/(mLen[m]*(mLen[m]-1)/2))
      tValue[m]=numerator/denominator      
    }
    
    # calculate pValue
    pValue=sapply(1:M, function(x) return( (sum(tNull[[mLen[x]]]>tValue[x])+1)/(B+1) ))
    
    # get Sk value
    Sk= -2 * sum(log(pValue))
    
    # get SkNull value
    SkNull=rep(0,times=B)
    for(b in 1:B){
      nullP=rep(0,times=M)
      for(m in 1:M){
        nullP[m]=sample(tNullp[[mLen[m]]],1)
      }
      SkNull[b]=-2 *sum(log(nullP))
    }
    
    pEQC=(sum(SkNull>Sk)+1)/(B+1)
    EQC[k]=-log(pEQC)
    print(EQC[k])
    
  }  # end k
  
  names(EQC)=names(DList)
  return(EQC)
  
}
