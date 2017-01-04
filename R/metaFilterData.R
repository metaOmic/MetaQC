metaFilterData <-
function (data,cutRatioByMean=0.4, cutRatioByVar=0.4){
  
  if(cutRatioByMean<0 | cutRatioByMean>1 | cutRatioByVar<0 | cutRatioByVar>1){
    warning("cutRatioByMean and cutRatioByVar have to range [0,1]. Set to be default value 0.4.")
  }
  
  meanKeepRate=1-cutRatioByMean
  sdKeepRate=1-cutRatioByVar
  
  
  # filtering
  #exclude rows with lower mean
  rowMean=sapply(data,rowMeans)
  meanRank=apply(-rowMean,2,rank)
  meanRankSum=apply(meanRank,1,sum)
  meanKeep=round(length(meanRankSum)*meanKeepRate)
  meanSumRank=rank(meanRankSum)
  Fdata=sapply(data, function(x) return(x[which(meanSumRank<meanKeep),]))
  
  #exclude rows with lower sd
  rowSD=sapply(Fdata,function(x) return(apply(x,1,sd)))
  sdRank=apply(-rowSD,2,rank)
  sdRankSum=apply(sdRank,1,sum)
  sdKeep=round(length(sdRankSum)*sdKeepRate)
  sdSumRank=rank(sdRankSum)
  Fdata=sapply(Fdata,function(x) return(x[which(sdSumRank<sdKeep),]))
  
  return(Fdata)
}
