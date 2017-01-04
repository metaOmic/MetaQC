getPathPKS <-
function(geneName,geneScore,GList,pvalAdjust){
  # geneName - full gene name
  # geneScore - gene score (e.g. p-value)
  # GList - pathway list
  # pvalAdjust - adjust p-value or not
  
  pathNum=length(GList)
  pValue=rep(0,times=pathNum)
  for(i in 1:pathNum){
    ind=geneName%in%GList[[i]]
    if(length(ind)==0){
      pValue[i]=1
    }
    pValue[i]=ks.test(x=geneScore[ind],y=geneScore[-ind],alternative="greater")$p.value
  }
  
  if(pvalAdjust==F){
    return(pValue)
  }else{
    p.adj=p.adjust(pValue,method="BH")
    return(p.adj)
  }
}
