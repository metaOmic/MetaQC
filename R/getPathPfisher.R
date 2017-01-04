getPathPfisher <-
function(DEgene,totalGene,GList,pvalAdjust){
  # method - "fisher" or "KS
  pathNum=length(GList)
  pValue=rep(0,times=pathNum)
  for(i in 1:pathNum){
    n11=length(intersect(DEgene,GList[[i]]))  # in DE, in path
    n12=length(DEgene)-n11  # in DE, out path
    n21=length(GList[[i]])-n11  # out DE, in path
    n22=length(totalGene)-length(DEgene)-n21  # out DE, out path
    pValue[i]=fisher.test(matrix(c(n11,n12,n21,n22),2,2),alternative="greater")$p.value
  }  
  if(pvalAdjust==F){
    return(pValue)
  }else{
    p.adj=p.adjust(pValue,method="BH")
    return(p.adj)
  }
}
