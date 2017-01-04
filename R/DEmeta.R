DEmeta <-
function(data,colLabel,pvalAdjust=TRUE,method="t"){
  # data - list of data matrix
  # colLabel - list of column label for data
  # pvalAdjust - adjust p-value or not
  # DE method - c("t","KS")
  
  S=length(data)
  pValue=matrix(0,nrow(data[[1]]),S)
  for(s in 1:S){
    pValue[,s]=DEsingle(data[[s]],colLabel[[s]],pvalAdjust=F,method)
  }
  
  Fisher=-2 * apply(log(pValue),1,sum)
  p.Fisher=pchisq(Fisher,df=2*S,lower.tail=FALSE)
  
  if(pvalAdjust==FALSE){
    return(p.Fisher)
  }else{
    p.adj=p.adjust(p.Fisher,method="BH")
    return(p.adj)
  }
  
}
