DEsingle <-
function(data,colLabel,pvalAdjust=TRUE,method="t"){
  # data - data matrix
  # colLabel - column label for data
  # pvalAdjust - adjust p-value or not
  # DE method - c("t","KS")
  
  if(method=="KS"){
    ind0=which(colLabel==0)
    ind1=which(colLabel==1)
    pValue=sapply(1:nrow(data),function(i) return(ks.test(x=data[i,ind0],y=data[i,ind1])$p.value))
  }else{ # t test
    t.stat=mt.teststat(data,colLabel,test="t.equalvar")
    pValue=2*pt(abs(t.stat),df=(length(colLabel)-2),lower.tail=FALSE)
  }
  
  if(pvalAdjust==FALSE){
    return(pValue)
  }else{
    p.adj=p.adjust(pValue,method="BH")
    return(p.adj)
  }
}
