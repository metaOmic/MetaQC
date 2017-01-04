getPcorr <-
function(DList){
  Pcorr=list()
  for(s in 1:length(DList)){
    #print(s)
    Pcorr[[s]]=cor(t(DList[[s]]))
    rownames(Pcorr[[s]])=rownames(DList[[s]])
    colnames(Pcorr[[s]])=rownames(DList[[s]])
  }
  #names(Pcorr)=names(DList)
  
  return(Pcorr)
}
