metaOverlap <-
function(data){
  S=length(data)
  
  for(s in 1:S){
    rownames(data[[s]])=toupper(rownames(data[[s]]))
  }
  
  interG=rownames(data[[1]])
  for(s in 2:S){
    interG=intersect(interG,rownames(data[[s]]))
  }
  if(length(interG)<=5){
    stop(paste("Error: Number of overlap genes: ",interG,", less than 5.",sep=""))
  }
  for(s in 1:S){
    data[[s]]=data[[s]][interG,]
  }
  
  return(data)
}
