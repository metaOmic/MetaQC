getSMR <-
function(scoreTable){
  # scoreTable - studies by measurements
  
  rankTable=sapply(1:ncol(scoreTable),function(x) rank(-scoreTable[,x]))
  SMR=apply(rankTable,1,sum)/ncol(scoreTable)
  names(SMR)=rownames(rankTable)
  return(SMR)
}
