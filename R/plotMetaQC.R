plotMetaQC <-
function(scoreTable){
  
  rmInd=c()
  for(j in 1:ncol(scoreTable)){
    if(sd(scoreTable[,j])==0){
      rmInd=c(rmInd,j)
    }
  }
  if(length(rmInd)>0){
    scoreTable=scoreTable[,-rmInd]
  }

  data.prcomp=prcomp(scoreTable, scale=T)
  a=(diag(ncol(scoreTable))*sqrt(10))%*%data.prcomp$rotation
  
  xl=range(c(data.prcomp$x[,1],a[,1]))*1.2
  yl=range(c(data.prcomp$x[,2],a[,2]))*1.2
  
  # study scatter
  plot(data.prcomp$x[,1],data.prcomp$x[,2],pch=1,cex=2.5,lwd=2,
       xlab="First Principle Component",ylab="Second Principle Component",
       xlim=xl,ylim=yl)
  text(data.prcomp$x[,1],data.prcomp$x[,2],cex=1)
  
  # dash line
  abline(v=0, lty=2, lwd=2)
  abline(h=0, lty=2, lwd=2)
  
  # measurement arrow
  arrows(0,0,a[,1],a[,2])
  text(a[,1]*1.2,a[,2]*1.2,label=colnames(scoreTable))
  
}
