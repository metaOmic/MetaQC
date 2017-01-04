MetaQC <-
function(DList=NULL,
                 colLabel=NULL,
                 GList=NULL,
                 #overlapGenes=TRUE,
                 filterGenes=FALSE,
                 cutRatioByMean=0.3,
                 cutRatioByVar=0.3,
                 pvalCutGene=0.05,
                 pvalAdjustGene=TRUE,
                 pvalCutPath=0.05,
                 pvalAdjustPath=TRUE,
                 #filterPathway=TRUE,
                 minNumGenes=5,
                 maxNumGenes=200,
                 B=100){
  
  overlapGenes=TRUE
  filterPathway=TRUE
  
  print("Begin IQC")
  IQC=getIQC(DList,overlapGenes,filterGenes,cutRatioByMean, cutRatioByVar)
  print(IQC)
  print("End IQC")
  
  print("Begin EQC")
  print("Warning: EQC to be updated using some temporary values")
  EQC=abs(IQC+rnorm(length(IQC),sd=0.5))
  print("End EQC")
    
  print("Begin AQCg")
  AQCg=getAQCg(DList,colLabel,overlapGenes,filterGenes,cutRatioByMean, cutRatioByVar,pvalCutGene, pvalAdjustGene)
  print(AQCg)
  print("End AQCg")
    
  print("Begin AQCp")
  AQCp=getAQCp(DList,colLabel,GList,overlapGenes,filterGenes,cutRatioByMean, cutRatioByVar,
               pvalAdjustGene,pvalCutPath,pvalAdjustPath,filterPathway,minNumGenes,maxNumGenes)
  print(AQCp)
  print("End AQCg")
  
  print("Begin CQCg")
  CQCg=getCQCg(DList,colLabel,overlapGenes,filterGenes,cutRatioByMean, cutRatioByVar)
  print("End CQCg")
  print("End CQCg")

  print("Begin CQCp")
  CQCp=getCQCp(DList,colLabel,GList,overlapGenes,filterGenes,cutRatioByMean, cutRatioByVar,
               filterPathway,minNumGenes,maxNumGenes)
  print(CQCp)
  print("End CQCp")
   
  scoreTable=cbind(IQC,EQC,AQCg,AQCp,CQCg,CQCp)
  rownames(scoreTable)=names(DList)
  SMR=getSMR(scoreTable)
  
  return(list(scoreTable=scoreTable, SMR=SMR))
}
