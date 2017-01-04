\name{MetaQC-package}
\alias{MetaQC-package}
\alias{MetaQC}
\docType{package}
\title{
MetaQC R package for quality control on multiple genomics data
}
\description{
MetaQC R package for quality control on multiple genomics data
}
\details{
\tabular{ll}{
Package: \tab MetaQC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-12-20\cr
License: \tab License\cr
}
}
\author{
Silvia Liu <silvia.shuchang.liu@gmail.com>
}
\references{
Kang, Dongwan D., et al. MetaQC: objective quality control and inclusion/exclusion criteria for genomic meta-analysis. Nucleic acids research 40.2 (2012): e15-e15.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%~~ documentation directory ~~
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
data(prostate8)
data(pathway)
DList=prostate8$data
colLabel=prostate8$dataLabel
GList=pathway[[1]]
filterGenes=TRUE
cutRatioByMean=0.3
cutRatioByVar=0.3

#QCresult=MetaQC(DList, colLabel, GList,filterGenes,cutRatioByMean,cutRatioByVar)
#plotMetaQC(QCresult$scoreTable)
}
