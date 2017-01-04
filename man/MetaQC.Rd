\name{MetaQC}
\alias{MetaQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta QC main function to calculate score table and SMR value
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MetaQC(DList = NULL, colLabel = NULL, GList = NULL, filterGenes = FALSE, 
  cutRatioByMean = 0.3, cutRatioByVar = 0.3, pvalCutGene = 0.05, 
  pvalAdjustGene = TRUE, pvalCutPath = 0.05, pvalAdjustPath = TRUE, 
  minNumGenes = 5, maxNumGenes = 200, B = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{
  DList}{a list of study matrix
%%     ~~Describe \code{DList} here~~
}
  \item{colLabel}{
  a list of column label for each study
%%     ~~Describe \code{colLabel} here~~
}
  \item{GList}{
  a list for gene sets (pathway)
%%     ~~Describe \code{GList} here~~
}
  \item{filterGenes}{
  a bool variable to filter genes or not
%%     ~~Describe \code{filterGenes} here~~
}
  \item{cutRatioByMean}{
  if filterGenes==TRUE, the percentage of filtering genes by low mean value
%%     ~~Describe \code{cutRatioByMean} here~~
}
  \item{cutRatioByVar}{
  if filterGenes==TRUE, the percentage of filtering genes by low variance value
%%     ~~Describe \code{cutRatioByVar} here~~
}
  \item{pvalCutGene}{
  p-value cutoff for gene to decide significant DE genes
%%     ~~Describe \code{pvalCutGene} here~~
}
  \item{pvalAdjustGene}{
  a bool variable to adjust gene p-value or not
%%     ~~Describe \code{pvalAdjustGene} here~~
}
  \item{pvalCutPath}{
  p-value cutoff for gene set enrichment analysis to decide significant pathway 
%%     ~~Describe \code{pvalCutPath} here~~
}
  \item{pvalAdjustPath}{
  a bool variable to adjust pathway p-value or not
%%     ~~Describe \code{pvalAdjustPath} here~~
}
  \item{minNumGenes}{
  minimum number of genes inside a pathway
%%     ~~Describe \code{minNumGenes} here~~
}
  \item{maxNumGenes}{
  maximum number of genes inside a pathway
%%     ~~Describe \code{maxNumGenes} here~~
}
  \item{B}{
  permutation times for EQC
%%     ~~Describe \code{B} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{scoreTable}{scoreTable for each measurement}
\item{SMR}{SMR value for each study}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kang, Dongwan D., et al. MetaQC: objective quality control and inclusion/exclusion criteria for genomic meta-analysis. Nucleic acids research 40.2 (2012): e15-e15.
%% ~put references to the literature/web site here ~
}
\author{
Silvia Liu <silvia.shuchang.liu@gmail.com>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
