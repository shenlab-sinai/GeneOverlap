\name{newGeneOverlap}
\alias{newGeneOverlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructor for the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
Use this function to create objects of the GeneOverlap class.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
newGeneOverlap(listA, listB, genome.size = NULL, 
    spec = c("mm9.gene", "hg19.gene", "rn4.gene"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listA}{
  Gene list A. This should be a character vector or a factor.
}
  \item{listB}{
  Gene list B. This should be a character vector or a factor.
}
  \item{genome.size}{
  An integer represents the number of genes on the genome. If not specified, it will use the preset number based on "spec".
}
  \item{spec}{
  A character string of the genome name. Currently choose one of: mm9.gene, hg19.gene, rn4.gene. The gene numbers are based on protein coding genes.
}
}
\value{
A GeneOverlap object.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\note{
Both listA and listB will be converted to unique character vectors before testing, that means, the duplicated gene names are removed and therefore not counted.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(GeneOverlap)
go.obj <- newGeneOverlap(hESC.ChIPSeq.list$H3K4me3, 
                         hESC.ChIPSeq.list$H3K9me3, 
                         gs.RNASeq)
print(go.obj)  # not tested yet.
go.obj <- testGeneOverlap(go.obj)
print(go.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
