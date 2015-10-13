\name{newGOM}
\alias{newGOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructor for the GeneOverlapMatrix class
%%  ~~function to do ... ~~
}
\description{
Use this function to create objects of the GeneOverlapMatrix class.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
newGOM(gsetA, gsetB=list(), genome.size=NULL, 
    spec=c('mm9.gene', 'hg19.gene', 'rn4.gene'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsetA}{
  Gene set A as a named list with each element being a vector/factor of gene names.
}
  \item{gsetB}{
  Gene set B as a named list with each element being a vector/factor of gene names.
}
  \item{genome.size}{
  The number of genes in the genome as an integer.
}
  \item{spec}{
  A string description of the genome to use. There are a few presetted genome sizes to choose from for a user's convenience.
}
%%     ~~Describe \code{x} here~~
}
\details{
This will create a matrix so that each grid represents the comparison between the two gene lists from the two gene sets. Given two gene sets A and B, the matrix will represent all comparisons bewteen gene lists in A vs. gene lists in B. The set A will be shown as rows and the set B will be shown as columns. If only gene set A is given, the matrix will represent all comparisons between gene lists within the gene set and only the upper triangular matrix will be used.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A GeneOverlapMatrix object. Use accesssors to access to its internal structure and members. Use show or print to obtain summarized information. Use \code{drawHeatmap} to visualize it.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlapMatrix-class}}, \code{\link{GeneOverlap-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
gom.obj <- newGOM(hESC.ChIPSeq.list, hESC.RNASeq.list, gs.RNASeq)
gom.obj
drawHeatmap(gom.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
