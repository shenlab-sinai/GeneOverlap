\name{testGeneOverlap}
\alias{testGeneOverlap}
\alias{testGeneOverlap,GeneOverlap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Test function for the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
    Perform Fisher's exact test based on the information supplied in the 
    GeneOverlap object, i.e. gene list A, B and genome size. This function 
    also calculates the Jaccard index. Will set the tested Boolean label after 
    done.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\S4method{testGeneOverlap}{GeneOverlap}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A GeneOverlap object.}
}
\value{
    A GeneOverlap object with valid p-value, odds ratio, Jaccard index and 
    contingency table. The tested Boolean label is set to true. Use show or 
    print to display a summary of the object. Use accessors to get information 
    of each slot.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlap-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
go.obj <- newGeneOverlap(hESC.ChIPSeq.list$H3K4me3, 
                         hESC.ChIPSeq.list$H3K27me3, 
                         genome.size=gs.RNASeq)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
print(go.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
