\name{getGenomeSize}
\alias{getGenomeSize}
\alias{getGenomeSize,GeneOverlap-method}
\alias{setGenomeSize<-}
\alias{setGenomeSize<-,GeneOverlap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Accessors for the "genome.size" slot of the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
    The genome.size slot contains the number of genes in the genome as an 
    integer.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\S4method{getGenomeSize}{GeneOverlap}(object)
\S4method{setGenomeSize}{GeneOverlap}(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A GeneOverlap object.}
  \item{value}{An integer representing genomic background.}
}
\details{
    After setGenomeSize function is called, the tested Boolean label will be 
    reset to false.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    An integer representing the genome size.
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
go.obj <- newGeneOverlap(hESC.ChIPSeq.list$H3K27me3, 
                         hESC.RNASeq.list$"Exp Medium", 
                         genome.size=gs.RNASeq)
getGenomeSize(go.obj)
v.gs <- c(12e3, 14e3, 16e3, 18e3, 20e3)
setNames(sapply(v.gs, function(g) {
    setGenomeSize(go.obj) <- g
    go.obj <- testGeneOverlap(go.obj)
    getPval(go.obj)
}), v.gs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
