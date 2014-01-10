\name{getList}
\alias{getListA}
\alias{getListA,GeneOverlap-method}
\alias{setListA<-}
\alias{setListA<-,GeneOverlap-method}
\alias{getListB}
\alias{getListB,GeneOverlap-method}
\alias{setListB<-}
\alias{setListB<-,GeneOverlap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Accessors for the "listA" and "listB" slots of GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
    The listA and listB slots hold the gene lists A and B as character vectors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\S4method{getListA}{GeneOverlap}(object)
\S4method{setListA}{GeneOverlap}(object) <- value
\S4method{getListB}{GeneOverlap}(object)
\S4method{setListB}{GeneOverlap}(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A GeneOverlap object.}
  \item{value}{A character vector representing gene names.}
}
\details{
    After setListX function is called, the tested Boolean label will be reset 
    to false.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A character vector representing gene list A/B.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlap-class}}, \code{\link{newGeneOverlap}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
go.obj <- newGeneOverlap(hESC.ChIPSeq.list$H3K4me3, 
                         hESC.ChIPSeq.list$H3K27me3, 
                         genome.size=gs.RNASeq)
go.obj <- testGeneOverlap(go.obj)
head(getListB(go.obj))
getTested(go.obj)  # true.
setListB(go.obj) <- hESC.ChIPSeq.list$H3K9me3
getTested(go.obj) # false.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
