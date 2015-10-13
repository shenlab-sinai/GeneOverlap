\name{getReadonlyMatrix}
\alias{getGsetA}
\alias{getGsetA,GeneOverlapMatrix-method}
\alias{getGsetB}
\alias{getGsetB,GeneOverlapMatrix-method}
\alias{getSelfCompare}
\alias{getSelfCompare,GeneOverlapMatrix-method}
\alias{getMatrix}
\alias{getMatrix,GeneOverlapMatrix-method}
\alias{getNestedList}
\alias{getNestedList,GeneOverlapMatrix-method}
\alias{[}
\alias{[,GeneOverlapMatrix-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Read-only accessors for the various slots of the GeneOverlapMatrix class
%%  ~~function to do ... ~~
}
\description{
    The gsetA and gsetB slots contain the gene set A and B as named lists. The
    self.compare slot contains the Boolean label for whether self-comparison is
    performed. Use getMatrix to retrieve intersection, union, Jaccard index, 
    p-value and odds ratio as a matrix. Use getNestedList to retrieve the 
    intersection and union gene lists and contingency tables as a nested list 
    (outer list represents columns and inner list represents rows). Use 
    brackets [] to retrieve a particular GeneOverlap object within the 
    GeneOverlapMatrix object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\S4method{getGsetA}{GeneOverlapMatrix}(object)
\S4method{getGsetB}{GeneOverlapMatrix}(object)
\S4method{getSelfCompare}{GeneOverlapMatrix}(object)
\S4method{getMatrix}{GeneOverlapMatrix}(object, name=c("pval", 
    "odds.ratio", "intersection", "union", "Jaccard"))
\S4method{getNestedList}{GeneOverlapMatrix}(object, name=c(
    "intersection", "union", "cont.tbl"))
\S4method{[}{GeneOverlapMatrix}(x, i, j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A GeneOverlapMatrix object.}
  \item{x}{A GeneOverlapMatrix object.}
  \item{name}{A string description of the information to retrive. Use pval and
  odds.ratio to retrieve p-values and odds ratios. Use Jaccard to retrieve 
  Jaccard indices. In the context of matrix retrieval, intersection and union 
  will return the numbers of genes. While in the case of nested list retrieval,
  intersection and union will return the actual gene lists. Use cont.tbl to 
  retrieve the contingency tables.}
  \item{i, j}{Integer or character indices to retrieve GeneOverlap objects from
  the matrix.}
}
\details{
    When character indices are used, they should match the names of gsetA or 
    gsetB.
%%  ~~ If necessary, more details than the description above ~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlapMatrix-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
gom.obj <- newGOM(hESC.ChIPSeq.list, hESC.RNASeq.list, gs.RNASeq)
getMatrix(gom.obj, "odds.ratio")
inter.nl <- getNestedList(gom.obj, "intersection")
str(inter.nl)
go.k4.high <- gom.obj[1, 1]
go.k4.high
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
