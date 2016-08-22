\name{GeneOverlapMatrix}
\alias{GeneOverlapMatrix-class}
\alias{show,GeneOverlapMatrix-method}
\alias{print,GeneOverlapMatrix-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Matrix representation of the pairwise overlaps between two gene sets
}
\description{
    Given one or two gene sets each contains one or more gene lists, create a 
    matrix to represent all pairwise comparisons between the gene lists. This 
    class also provides functions to visualize the matrix.
}
\usage{
\S4method{show}{GeneOverlapMatrix}(object)
\S4method{print}{GeneOverlapMatrix}(x, ...)
}
\arguments{
  \item{object}{A GeneOverlapMatrix object.}
  \item{x}{A GeneOverlapMatrix object.}
  \item{...}{They are not used.}
}
\details{
    The problem is stated as the representation of all pairwise comparisons 
    between two gene sets each contains one or more gene lists. This is 
    represented as a matrix where the rows correspond to one gene set and the 
    columns correspond to the other. Each grid represents the overlap between 
    two gene lists of the corresponding row and column. This class calls the
    GeneOverlap constructor to create objects that represent the overlapping
    information. When there is only one gene set, the matrix represents the 
    self-comparison within the gene set and only the upper triangular matrix is 
    used.

    The significance of gene overlap is characterized by two pieces of 
    information: odds ratio and p-value. This class provides functions to 
    visualize these information as a heatmap. The color gradients of each 
    grid represents the odds ratio while the texts superimposed on the grids 
    state the p-values. It is also possible to let the color gradients 
    represent Jaccard index - a measurement of similarity between two sets.
%%  ~~ If necessary, more details than the description above ~~
}
\author{
Li Shen <\email{li.shen@mssm.edu}>

Mount Sinai profile:\url{http://www.mountsinai.org/profiles/li-shen}
    
Personal:\url{http://www.linkedin.com/in/lshen/}
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlap-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
gom.obj <- newGOM(hESC.ChIPSeq.list, hESC.RNASeq.list, gs.RNASeq)
gom.obj
print(gom.obj)
drawHeatmap(gom.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
