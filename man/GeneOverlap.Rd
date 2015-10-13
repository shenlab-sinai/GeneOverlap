\name{GeneOverlap}
\alias{GeneOverlap-class}
\alias{show,GeneOverlap-method}
\alias{print,GeneOverlap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Test overlap between two gene lists using Fisher's exact test.
}
\description{
    Given two gene lists, tests the significance of their overlap in comparison
    with a genomic background. The null hypothesis is that the odds ratio is no
    larger than 1. The alternative is that the odds ratio is larger than 1.0. 
    It returns the p-value, estimated odds ratio and intersection.
}
\usage{
\S4method{show}{GeneOverlap}(object)
\S4method{print}{GeneOverlap}(x, ...)
}
\arguments{
  \item{object}{A GeneOverlap object.}
  \item{x}{A GeneOverlap object.}
  \item{...}{They are not used.}
}
\details{
    The problem of gene overlap testing can be described by a hypergeometric
    distribution where one gene list A defines the number of white balls in the 
    urn and the other gene list B defines the number of white balls in the 
    draw. Assume the total number of genes is \var{n}, the number of genes in A 
    is \var{a} and the number of genes in B is \var{b}. If the intersection 
    between A and B is \var{t}, the probability density of seeing \var{t} can 
    be calculated as:

    \code{dhyper(t, a, n - a, b)}

    without loss of generality, we can assume \var{b} <= \var{a}. So the 
    largest possible value for \var{t} is \var{b}. Therefore, the p-value of
    seeing intersection \var{t} is:

    \code{sum(dhyper(t:b, a, n - a, b))}

    The Fisher's exact test forms this problem slightly different but its
    calculation is also based on the hypergeometric distribution. It starts by
    constructing a contingency table:

    \code{matrix(c(n - union(A,B), setdiff(A,B), 
        setdiff(B,A), intersect(A,B)), 
        nrow=2)}

    It therefore tests the independence between A and B and is conceptually 
    more straightforward. The GeneOverlap class is implemented using Fisher's 
    exact test. 

    It is better to illustrate a concept using some example. Let's assume we 
    have a genome of size 200 and two gene lists with 70 and 30 genes each. If 
    the intersection between the two is 10, the hypergeometric way to calculate 
    the p-value is:

    sum(dhyper(10:30, 70, 130, 30))

    which gives us p-value 0.6561562. If we use Fisher's exact test, we should 
    do:

    \code{fisher.test(matrix(c(110, 20, 60, 10), nrow=2), 
        alternative="greater")}

    which gives exactly the same p-value. In addition, the Fisher's test 
    function also provides an estimated odds ratio, confidence interval, etc.
    
    The Jaccard index is a measurement of similarity between two sets. It is 
    defined as the number of intersections over the number of unions.
%%  ~~ If necessary, more details than the description above ~~
}
\references{
\url{http://en.wikipedia.org/wiki/Fisher's_exact_test}

\url{http://en.wikipedia.org/wiki/Jaccard_index}
%% ~put references to the literature/web site here ~
}
\author{
Li Shen <\email{shenli.sam@gmail.com}>

Lab:\url{http://shenlab-sinai.github.io/shenlab-sinai/}
    
Personal:\url{http://www.linkedin.com/in/lshen/}
%%  ~~who you are~~
}
\note{
    Although Fisher's exact test is chosen for implementation, it should be 
    noted that the R implementation of Fisher's exact test is slower than using 
    \code{dhyper} directly. As an example, run:

    \code{system.time(sum(dhyper(10e3:30e3, 70e3, 130e3, 30e3)))}

    uses around 0.016s to finish. While run:

    \code{system.time(fisher.test(matrix(c(110e3, 20e3, 60e3, 10e3), nrow=2), 
        alternative="greater"))}

    uses around 0.072s. In practice, this time difference can often be ignored.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlapMatrix-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(GeneOverlap)
go.obj <- newGeneOverlap(hESC.ChIPSeq.list$H3K4me3, 
                         hESC.ChIPSeq.list$H3K9me3, 
                         gs.RNASeq)
go.obj <- testGeneOverlap(go.obj)
go.obj  # show.
print(go.obj)  # more details.
getContbl(go.obj)  # contingency table.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
