\name{GeneOverlap-package}
\alias{GeneOverlap-package}
\docType{package}
\title{
    Test and visualize overlaps between gene lists
}
\description{
    Given two sets of gene lists, this package calculates the overlaps between 
    all pairs of lists from the two sets. Fisher's exact test is used to 
    determine the p-value and odds ratio in comparison to a genomic background.
    Plotting functions are provided to visualize the results.
}
\details{
\tabular{ll}{
Package: \tab GeneOverlap\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-08-22\cr
License: \tab GPL-3\cr
}
    To use the package, construct one or two named lists each representing a 
    gene set. Each list should contain one or more vectors of gene names. Then 
    use GeneOverlapMatrix to perform pairwise comparisons. It will return an 
    object that can be used for visualization. The GeneOverlapMatrix calls 
    GeneOverlap internally to perform comparison between two gene lists.
}
\author{
Li Shen <\email{li.shen@mssm.edu}>

Mount Sinai profile:\url{http://www.mountsinai.org/profiles/li-shen}
    
Personal:\url{http://www.linkedin.com/in/lshen/}
}
\keyword{ htest }
\keyword{ hplot }
\keyword{ graphs }
