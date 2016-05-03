#### Display methods ####
setMethod(
    "show", "GeneOverlap",
    function(object) {
        cat("GeneOverlap object:\n")
        cat(sprintf("listA size=%d\n", length(object@listA)))
        cat(sprintf("listB size=%d\n", length(object@listB)))
        cat(sprintf("Intersection size=%d\n", length(object@intersection)))
        if(object@is.tested) {
            cat(sprintf("Overlapping p-value=%s\n", 
                        ifelse(object@pval < .01, 
                               format(object@pval, scientific=T, digits=2),
                               format(object@pval, digits=2)
                        )))
            cat(sprintf("Jaccard Index=%.1f\n", object@Jaccard))
        } else {
            cat("Overlap testing has not been performed yet.\n")
        }
    }
)

setMethod(
    "print", "GeneOverlap",
    function(x, ...) {
        cat("Detailed information about this GeneOverlap object:\n")
        cat(sprintf("listA size=%d, e.g. %s\n", 
                    length(x@listA), 
                    paste(head(x@listA, n=3), collapse=" ")))
        cat(sprintf("listB size=%d, e.g. %s\n", 
                    length(x@listB), 
                    paste(head(x@listB, n=3), collapse=" ")))
        cat(sprintf("Intersection size=%d, e.g. %s\n", 
                    length(x@intersection),
                    paste(head(x@intersection, n=3), collapse=" ")))
        cat(sprintf("Union size=%d, e.g. %s\n", 
                    length(x@union),
                    paste(head(x@union, n=3), collapse=" ")))
        cat(sprintf("Genome size=%d\n", x@genome.size))
        if(x@is.tested) {
            cat("# Contingency Table:\n")
            print(x@cont.tbl)
            cat(sprintf("Overlapping p-value=%s\n", 
                        ifelse(x@pval < .01, 
                               format(x@pval, scientific=T, digits=2),
                               format(x@pval, digits=2)
                        )))
            cat(sprintf("Odds ratio=%.1f\n", x@odds.ratio))
            cat("Overlap tested using Fisher's exact test (alternative=greater)\n")
            cat(sprintf("Jaccard Index=%.1f\n", x@Jaccard))
        } else {
            cat("Overlap has not been tested yet. Use testGeneOverlap method.\n")
        }
    }
)

#### Test method ####
setMethod(
    "testGeneOverlap", "GeneOverlap",
    function(object) {
        # Configure contingency table.
        sizeA <- length(object@listA)
        sizeB <- length(object@listB)
        object@cont.tbl <- matrix(c(object@genome.size - length(object@union), 
                                    sizeB - length(object@intersection), 
                                    sizeA - length(object@intersection), 
                                    length(object@intersection)), 
                                  ncol=2)
        rownames(object@cont.tbl) <- c('notB', 'inB')
        colnames(object@cont.tbl) <- c('notA', 'inA')
        
        # Perform Fisher's exact test.
        res.fisher <- try(fisher.test(object@cont.tbl, alternative='greater'), 
                          silent=TRUE)
        if(is.list(res.fisher)) {
            object@odds.ratio <- setNames(res.fisher$estimate, NULL)
            object@pval <- res.fisher$p.value
        } else {
            object@odds.ratio <- .0
            object@pval <- 1.
        }
        
        # Calculate Jaccard index.
        object@Jaccard <- ifelse(length(object@union) == 0, 0, 
                                 length(object@intersection) / 
                                     length(object@union)
        )
        
        object@is.tested <- T
        object
    }
)
