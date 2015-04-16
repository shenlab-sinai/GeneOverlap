#### Read-only ####
setMethod("getListA", "GeneOverlap", function(object) { object@listA } )
setMethod("getListB", "GeneOverlap", function(object) { object@listB } )
setMethod("getIntersection", "GeneOverlap", 
          function(object) { object@intersection } )
setMethod("getUnion", "GeneOverlap", function(object) { object@union } )
setMethod("getGenomeSize", "GeneOverlap", 
          function(object) { object@genome.size } )
setMethod("getTested", "GeneOverlap", function(object) { object@is.tested } )
setMethod(
    "getContbl", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@cont.tbl
        } else {
            warning("Test has not been performed yet.\n")
            matrix(nrow=0, ncol=0)
        }
    }
)
setMethod(
    "getPval", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@pval
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)
setMethod(
    "getOddsRatio", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@odds.ratio
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)
setMethod(
    "getJaccard", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@Jaccard
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)

#### Writable methods ####
setReplaceMethod(
    "setListA", "GeneOverlap",
    function(object, value) {
        object@listA <- as.character(value)
        object@intersection <- intersect(object@listA, object@listB)
        object@union <- union(object@listA, object@listB)
        object@is.tested <- F
        validObject(object)
        
        object
    }
)
setReplaceMethod(
    "setListB", "GeneOverlap",
    function(object, value) {
        object@listB <- as.character(value)
        object@intersection <- intersect(object@listA, object@listB)
        object@union <- union(object@listA, object@listB)
        object@is.tested <- F
        validObject(object)
        
        object
    }
)
setReplaceMethod(
    "setGenomeSize", "GeneOverlap",
    function(object, value) {
        object@genome.size <- value
        object@is.tested <- F
        validObject(object)
        
        object
    }
)



