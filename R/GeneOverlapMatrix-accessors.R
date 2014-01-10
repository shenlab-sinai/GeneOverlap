setMethod("getGsetA", "GeneOverlapMatrix", function(object) { object@gsetA } )
setMethod("getGsetB", "GeneOverlapMatrix", function(object) { object@gsetB } )
setMethod("getSelfCompare", "GeneOverlapMatrix",
          function(object) {
              object@self.compare
          }
)
setMethod(
    "getMatrix", "GeneOverlapMatrix",
    function(object, name=c("pval", "odds.ratio", "intersection", "union", 
                            "Jaccard")) {
        name <- match.arg(name)
        sapply(object@go.nested.list, function(ci) {
            sapply(ci, function(ri) {
                switch(name, 
                       pval=getPval(ri), 
                       odds.ratio=getOddsRatio(ri),
                       intersection=length(getIntersection(ri)),
                       union=length(getUnion(ri)),
                       Jaccard=getJaccard(ri)
                       )
            })
        })
    }
)
setMethod(
    "getNestedList", "GeneOverlapMatrix",
    function(object, name=c("intersection", "union", "cont.tbl")) {
        name <- match.arg(name)
        lapply(object@go.nested.list, function(ci) {
            lapply(ci, function(ri) {
                switch(name, 
                       intersection=getIntersection(ri),
                       union=getUnion(ri),
                       cont.tbl=getContbl(ri))
            })
        })
    }
)
setMethod(
    "[", "GeneOverlapMatrix",
    function(x, i, j) {
        stopifnot(is.numeric(i) || is.character(i))
        stopifnot(is.numeric(j) || is.character(j))
        if(is.numeric(j)) {
            j <- as.integer(j)
            stopifnot(abs(j) <= length(x@go.nested.list))
        } else {
            stopifnot(j %in% names(x@go.nested.list))
        }
        gom.col <- x@go.nested.list[[j]]
        if(is.numeric(i)) {
            i <- as.integer(i)
            stopifnot(abs(i) <= length(gom.col))
        } else {
            stopifnot(i %in% names(gom.col))
        }
        gom.col[[i]]
    }
)

















