#### GeneOverlap ####
setGeneric("getListA", function(object) { standardGeneric("getListA")})
setGeneric("getListB", function(object) { standardGeneric("getListB")})
setGeneric("getIntersection", function(object) { 
    standardGeneric("getIntersection")})
setGeneric("getUnion", function(object) { 
    standardGeneric("getUnion")})
setGeneric("getGenomeSize", function(object) { 
    standardGeneric("getGenomeSize")})
setGeneric("getTested", function(object) { standardGeneric("getTested")})
setGeneric("getContbl", function(object) { standardGeneric("getContbl")})
setGeneric("getPval", function(object) { standardGeneric("getPval")})
setGeneric("getOddsRatio", function(object) { standardGeneric("getOddsRatio")})
setGeneric("getJaccard", function(object) { standardGeneric("getJaccard")})
setGeneric("setListA<-", function(object, value) { 
    standardGeneric("setListA<-") })
setGeneric("setListB<-", function(object, value) { 
    standardGeneric("setListB<-") })
setGeneric("setGenomeSize<-", function(object, value) { 
    standardGeneric("setGenomeSize<-") })
setGeneric("testGeneOverlap", function(object) { 
    standardGeneric("testGeneOverlap") })


#### GeneOverlapMatrix ####
setGeneric("getGsetA", function(object) { standardGeneric("getGsetA") } )
setGeneric("getGsetB", function(object) { standardGeneric("getGsetB") } )
setGeneric("getSelfCompare", function(object) { 
    standardGeneric("getSelfCompare")})
setGeneric("getMatrix", 
           function(object, name=c("pval", "odds.ratio", "intersection", 
                                   "union", "Jaccard")) { 
               standardGeneric("getMatrix")
           }
)
setGeneric("getNestedList", 
           function(object, name=c("intersection", "union", 
                                   "cont.tbl")) { 
               standardGeneric("getNestedList")
           }
)
setGeneric("drawHeatmap", 
           function(object, what=c("odds.ratio", "Jaccard"), log.scale=F, 
                    adj.p=F, cutoff=.05, ncolused=9, 
                    grid.col=c("Greens", "Blues", "Greys", 
                               "Oranges", "Purples", "Reds"),
                    note.col="red") { 
               standardGeneric("drawHeatmap") 
           }
)














