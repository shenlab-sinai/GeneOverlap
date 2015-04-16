#### GeneOverlap ####
setClass(
    "GeneOverlap", 
    representation(listA="character",
                   listB="character",
                   intersection="character",
                   union="character",
                   genome.size="numeric",
                   cont.tbl="matrix",
                   odds.ratio="numeric",
                   pval="numeric",
                   Jaccard="numeric",
                   is.tested="logical"),
    validity=function(object) {
        if(length(object@listA) > 0 && is.na(object@listA)) {
            stop("listA cannot be NA. Check your input.")
        }
        if(length(object@listB) > 0 && is.na(object@listB)) {
            stop("listB cannot be NA. Check your input.")
        }
        
        if(length(object@union) > object@genome.size) {
            stop("Union must NOT be larger than genome size.")
        }
    }
)

# Constructor
newGeneOverlap <- function(listA, listB, genome.size=NULL, 
                           spec=c('mm9.gene', 'hg19.gene', 'rn4.gene')) {
    listA <- unique(as.character(listA))
    listB <- unique(as.character(listB))
    listA <- listA[!is.na(listA)]
    listB <- listB[!is.na(listB)]
    
    # Setup genome size.
    if(is.null(genome.size)){
        spec <- match.arg(spec)
        genome.size <- switch(spec, 
                              mm9.gene=23000, 
                              hg19.gene=25000, 
                              rn4.gene=17000)
    }
    genome.size <- as.integer(genome.size)
    
    new("GeneOverlap", listA=listA, listB=listB, 
        intersection=intersect(listA, listB),
        union=union(listA, listB), 
        genome.size=genome.size, is.tested=F)
}


#### GeneOverlapMatrix ####
setClass(
    "GeneOverlapMatrix", 
    representation(gsetA="list",
                   gsetB="list",
                   self.compare="logical",
                   go.nested.list="list"),
    validity=function(object) {
        if(length(object@gsetB) == 0) {
            stopifnot(length(object@gsetA) > 1 && object@self.compare)
        } else {
            stopifnot(length(object@gsetA) > 0 && !object@self.compare)
        }
        
    }
)

# Constructor.
newGOM <- function(gsetA, gsetB=list(), genome.size=NULL, 
                   spec=c('mm9.gene', 'hg19.gene', 'rn4.gene')) {
    stopifnot(is.list(gsetA) && is.list(gsetB))
    # Construct GeneOverlap objects for all pairwise comparisons.
    if(length(gsetB) == 0) {
        stopifnot(length(gsetA) >= 2)
        self.compare <- T
        row.iter <- 1:(length(gsetA) - 1)
        col.iter <- 2:length(gsetA)
        go.nested.list <- 
            lapply(col.iter, function(ci) {
                this.col <- lapply(row.iter, function(ri) {
                    if(ri >= ci) {
                        go.obj <- newGeneOverlap(NULL, NULL)  # same list.
                        testGeneOverlap(go.obj)
                    } else {
                        go.obj <- newGeneOverlap(gsetA[[ri]], gsetA[[ci]], 
                                                 genome.size, spec)
                        testGeneOverlap(go.obj)
                    }
                })
                names(this.col) <- names(gsetA)[row.iter]
                this.col
            })
        names(go.nested.list) <- names(gsetA)[col.iter]
    } else {
        stopifnot(length(gsetA) >= 1 && length(gsetB) >= 1)
        self.compare <- F
        go.nested.list <- 
            lapply(gsetB, function(b) {
                this.col <- lapply(gsetA, function(a) {
                    go.obj <- newGeneOverlap(a, b, genome.size, spec)
                    testGeneOverlap(go.obj)
                })
                names(this.col) <- names(gsetA)
                this.col
            })
        names(go.nested.list) <- names(gsetB)
    }
    
    new("GeneOverlapMatrix", gsetA=gsetA, gsetB=gsetB, 
        self.compare=self.compare, go.nested.list=go.nested.list)
}

















