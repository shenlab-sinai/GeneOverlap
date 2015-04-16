#### show and print functions ####
setMethod("show", "GeneOverlapMatrix",
          function(object) {
              gom.dim <- dim(getMatrix(object, "pval"))
              cat(sprintf("A <%d x %d> GeneOverlapMatrix object\n", 
                          gom.dim[1], gom.dim[2]))
              gsetA <- getGsetA(object)
              gsetB <- getGsetB(object)
              cat("Geneset A sizes:\n")
              print(sapply(gsetA, length))
              if(getSelfCompare(object)) {
                  cat("Matrix is based on self-comparison of geneset A.\n")
              } else {
                  cat("Geneset B sizes:\n")
                  print(sapply(gsetB, length))
              }
          }
)

setMethod("print", "GeneOverlapMatrix",
          function(x, ...) {
              cat("A GeneOverlapMatrix object:\n")
              
              int.mat <- getMatrix(x, "intersection")
              cat("###### Intersection ######\n")
              print(int.mat)
              
              pval.mat <- getMatrix(x, "pval")
              cat("###### P-value ######\n")
              print(pval.mat)
              
              or.mat <- getMatrix(x, "odds.ratio")
              cat("###### Odds Ratio ######\n")
              print(or.mat)
              
              ja.mat <- getMatrix(x, "Jaccard")
              cat("###### Jaccard Index ######\n")
              print(ja.mat)
          }
)

#### Draw function ####
setMethod(
    "drawHeatmap", "GeneOverlapMatrix",
    function(object, what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, 
             cutoff=.05, ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                                "Oranges", "Purples", "Reds"),
             note.col="red") {
        
        # Arguments setup.
        stopifnot(cutoff > 0 && cutoff <= 1)
        what <- match.arg(what)
        grid.col <- match.arg(grid.col)
        
        # Matrix values.
        pv.mat <- getMatrix(object, "pval")
        plot.mat <- switch(what, 
                           odds.ratio=getMatrix(object, "odds.ratio"),
                           Jaccard=getMatrix(object, "Jaccard")
        )
        if(what == "odds.ratio" && log.scale) {
            plot.mat <- log2(plot.mat)
        }
        
        # Adjust p-values if needed.
        pv.mask <- NULL
        if(object@self.compare) {
            pv.mask <- sapply(1:ncol(pv.mat), function(j) {
                c(rep(T, j), rep(F, nrow(pv.mat) - j))
            })
        }
        if(adj.p) {
            if(object@self.compare) {
                pv.mat[pv.mask] <- p.adjust(pv.mat[pv.mask], method='BH')
            } else {
                pv.mat <- matrix(p.adjust(pv.mat, method='BH'), 
                                 nrow=nrow(pv.mat))
            }
        }
        
        # Marker value of insignificant events.
        insig.val <- 1
        if(what == "odds.ratio" && log.scale || what == "Jaccard") {
            insig.val <- 0
        }
        
        # Use p-value cutoff to mask insignificant cells.
        plot.mat[ pv.mat >= cutoff ] <- insig.val
        
        # Cell notes.
        note.mat <- format(pv.mat, digits=1)
        note.mat[pv.mat < .01] <- format(pv.mat, digits=1, 
                                         scientific=T)[pv.mat < .01]
        note.mat[plot.mat == insig.val] <- "N.S."
        if(object@self.compare) { note.mat[ !pv.mask ] <- "--" }
        
        # Configure heatmap graphic properties.
        row_sep <- 1:(nrow(plot.mat) - 1)
        col_sep <- 1:(ncol(plot.mat) - 1)
        longedge <- max(nrow(plot.mat), ncol(plot.mat))
        row_cexrc <- 0.4 + 1/log10(longedge + 2)
        col_cexrc <- row_cexrc
        key_size <- 0.2 + 1 / log10(longedge + 4)
        margins_use <- c(max(nchar(colnames(plot.mat))) * 0.8 + 5, 
                         max(nchar(rownames(plot.mat))) * 0.8 + 5)
        main.txt <- switch(what, 
                           odds.ratio=ifelse(log.scale, "log2(Odds Ratio)", 
                                             "Odds Ratio"), 
                           Jaccard="Jaccard Index")
        footnote <- "N.S.: Not Significant; --: Ignored"
        # sidenote <- sprintf("Log Scale=%s", log.scale)
        
        # Draw the heatmap!
        heatmap.2(plot.mat, cellnote=note.mat, 
                  main=main.txt, xlab=footnote, # ylab=sidenote,
                  col=brewer.pal(ncolused, grid.col), notecol=note.col, 
                  margins=margins_use, colsep=col_sep, rowsep=row_sep, 
                  key=T, keysize=key_size,
                  cexRow=row_cexrc, cexCol=col_cexrc, 
                  scale='none', Colv=NA, Rowv=NA, trace='none', 
                  dendrogram='none', density.info='none', 
                  sepcolor='white', sepwidth=c(0.002,0.002),
                  notecex=1.6)
    }
)










