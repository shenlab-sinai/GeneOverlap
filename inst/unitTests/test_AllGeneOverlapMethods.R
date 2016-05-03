test_GeneOverlap <- function() {
    # Vanilla test.
    listA <- c("A", "B", "A")
    listB <- c("B", "C", "B")
    manual.contbl <- matrix(c(7, 1, 1, 1), nrow=2)
    fish.res <- fisher.test(manual.contbl, alternative="greater")
    go.obj <- newGeneOverlap(listA, listB, genome.size=10)
    checkEquals(getListA(go.obj), c("A", "B"))
    checkEquals(getListB(go.obj), c("B", "C"))
    checkEqualsNumeric(getGenomeSize(go.obj), 10)
    checkEqualsNumeric(length(getIntersection(go.obj)), 1)
    checkEqualsNumeric(length(getUnion(go.obj)), 3)
    checkEqualsNumeric(getPval(go.obj), NA)
    checkEqualsNumeric(getOddsRatio(go.obj), NA)
    go.obj <- testGeneOverlap(go.obj)
    checkEqualsNumeric(getPval(go.obj), fish.res$p.value)
    checkEqualsNumeric(getOddsRatio(go.obj), fish.res$estimate)
    
    # Gene lists contain invalid entries.
    listC <- c("A", NA, "B")
    go.obj <- newGeneOverlap(listC, listB, genome.size=10)
    checkEquals(getListA(go.obj), c("A", "B"))
    go.obj <- testGeneOverlap(go.obj)
    checkEqualsNumeric(getPval(go.obj), fish.res$p.value)
    checkEqualsNumeric(getOddsRatio(go.obj), fish.res$estimate)
    
    # Test NO overlap.
    listA <- c("A", "B")
    listB <- c("C", "D")
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    
    # Test absolute overlap (p-value=0).
    listA <- LETTERS
    listB <- LETTERS
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=100))
    checkEqualsNumeric(getPval(go.obj), 0)
    checkEqualsNumeric(getOddsRatio(go.obj), Inf)
    
    # Test empty gene list.
    go.obj <- testGeneOverlap(newGeneOverlap("A", NULL, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, "B", genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, NULL, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    
    # Unknown species.
    checkException(newGeneOverlap("A", "B", spec="speciesdoesnotexist"))

    # Genome smaller than gene lists combined.
    checkException(newGeneOverlap("A", "B", genome.size=1))
}

test_GeneOverlapMatrix <- function() {
    # Vanilla test.
    gv1 <- c("A", "B")
    gv2 <- c("B", "C")
    manual.contbl <- matrix(c(7, 1, 1, 1), nrow=2)
    fish.res <- fisher.test(manual.contbl, alternative="greater")
    gsetA <- list(A=gv1)
    gsetB <- list(B=gv2)
    gom.obj <- newGOM(gsetA, gsetB, genome.size=10)
    checkEquals(getGsetA(gom.obj), gsetA)
    checkEquals(getGsetB(gom.obj), gsetB)
    checkEquals(getSelfCompare(gom.obj), F)
    checkEqualsNumeric(getMatrix(gom.obj, "pval"), fish.res$p.value)
    checkEqualsNumeric(getMatrix(gom.obj, "odds.ratio"), fish.res$estimate)
    checkEqualsNumeric(getMatrix(gom.obj, "intersection"), 1)
    checkEqualsNumeric(getMatrix(gom.obj, "union"), 3)
    checkEquals(getNestedList(gom.obj, "intersection")[[1]][[1]], "B")
    checkEquals(getNestedList(gom.obj, "union")[[1]][[1]], c("A", "B", "C"))
    go.obj <- testGeneOverlap(newGeneOverlap(gv1, gv2, genome.size=10))
    checkEquals(gom.obj[1, 1], go.obj)
    checkEquals(gom.obj["A", "B"], go.obj)
    
    # Wrong inputs.
    checkException(newGOM(gsetA, matrix(c(1, 2, 3, 4), nrow=2)))
    
    # gsetA self-comparison not enough size.
    checkException(newGOM(gsetA))
    checkException(newGOM(list()))
    
    # gsetA cannot be empty.
    checkException(newGOM(list(), gsetB))
    
    # Accessing index out of boundary.
    checkException(gom.obj[2, 1])
    checkException(gom.obj[1, -2])
    checkException(gom.obj["C", 1])
    checkException(gom.obj[1, "C"])
    
}













