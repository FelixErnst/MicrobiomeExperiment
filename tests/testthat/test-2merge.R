context("merge")
test_that("merge", {
    f <- factor(c(rep("a",3),rep("b",3)))
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    me <- MicrobiomeExperiment(assays = list(mat = mat),
                               rowRanges = unname(grl))
    FUN_check_x <- function(x,archetype=1){
      actual <- mergeRows(x, f, archetype)
      expect_s4_class(actual,class(x))
      expect_equal(dim(actual),c(2,10))
    }
    lapply(list(xtse,me),FUN_check_x)
    lapply(list(xtse,me),FUN_check_x,archetype=2)
})
