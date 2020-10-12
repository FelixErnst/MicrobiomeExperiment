context("merge")
test_that("merge", {
    # .merge_rows
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    expect_error(MicrobiomeExperiment:::.merge_rows(),
                 'argument "f" is missing')
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    me <- MicrobiomeExperiment(assays = list(mat = mat),
                               rowRanges = unname(grl))
    expect_error(MicrobiomeExperiment:::.merge_rows(x),
                 'argument "f" is missing')
    FUN_check_x <- function(x,archetype=1){
      actual <- mergeRows(x, f, archetype)
      expect_s4_class(actual,class(x))
      expect_equal(dim(actual),c(2,10))
    }
    lapply(list(x,xr,xrl,xsce,xtse,me),FUN_check_x)
    lapply(list(x,xr,xrl,xsce,xtse,me),FUN_check_x,archetype=2)
})
