context("merge")
test_that("merge", {
    # .check_f
    expect_error(MicrobiomeExperiment:::.check_f(),
                 'argument "i" is missing')
    expect_error(MicrobiomeExperiment:::.check_f(6),
                 'argument "f" is missing')
    expect_error(MicrobiomeExperiment:::.check_f(6,5),
                 "'f' must have the same number of rows")
    f <- factor(c(rep("a",3),rep("b",3)))
    expect_null(MicrobiomeExperiment:::.check_f(6,f))
    # .check_archetype
    expect_error(MicrobiomeExperiment:::.check_archetype(),
                 'argument "archetype" is missing')
    expect_error(MicrobiomeExperiment:::.check_archetype(f),
                 'argument "archetype" is missing')
    expect_error(MicrobiomeExperiment:::.check_archetype(f),
                 'argument "archetype" is missing')
    expect_null(MicrobiomeExperiment:::.check_archetype(f, 1))
    expect_null(MicrobiomeExperiment:::.check_archetype(f, c(1,2)))
    expect_error(MicrobiomeExperiment:::.check_archetype(f, c(1,2,3)),
                 "length of 'archetype' must have the same length as levels")
    expect_error(MicrobiomeExperiment:::.check_archetype(f, c(5)),
                 "'archetype' out of bounds for some levels of 'f'")
    # .norm_archetype
    expect_error(MicrobiomeExperiment:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(MicrobiomeExperiment:::.norm_archetype(f),
                 'argument "archetype" is missing')
    actual <- MicrobiomeExperiment:::.norm_archetype(f, c(1,2))
    expect_equal(actual, c(1,2))
    actual <- MicrobiomeExperiment:::.norm_archetype(f, c(1))
    expect_equal(actual, c(1,1))
    # .get_element_pos
    expect_error(MicrobiomeExperiment:::.get_element_pos(),
                 'argument "f" is missing')
    actual <- MicrobiomeExperiment:::.get_element_pos(f)
    expect_equal(actual,c(a = 1, b = 4))
    actual <- MicrobiomeExperiment:::.get_element_pos(f, archetype = 2)
    expect_equal(actual,c(a = 2, b = 5))
    actual <- MicrobiomeExperiment:::.get_element_pos(f, archetype = c(2,1))
    expect_equal(actual,c(a = 2, b = 4))
    # .merge_rows
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    expect_error(MicrobiomeExperiment:::.merge_rows(),
                 'argument "f" is missing')
    x <- SummarizedExperiment(assays = list(mat = mat))
    xr <- SummarizedExperiment(assays = list(mat = mat),
                               rowRanges = gr)
    xrl <- SummarizedExperiment(assays = list(mat = mat),
                                rowRanges = unname(grl))
    xsce <- SingleCellExperiment(assays = list(mat = mat),
                                 rowRanges = unname(grl))
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
