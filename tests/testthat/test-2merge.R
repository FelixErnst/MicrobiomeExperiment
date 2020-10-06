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
  # .merge_assay
  expect_error(MicrobiomeExperiment:::.merge_assay(),
               'argument "assay" is missing')
  mat <- matrix(1:60, nrow = 6)
  expect_error(MicrobiomeExperiment:::.merge_assay(mat),
               'argument "f" is missing')
  actual_mat <- actual <- MicrobiomeExperiment:::.merge_assay(mat, f)
  expect_true(is.matrix(actual))
  expect_type(actual,"double")
  expect_equal(nrow(actual),2)
  expect_equal(ncol(actual),ncol(mat))
  FUN_check_mat_values <- function(col,mat,actual){
    expect_equivalent(sum(mat[1:3,col]), actual[1,col])
    expect_equivalent(sum(mat[4:6,col]), actual[2,col])
  }
  lapply(seq_len(ncol(actual)), FUN_check_mat_values, mat, actual)
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
  expect_error(MicrobiomeExperiment:::.merge_rows(x),
               'argument "f" is missing')
  FUN_check_x <- function(x,actual_mat,archetype=1){
    actual <- MicrobiomeExperiment:::.merge_rows(x, f, archetype)
    expect_s4_class(actual,class(x))
    expect_equal(dim(actual),c(2,10))
    expect_equal(assays(actual)$mat,actual_mat)
  }
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_x,actual_mat)
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_x,actual_mat,archetype=2)
  #
  expect_equal(MicrobiomeExperiment:::.merge_rows(x, f), mergeRows(x, f))
  # useScater = TRUE
  expect_equal(mergeRows(x, f), mergeRows(x, f, useScater = TRUE))
  #
  f <- factor(c(rep("a",5),rep("b",5)))
  mat <- t(mat)
  actual_mat <- actual <- MicrobiomeExperiment:::.merge_assay(mat, f)
  FUN_check_mat_values <- function(col,mat,actual){
      expect_equivalent(sum(mat[1:5,col]),actual[1,col])
      expect_equivalent(sum(mat[6:10,col]),actual[2,col])
  }
  lapply(seq_len(ncol(actual)), FUN_check_mat_values, mat, actual)
  FUN_check_cols_x <- function(x,actual_mat,archetype=1){
    actual <- MicrobiomeExperiment:::.merge_cols(x, f, archetype)
    expect_s4_class(actual,class(x))
    expect_equal(dim(actual),c(6,2))
    expect_equal(assays(actual)$mat,actual_mat)
  }
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_cols_x,t(actual_mat))
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_cols_x,t(actual_mat),archetype=2)
  #
  expect_equal(MicrobiomeExperiment:::.merge_cols(x, f), mergeCols(x, f))
  # useScater = TRUE
  expect_equal(mergeCols(x, f), mergeCols(x, f, useScater = TRUE))
})
