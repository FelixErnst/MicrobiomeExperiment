context("merge")
test_that("merge", {
  # .check_f
  expect_error(Mia:::.check_f(),
               'argument "i" is missing')
  expect_error(Mia:::.check_f(6),
               'argument "f" is missing')
  expect_error(Mia:::.check_f(6,5),
               "'f' must have the same number of rows")
  f <- factor(c(rep("a",3),rep("b",3)))
  expect_null(Mia:::.check_f(6,f))
  # .check_archetype
  expect_error(Mia:::.check_archetype(),
               'argument "archetype" is missing')
  expect_error(Mia:::.check_archetype(f),
               'argument "archetype" is missing')
  expect_error(Mia:::.check_archetype(f),
               'argument "archetype" is missing')
  expect_null(Mia:::.check_archetype(f, 1))
  expect_null(Mia:::.check_archetype(f, c(1,2)))
  expect_error(Mia:::.check_archetype(f, c(1,2,3)),
               "length of 'archetype' must have the same length as levels")
  expect_error(Mia:::.check_archetype(f, c(5)),
               "'archetype' out of bounds for some levels of 'f'")
  # .norm_archetype
  expect_error(Mia:::.norm_archetype(),
               'argument "archetype" is missing')
  expect_error(Mia:::.norm_archetype(f),
               'argument "archetype" is missing')
  actual <- Mia:::.norm_archetype(f, c(1,2))
  expect_equal(actual, c(1,2))
  actual <- Mia:::.norm_archetype(f, c(1))
  expect_equal(actual, c(1,1))
  # .merge_assay
  expect_error(Mia:::.merge_assay(),
               'argument "assay" is missing')
  mat <- matrix(1:60, nrow = 6)
  expect_error(Mia:::.merge_assay(mat),
               'argument "f" is missing')
  actual_mat <- actual <- Mia:::.merge_assay(mat, f)
  expect_true(is.matrix(actual))
  expect_type(actual,"double")
  expect_equal(nrow(actual),2)
  expect_equal(ncol(actual),ncol(mat))
  FUN_check_mat_values <- function(col,mat,actual){
    expect_equal(sum(mat[1:3,col]),actual[1,col])
    expect_equal(sum(mat[4:6,col]),actual[2,col])
  }
  lapply(seq_len(ncol(actual)), FUN_check_mat_values, mat, actual)
  # .merge_row_or_col_data
  gr <- GRanges("chr1",rep("1-6",6))
  df <- DataFrame(n = c(1:6))
  mcols(gr) <- df
  grl <- splitAsList(gr,1:6)
  expect_error(Mia:::.merge_row_or_col_data(),
               'argument "f" is missing')
  expect_error(Mia:::.merge_row_or_col_data(df),
               'argument "f" is missing')
  actual <- Mia:::.merge_row_or_col_data(df, f)
  expect_s4_class(actual,"DataFrame")
  expect_equal(actual$n,c(1,4))
  actual <- Mia:::.merge_row_or_col_data(df, f, archetype = 2)
  expect_equal(actual$n,c(2,5))
  actual <- Mia:::.merge_row_or_col_data(df, f, archetype = c(3,1))
  expect_equal(actual$n,c(3,4))
  actual <- Mia:::.merge_row_ranges(gr, f)
  expect_s4_class(actual,"GRanges")
  expect_equal(mcols(actual)$n,c(1,4))
  actual <- Mia:::.merge_row_ranges(gr, f, archetype = 2)
  expect_equal(mcols(actual)$n,c(2,5))
  actual <- Mia:::.merge_row_ranges(gr, f, archetype = c(3,1))
  expect_equal(mcols(actual)$n,c(3,4))
  actual <- Mia:::.merge_row_ranges(grl, f)
  expect_s4_class(actual,"GRangesList")
  expect_equal(unlist(mcols(actual,level="within"))$n,c(1,4))
  actual <- Mia:::.merge_row_ranges(grl, f, archetype = 2)
  expect_equal(unlist(mcols(actual,level="within"))$n,c(2,5))
  actual <- Mia:::.merge_row_ranges(grl, f, archetype = c(3,1))
  expect_equal(unlist(mcols(actual,level="within"))$n,c(3,4))
  # .merge_rows
  expect_error(Mia:::.merge_rows(),
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
  expect_error(Mia:::.merge_rows(x),
               'argument "f" is missing')
  FUN_check_x <- function(x,actual_mat,archetype=1){
    actual <- Mia:::.merge_rows(x, f, archetype)
    expect_s4_class(actual,class(x))
    expect_equal(dim(actual),c(2,10))
    expect_equal(assays(actual)$mat,actual_mat)
  }
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_x,actual_mat)
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_x,actual_mat,archetype=2)
  #
  expect_equal(Mia:::.merge_rows(x, f), mergeRows(x, f))
  #
  f <- factor(c(rep("a",5),rep("b",5)))
  mat <- t(mat)
  actual_mat <- actual <- Mia:::.merge_assay(mat, f)
  FUN_check_mat_values <- function(col,mat,actual){
    expect_equal(sum(mat[1:5,col]),actual[1,col])
    expect_equal(sum(mat[6:10,col]),actual[2,col])
  }
  lapply(seq_len(ncol(actual)), FUN_check_mat_values, mat, actual)
  FUN_check_cols_x <- function(x,actual_mat,archetype=1){
    actual <- Mia:::.merge_cols(x, f, archetype)
    expect_s4_class(actual,class(x))
    expect_equal(dim(actual),c(6,2))
    expect_equal(assays(actual)$mat,actual_mat)
  }
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_cols_x,t(actual_mat))
  lapply(list(x,xr,xrl,xsce,xtse),FUN_check_cols_x,t(actual_mat),archetype=2)
  #
  expect_equal(Mia:::.merge_cols(x, f), mergeCols(x, f))
})
