context("MicrobiomeExperiment subset methods")

test_that("[ subsetting works", {
    ## TODO: add toy example for testing here
    counts <- matrix(rbinom(20, 5, .5), ncol = 4)
    rowData <- DataFrame()
    me <- MicrobiomeExperiment(SimpleList(counts = counts), rowData = rowData)

    expect_identical(me[], me)
    expect_identical(me[,], me)
    expect_identical(me[TRUE, ], me)
    expect_identical(me[, TRUE], me)
    expect_identical(me[TRUE, TRUE], me)
})
