context("MicrobiomeExperiment-class")

test_that("constructors work for basic cases", {
    expect_true(validObject(MicrobiomeExperiment()))
    expect_true(validObject(MicrobiomeExperiment(SimpleList())))

})
