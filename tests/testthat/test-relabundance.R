context("relabundance")
test_that("relabundance", {
    mat <- matrix(1:60, nrow = 6)
    df <- DataFrame(n = c(1:6))
    expect_error(relAbundanceCounts(SummarizedExperiment(assays = list(mat = mat),
                                                         rowData = df)),
                 "'abund_values' must be a valid name of assays")
    se <- SummarizedExperiment(assays = list(counts = mat),
                               rowData = df)
    actual <- relAbundanceCounts(se)
    expect_named(assays(actual), c("counts", "relabundance"))
    expect_equal(assay(actual,"relabundance")[,1],
                 seq.int(1,6)/21)
})
