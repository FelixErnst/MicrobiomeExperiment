context("sampleConcentrations")

test_that("sampleConcentrations", {
    data(soilrep)
    soilrep <- as(soilrep,"MicrobiomeExperiment")
    expect_null(sampleConcentrations(soilrep))
    expect_error(sampleConcentrations(soilrep) <- rep("a",ncol(soilrep)),
                 "sampleConcentrations\\(\\)<- accepts only numeric values")
    data <- runif(ncol(soilrep))
    sampleConcentrations(soilrep) <- data
    expect_equal(data,sampleConcentrations(soilrep))
    sampleConcentrations(soilrep) <- NULL
    expect_null(sampleConcentrations(soilrep))
})
