context("MicrobiomeExperiment-class")

test_that("MicrobiomeExperiment-class", {
    expect_true(validObject(MicrobiomeExperiment()))
    expect_true(validObject(MicrobiomeExperiment(SimpleList())))
    data(taxa)
    sampleNames <- letters[1:4]
    pd <- DataFrame(a=letters[1:4], b=1:4)
    numcounts <- nrow(taxa) * 4
    counts <- matrix(sample(1:1000, numcounts, replace=TRUE),
                     nr = nrow(taxa), nc = 4)
    refSeq <- DNAStringSetList(one = DNAStringSet(c("A","A","A","A","A")),
                               two = DNAStringSet(c("A","A","A","A","A")))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = pd,
                                      referenceSeq = refSeq[list(1:5,1:4)]))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = pd,
                                      referenceSeq = refSeq[list(1:4,1:4)]))
    expect_error(MicrobiomeExperiment(assays = SimpleList(counts = counts),
                                      rowData = taxa,
                                      colData = pd,
                                      referenceSeq = refSeq[[1L]][1:4]))
    me <- MicrobiomeExperiment(assays = SimpleList(counts = counts),
                               rowData = taxa,
                               colData = pd,
                               referenceSeq = refSeq)
    expect_error(referenceSeq(me) <- refSeq[list(1:5,1:4)])
    expect_error(referenceSeq(me) <- refSeq[list(1:4,1:4)])
    expect_error(referenceSeq(me) <- refSeq[[1L]][1:4])
    #
    data(soilrep)
    # expect_output(as(as(soilrep,"SummarizedExperiment"),
    #                  "MicrobiomeExperiment"))
    expect_s4_class(as(as(soilrep,"SingleCellExperiment"),
                       "MicrobiomeExperiment"),
                    "MicrobiomeExperiment")
    expect_s4_class(as(soilrep, "MicrobiomeExperiment"),
                    "MicrobiomeExperiment")
    referenceSeq(me) <- refSeq
    expect_s4_class(referenceSeq(me),"DNAStringSetList")
    referenceSeq(me) <- as.list(refSeq)
    expect_s4_class(referenceSeq(me),"DNAStringSetList")
    referenceSeq(me) <- refSeq[[1L]]
    expect_s4_class(referenceSeq(me),"DNAStringSet")
    referenceSeq(me) <- as.character(refSeq[[1L]])
    expect_s4_class(referenceSeq(me),"DNAStringSet")
})
