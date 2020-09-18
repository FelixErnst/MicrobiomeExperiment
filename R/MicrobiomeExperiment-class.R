## Class constructor
.MicrobiomeExperiment <- setClass("MicrobiomeExperiment",
    contains="TreeSummarizedExperiment",
    representation(
        rowData="MicrobiomeFeatures"
    )
)

#' The MicrobiomeExperiment class
#'
#' SummarizedExperiment-like class for microbiome data. rowData is
#' a MicrobiomeFeatures object LINK so it includes: a taxonomy table
#' (DataFrame), #' optional phylogentic tree (phylo object), and a sequence
#' database.
#'
#' It supports (most) of the interface to phyloseq objects
#'
#' @include MicrobiomeFeatures-class.R
#' @importClassesFrom metagenomeFeatures mgFeatures
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#'
#' library(metagenomeFeatures)
#' data(mock_mgF)
#'
#' sampleNames <- letters[1:4]
#' pd <- DataFrame(a=letters[1:4], b=1:4)
#' numcounts <- nrow(mock_mgF) * 4
#' counts <- matrix(sample(1:1000,numcounts,replace=TRUE), nr=nrow(mock_mgF), nc=4)
#'
#' MicrobiomeExperiment(assays=SimpleList(counts=counts),
#' rowData=mock_mgF,
#' colData=pd
#' )
#'
#' @name MicrobiomeExperiment-class
#' @export
MicrobiomeExperiment <-
    function(assays = SimpleList(),
             rowData = MicrobiomeFeatures(), ...) {

    if (is.data.frame(rowData) || is(rowData, "DataFrame")){
        rowData <- as(rowData, "MicrobiomeFeatures")
    }

    if (!is(assays, "SummarizedExperiment")) {
        se <- SummarizedExperiment(assays = assays, ...)
    } else {
        se <- assays
    }
    .MicrobiomeExperiment(se, rowData = rowData)
}
