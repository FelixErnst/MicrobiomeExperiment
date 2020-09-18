#' Loading a biom file
#'
#' For convenciance a few functions are available to convert data from a
#' \sQuote{biom} file or object into a
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'
#' @param file biom file location
#'
#' @return An object of class
#'   \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'
#' @name biom-load
#'
#' @examples
#' if(requireNameSpace("biomformat")) {
#'   # load from file
#'   rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
#'                                  package = "biomformat")
#'   loadFromBiom(rich_dense_file)
#'   # load from object
#'   rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
#'                                  package = "biomformat")
#'   x1 <- read_biom(rich_dense_file)
#'   as(x1, "MicrobiomeExperiment")
#' }
NULL

#' @rdname biom-load
#'
#' @importFrom biomformat read_biom
#'
#' @export
loadFromBiom <- function(file) {
    biom <- biomformat::read_biom(file)
    as(biom, "MicrobiomeExperiment")
}

#' @rdname biom-load
#'
#' @param obj object of type \code{\link[biomformat:read_biom]{biom}}
#'
#' @export
makeMicrobiomeExperimentFromBiom <- function(obj){
    if(!is(obj,"biom")){
        stop("'obj' must be a 'biom' object")
    }
    counts <- as(biom_data(obj), "matrix")
    sample_data <- sample_metadata(obj)
    feature_data <- observation_metadata(obj)

    MicrobiomeExperiment(assays=list(counts=counts), colData=sample_data,
                         rowData=MicrobiomeFeatures(taxa=feature_data))
}

setAs(from="biom", to="MicrobiomeExperiment", function(from) {
    makeMicrobiomeExperimentFromBiom(from)
})
