## Class constructor
.MicrobiomeExperiment <- setClass("MicrobiomeExperiment",
    contains="TreeSummarizedExperiment",
    representation(
        microbiomeData = "MicrobiomeFeatures"
    )
)

#' The MicrobiomeExperiment class
#'
#' SummarizedExperiment-like class for microbiome data.
#'
#' @param ... Arguments passed to \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' @param microbiomeData a \code{MicrobiomeFeatures} object or some object
#'   coercible to a \code{MicrobiomeFeatures} object. See
#'   \code{\link[=MicrobiomeFeatures-class]{MicrobiomeFeatures}} for more
#'   details.
#'
#'
#' @name MicrobiomeExperiment-class
#'
#' @include MicrobiomeFeatures-class.R
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @importClassesFrom metagenomeFeatures mgFeatures
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @examples
#' library(metagenomeFeatures)
#' data(mock_mgF)
#'
#' sampleNames <- letters[1:4]
#' pd <- DataFrame(a=letters[1:4], b=1:4)
#' numcounts <- nrow(mock_mgF) * 4
#' counts <- matrix(sample(1:1000, numcounts, replace=TRUE),
#'                  nr = nrow(mock_mgF), nc = 4)
#'
#' MicrobiomeExperiment(assays = SimpleList(counts = counts),
#'                      rowData = mock_mgF,
#'                      colData = pd )
NULL

################################################################################
# constructor

#' @rdname MicrobiomeExperiment-class
#' @export
MicrobiomeExperiment <- function(..., microbiomeData = list()) {
    tse <- TreeSummarizedExperiment(...)
    .tse_to_me(tse, microbiomeData)
}

.tse_to_me <- function(tse, microbiomeData = MicrobiomeFeatures()){
    out <- new("MicrobiomeExperiment", tse)
    microbiomeData(out) <- microbiomeData
    out
}


################################################################################
# coecrion

.from_SE_to_ME <- function(from){
    from <- as(from,"TreeSummarizedExperiment")
    .tse_to_me(from)
}

setAs("SummarizedExperiment","MicrobiomeExperiment",
      .from_SE_to_ME)


################################################################################
# accessors

#' Microbiome data methods
#'
#' Methods to get or set microbiome data on a
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}} object.
#'
#' @param x a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object
#'
#' @param value a a \code{\link[=MicrobiomeFeatures-class]{MicrobiomeFeatures}}
#'   object or an object coercible to one.
#'
#' @rdname MicrobiomeExperiment-class
#' @export
setGeneric("microbiomeData", signature = c("x"),
           function(x) standardGeneric("microbiomeData"))
#' @rdname MicrobiomeExperiment-class
#' @export
setMethod("microbiomeData", signature = c(x = "MicrobiomeExperiment"),
    function(x){
        x@microbiomeData
    }
)

#' @rdname MicrobiomeExperiment-class
#' @export
setGeneric("microbiomeData<-", signature = c("x"),
           function(x, value) standardGeneric("microbiomeData<-"))
#' @rdname MicrobiomeExperiment-class
#' @export
setReplaceMethod("microbiomeData", signature = c(x = "MicrobiomeExperiment"),
    function(x, value){
        if(!is(value,"MicrobiomeFeatures")){
          value <- as(value,"MicrobiomeFeatures")
        }
        .set_microbiomeData(x, value)
    }
)

.set_microbiomeData <- function(x, value){
    x@microbiomeData <- value
    x
}
