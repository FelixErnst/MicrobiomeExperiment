#' Sample concentration methods
#'
#' Getter and setter for information on sample concentrations in a
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}} object.
#'
#' @param object A \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object.
#'
#' @param value A numeric vector of length equal to \code{ncol(object)},
#'   containing sample concentrations.
#'
#' @param ... Additional arguments, currently not used.
#'
#' @details
#' ToDO
#'
#' @return
#' For \code{sampleConcentrations} a numeric vector or \code{NULL}, if no
#' information are set.
#'
#' For \code{sampleConcentrations<-}, a modified \code{object} is returned with
#' sample concentrations stored as a column in \code{\link{colData}}.
#'
#' @seealso
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'
#' @author Felix G.M. Ernst
#'
#' @name sampleConcentrations
#'
#' @examples
#' data(soilrep)
#' soilrep <- as(soilrep,"MicrobiomeExperiment")
#' sampleConcentrations(soilrep) <- runif(ncol(soilrep))
#' sampleConcentrations(soilrep)
NULL

SAMPLE_CONCENTRATION_FIELD <- "sampleConcentration"

#' @rdname sampleConcentrations
#' @export
setGeneric("sampleConcentrations",
           function(object, ...)
               standardGeneric("sampleConcentrations"))
#' @rdname sampleConcentrations
#' @export
setGeneric("sampleConcentrations<-",
           function(object, ..., value)
               standardGeneric("sampleConcentrations<-"))

#' @rdname sampleConcentrations
#'
#' @importFrom SummarizedExperiment colData
#'
#' @export
setMethod("sampleConcentrations",
          signature = c(object = "MicrobiomeExperiment"),
    function(object) {
        colData(object)[[SAMPLE_CONCENTRATION_FIELD]]
    }
)

#' @rdname sampleConcentrations
#'
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
setReplaceMethod("sampleConcentrations",
                 signature = c(object = "MicrobiomeExperiment"),
    function(object, ..., value) {
        # input check
        if(!is.null(value) && !is.numeric(value)){
            stop("sampleConcentrations()<- accepts only numeric values.",
                 call. = FALSE)
        }
        #
        colData(object)[[SAMPLE_CONCENTRATION_FIELD]] <- value
        object
    }
)
