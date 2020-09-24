#' Getter / setter for relative abundance data
#'
#' \code{relabundance} is a getter/setter for relative abundance stored in the
#' assay slot \sQuote{relabundance} of a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @param x a  \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for calculating the relative abundance.
#' @param name A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to stor the calculated relative abundance.
#' @param value a \code{matrix} to store as the the \sQuote{relabundance} assay
#' @param ... optional arguments not used currently.
#'
#' @name relabundance
#'
#' @importFrom SummarizedExperiment assays assays<-
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' me <- as(GlobalPatterns,"MicrobiomeExperiment")
#' me <- relAbundanceCounts(me)
#' relabundance(me)
NULL

#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))
#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))
#' @rdname relabundance
setGeneric("relAbundanceCounts", signature = c("x"),
           function(x, abund_values = "counts", name = "relabundance")
               standardGeneric("relAbundanceCounts"))

#' @rdname relabundance
#' @export
setMethod("relabundance",signature = c(x = "MicrobiomeExperiment"),
  function(x){
    assays(x)[["relabundance"]]
  }
)

#' @rdname relabundance
#' @export
setReplaceMethod("relabundance", signature = c(x = "MicrobiomeExperiment"),
  function(x, value){
   assays(x)[["relabundance"]] <- value
   x
  }
)

#' @rdname relabundance
#' @export
setMethod("relAbundanceCounts",signature = c(x = "MicrobiomeExperiment"),
  function(x, abund_values = "counts", name = "relabundance"){
    # input check
    if(!.is_non_empty_string(name)){
        stop("'name' must be a single non-empty character value.",
             call. = FALSE)
    }
    .check_abund_values(abund_values, x)
    #
    data <- assays(x)[[abund_values]]
    value <- t(t(data)/colSums(data))
    dimnames(value) <- dimnames(data)
    assays(x)[[name]] <- value
    x
  }
)
