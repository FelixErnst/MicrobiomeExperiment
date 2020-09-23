#' Getter / setter for relative abundance data
#'
#' \code{relabundance} is a getter/setter for relative abundance stored in the
#' assay slot \sQuote{relabundance} of a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name relabundance
#'
#' @importFrom SummarizedExperiment assays assays<-
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatterns <- relAbundanceCounts(GlobalPatterns)
#' relabundance(GlobalPatterns)
NULL

#' @rdname relabundance
setGeneric("relabundance", signature = c("x"),
           function(x, ...) standardGeneric("relabundance"))
#' @rdname relabundance
setGeneric("relabundance<-", signature = c("x"),
           function(x, value) standardGeneric("relabundance<-"))
#' @rdname relabundance
setGeneric("relAbundanceCounts", signature = c("x"),
           function(x, abund_values = "counts", name)
               standardGeneric("relAbundanceCounts"))

#' @rdname relabundance
#' @export
setMethod("relabundance",signature = c(x = "TreeSummarizedExperiment"),
  function(x){
    assays(x)[["relabundance"]]
  }
)

#' @rdname relabundance
#' @export
setReplaceMethod("relabundance", signature = c(x = "TreeSummarizedExperiment"),
  function(x, value){
   assays(x)[["relabundance"]] <- value
   x
  }
)

#' @rdname relabundance
#' @export
setMethod("relAbundanceCounts",signature = c(x = "TreeSummarizedExperiment"),
  function(x, abund_values = "counts", name){
    # input check
    if(missing(name)){
      name <- "relabundance"
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
