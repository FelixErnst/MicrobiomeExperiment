#' @include MicrobiomeExperiment.R
NULL

#' The MicrobiomeExperiment class
#'
#' SummarizedExperiment-like class for microbiome data.
#'
#' @param ... Arguments passed to \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' @param referenceSeq a \code{DNAStringSet} object or some object
#'   coercible to a \code{DNAStringSet} object. See
#'   \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} for more
#'   details.
#'
#'
#' @name MicrobiomeExperiment-class
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom Biostrings DNAStringSet
#'
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importClassesFrom Biostrings DNAStringSet
#'
#' @examples
#' data(taxa)
#'
#' sampleNames <- letters[1:4]
#' pd <- DataFrame(a=letters[1:4], b=1:4)
#' numcounts <- nrow(taxa) * 4
#' counts <- matrix(sample(1:1000, numcounts, replace=TRUE),
#'                  nr = nrow(taxa), nc = 4)
#'
#' MicrobiomeExperiment(assays = SimpleList(counts = counts),
#'                      rowData = taxa,
#'                      colData = pd)
NULL

setClassUnion("DNAStringSet_OR_NULL", c("DNAStringSet","NULL"))

#' @rdname MicrobiomeExperiment-class
#' @export
setClass("MicrobiomeExperiment",
         contains = "TreeSummarizedExperiment",
         slots = c(referenceSeq = "DNAStringSet_OR_NULL"),
         prototype = list(referenceSeq = NULL)
)


setMethod("vertical_slot_names", "MicrobiomeExperiment",
          function(x) c("referenceSeq", callNextMethod())
)

################################################################################
# validity


################################################################################
# constructor

#' @rdname MicrobiomeExperiment-class
#' @export
MicrobiomeExperiment <- function(..., referenceSeq = NULL) {
    tse <- TreeSummarizedExperiment(...)
    .tse_to_me(tse, referenceSeq)
}

.tse_to_me <- function(tse, referenceSeq = NULL){
    me <- new("MicrobiomeExperiment",
              tse,
              referenceSeq = referenceSeq)
    me
}

################################################################################
# coercion

setAs("TreeSummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(from)
})

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setAs("RangedSummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setAs("SummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setAs("SingleCellExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

################################################################################
# accessors

#' Microbiome data methods
#'
#' Methods to get or set reference sequence data on a
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}} object.
#'
#' @param x a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object
#'
#' @param value a a \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}
#'   object or an object coercible to one.
#'
#' @name referenceSeq
#'
#' @export
setGeneric("referenceSeq", signature = c("x"),
           function(x) standardGeneric("referenceSeq"))
#' @rdname referenceSeq
#' @export
setMethod("referenceSeq", signature = c(x = "MicrobiomeExperiment"),
    function(x){
        x@referenceSeq
    }
)

#' @rdname referenceSeq
#' @export
setGeneric("referenceSeq<-", signature = c("x"),
           function(x, value) standardGeneric("referenceSeq<-"))
#' @rdname referenceSeq
#' @export
setReplaceMethod("referenceSeq", signature = c(x = "MicrobiomeExperiment"),
    function(x, value){
        if(!is(value,"DNAStringSet")){
          value <- as(value,"DNAStringSet")
        }
        .set_referenceSeq(x, value)
    }
)

.set_referenceSeq <- function(x, value){
    x@referenceSeq <- value
    x
}

################################################################################
# subsetting


setMethod("[", signature = c("MicrobiomeExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              if (!missing(i)) {
                  x@referenceSeq <- referenceSeq(x)[i]
              }

              callNextMethod()
          }
)

setReplaceMethod("[", signature = c("MicrobiomeExperiment", "ANY", "ANY", "MicrobiomeExperiment"),
                 function(x, i, j, ..., value) {
                     if (missing(i) && missing(j)) {
                         return(value)
                     }

                     if (!missing(i)) {
                         tmp <- referenceSeq(x)
                         tmp[i] <- referenceSeq(value)
                         x@referenceSeq <- tmp
                     }
                     callNextMethod()
                 }
)


################################################################################
# show

setMethod("show", signature = c(object = "MicrobiomeExperiment"),
    function(object){
        callNextMethod()
        referenceSeq <- object@referenceSeq
        if(!is.null(referenceSeq)){
            msg <- sprintf(paste0("referenceSeq: a ", class(referenceSeq),
                                  " (%s sequences)\n"),
                           length(referenceSeq))
        } else {
            msg <- "referenceSeq: NULL\n"
        }
        cat(msg)
    }
)
