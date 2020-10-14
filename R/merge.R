#' Merge a subset of the rows or cols
#'
#' In the context of microbiome analysis, it might be desirable to merge data
#' based on taxonomic levels. \code{mergeRows} implements the low level function
#' for more detail have a look at
#' \code{\link[=agglomerate-methods]{agglomerateByRank}}.
#'
#' @param x \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param f A factor for merging. Must be the same length as
#'   \code{nrow(x)/ncol(x)}. Row corresponding to the same level will be merged.
#'   If \code{length(levels(f)) == nrow(x)/ncol(x)}, \code{x} will be returned
#'   unchanged.
#'
#' @param archetype Of each level of \code{f}, which element should be regarded
#'   as the archetype and metadata in the columns or rows kept, while merging?
#'   This can be single interger value or an integer vector of the same length
#'   as \code{levels(f)}. (Default: \code{archetype = 1L}, which means the first
#'   element encountered will be kept)
#'
#' @param mergeTree \code{TRUE} or \code{FALSE}: should to
#'   \code{rowTree()} also be merged? (Default: \code{mergeTree = FALSE})
#'
#' @param ... optional arguments passed onto
#'   \code{\link[scater:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#'   except \code{subset_row} and \code{subset_col}
#'
#' @return an object with the same class \code{x} with the specified entries
#'   merged into one entry in all relevant components.
#'
#' @seealso \code{\link[=agglomerate-methods]{agglomerateByRank}}
#'
#' @name merge-methods
#'
#' @examples
#' data(esophagus)
#' esophagus
#' plot(rowTree(esophagus))
#' # get a factor for merging
#' f <- factor(regmatches(rownames(esophagus),
#'                        regexpr("^[0-9]*_[0-9]*",rownames(esophagus))))
#' merged <- mergeRows(esophagus,f)
#' plot(rowTree(merged))
#' #
#' data(GlobalPatterns)
#' GlobalPatterns
#' merged <- mergeCols(GlobalPatterns,colData(GlobalPatterns)$SampleType)
#' merged
NULL

#' @rdname merge-methods
#' @importFrom SEtup mergeRows
#' @importFrom TreeSummarizedExperiment rowTree convertNode rowLinks changeTree
#' @importFrom ape keep.tip
#' @export
setMethod("mergeRows", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, ...){
        # input check
        if(!.is_a_bool(mergeTree)){
          stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        x <- callNextMethod(x, f, archetype = 1L, ...)
        # optionally merge rowTree
        row_tree <- rowTree(x)
        if(!is.null(row_tree) && mergeTree){
          row_leaf <- convertNode(tree = row_tree, node = rowLinks(x)$nodeNum)
          row_tree <- ape::keep.tip(phy = row_tree, tip = row_leaf)
          x <- changeTree(x, rowTree = row_tree)
        }
        x
    }
)

#' @rdname merge-methods
#' @importFrom SEtup mergeCols
#' @importFrom TreeSummarizedExperiment colTree convertNode colLinks changeTree
#' @importFrom ape keep.tip
#' @export
setMethod("mergeCols", signature = c(x = "TreeSummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, ...){
        # input check
        if(!.is_a_bool(mergeTree)){
          stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        x <- callNextMethod(x, f, archetype = 1L, ...)
        # optionally merge colTree
        col_tree <- colTree(x)
        if(!is.null(col_tree) && mergeTree){
          col_leaf <- convertNode(tree = col_tree, node = colLinks(x)$nodeNum)
          col_tree <- ape::keep.tip(phy = col_tree, tip = col_leaf)
          x <- changeTree(x, colTree = col_tree)
        }
    }
)
