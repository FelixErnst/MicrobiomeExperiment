#' @name merge-methods
#'
#' @title Merge a subset of the species in \code{x} into one species/taxa/OTU.
#'
#' @description
#' In the context of microbiome analysis, it might be desirable to merge data
#' based on taxonomic levels. \code{mergeRows} uses
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' and sub-classes, which contains certain field in its \code{rowData} (See
#' \code{\link[=taxonomy-methods]{taxonomyRanks}} for more details).
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
#'   if \code{useScater = TRUE}, except \code{subset_row}, \code{subset_col} and
#'   \code{average}
#'
#' @return an object with the same class \code{x} with the specified entries
#'   merged into one entry in all relevant components.
#'
#' @seealso \code{\link[=agglomerate-methods]{agglomerateByRank}}
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

.check_f <- function(i, f){
    if(i != length(f)){
        stop("'f' must have the same number of rows/length as 'x'", call. = FALSE)
    }
}

.check_archetype <- function(f, archetype){
    if(length(archetype) > 1L){
        if(length(levels(f)) != length(archetype)){
          stop("length of 'archetype' must have the same length as levels('f')",
               call. = FALSE)
        }
    }
    if(any(min(table(f)) < archetype)){
        stop("'archetype' out of bounds for some levels of 'f'. The maximum of",
             " 'archetype' is defined as min(table('f'))", call. = FALSE)
    }
}

.norm_archetype <- function(f, archetype){
    if(length(archetype) == 1L){
        archetype <- rep(archetype,length(levels(f)))
    }
    archetype
}

#' @importFrom S4Vectors splitAsList
.get_element_pos <- function(f, archetype = 1L){
    archetype <- .norm_archetype(f, archetype)
    f_pos <- seq.int(1L, length(f))
    f_pos_split <- S4Vectors::splitAsList(f_pos, f)
    f_pos <- unlist(f_pos_split[as.list(archetype)])
    f_pos
}

#' @importFrom S4Vectors SimpleList
#' @importFrom TreeSummarizedExperiment rowTree colTree convertNode rowLinks
#'   changeTree
#' @importFrom scater sumCountsAcrossFeatures
#' @importFrom ape keep.tip
.merge_rows <- function(x, f, archetype = 1L, mergeTree = FALSE, ...){
    # input check
    if(length(levels(f)) == nrow(x)){
      return(x)
    }
    .check_f(nrow(x), f)
    .check_archetype(f, archetype)
    if(!.is_a_bool(mergeTree)){
      stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
    }
    # merge assays
    assays <- assays(x)
    assays <- S4Vectors::SimpleList(lapply(assays, scater::sumCountsAcrossFeatures,
                                           ids = f, subset_row = NULL,
                                           subset_col = NULL, average = FALSE,
                                           ...))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[.get_element_pos(f, archetype = archetype),]
    assays(x, withDimnames = FALSE) <- assays
    if(is(x,"TreeSummarizedExperiment")){
      # optionally merge tree
      row_tree <- rowTree(x)
      if(!is.null(row_tree) && mergeTree){
        row_leaf <- convertNode(tree = row_tree, node = rowLinks(x)$nodeNum)
        row_tree <- ape::keep.tip(phy = row_tree, tip = row_leaf)
        x <- changeTree(x, rowTree = row_tree)
      }
    }
    x
}

#' @importFrom S4Vectors SimpleList
#' @importFrom TreeSummarizedExperiment rowTree colTree convertNode colLinks
#'   changeTree
#' @importFrom scater sumCountsAcrossCells
#' @importFrom ape keep.tip
.merge_cols <- function(x, f, archetype = 1L, mergeTree = FALSE, ...){
    # input check
    if(length(levels(f)) == ncol(x)){
      return(x)
    }
    .check_f(ncol(x), f)
    .check_archetype(f, archetype)
    if(!.is_a_bool(mergeTree)){
      stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
    }
    # merge col data
    col_data <- colData(x)[.get_element_pos(f, archetype = archetype),,drop=FALSE]
    # merge assays
    assays <- assays(x)
    assays <- S4Vectors::SimpleList(lapply(assays, scater::sumCountsAcrossCells,
                                           ids = f, subset_row = NULL,
                                           subset_col = NULL, average = FALSE,
                                           ...))
    names(assays) <- names(assays(x))
    # merge to result
    x <- x[,.get_element_pos(f, archetype = archetype)]
    assays(x, withDimnames = FALSE) <- assays
    if(is(x,"TreeSummarizedExperiment")){
      # optionally merge tree
      col_tree <- colTree(x)
      if(!is.null(col_tree) && mergeTree){
        col_leaf <- convertNode(tree = col_tree, node = colLinks(x)$nodeNum)
        col_tree <- ape::keep.tip(phy = col_tree, tip = col_leaf)
        x <- changeTree(x, colTree = col_tree)
      }
    }
    x
}

.merge_tree <- function(x, f, archetype = 1L){
    # input check
    if(length(levels(f)) == (x$Nnode + 1L)){
      return(x)
    }
    .check_f(x$Nnode + 1L, f)
    .check_archetype(f, archetype)
    #
    archetype <- .norm_archetype(f, archetype)
    f_index <- match(as.character(f),levels(f))
    names(f_index) <- seq.int(1L,length(f))
    f_index <- IRanges::CharacterList(split(f_index,f))
    f_index <- f_index[IRanges::IntegerList(as.list(archetype))]
    f_index <- as.integer(names(unlist(f_index,use.names = FALSE)))
    ape::keep.tip(x, f_index)
}

#' @rdname merge-methods
setGeneric("mergeRows",
           signature = "x",
           function(x, f, ...)
             standardGeneric("mergeRows"))

#' @rdname merge-methods
#' @export
setMethod("mergeRows", signature = c(x = "SummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, ...){
        .merge_rows(x, f, archetype = archetype, mergeTree = mergeTree,  ...)
    }
)

#' @rdname merge-methods
setGeneric("mergeCols",
           signature = "x",
           function(x, f, ...)
             standardGeneric("mergeCols"))

#' @rdname merge-methods
#' @export
setMethod("mergeCols", signature = c(x = "SummarizedExperiment"),
    function(x, f, archetype = 1L, mergeTree = FALSE, ...){
        .merge_cols(x, f, archetype = archetype, mergeTree = mergeTree, ...)
    }
)
