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
#' @param ... optional arguments not used.
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

.merge_assay <- function(assay, f, type = "row", FUN = colSums){
  if(type == "col"){
    assay <- t(assay)
  }
  assay_split <- split(as.data.frame(assay), f)
  long_assay_split <- vapply(assay_split,nrow,numeric(1)) > 1L
  assay_split[long_assay_split] <- lapply(assay_split[long_assay_split],
                                          FUN)
  ans <- as.matrix(do.call(rbind,assay_split))
  if(type == "col"){
    ans <- t(ans)
  }
  ans
}

#' @importFrom IRanges SplitDataFrameList IntegerList
#' @importFrom S4Vectors splitAsList
.merge_row_or_col_data <- function(x, f, archetype = 1L){
  # input check
  if(length(levels(f)) == nrow(x)){
    return(x)
  }
  .check_f(nrow(x), f)
  .check_archetype(f, archetype)
  #
  archetype <- .norm_archetype(f, archetype)
  x_split <- S4Vectors::splitAsList(x,f)
  x_split <- IRanges::SplitDataFrameList(x_split)
  x_split <- x_split[IRanges::IntegerList(as.list(archetype))]
  unlist(x_split,use.names = FALSE)
}

#' @importFrom IRanges SplitDataFrameList IntegerList
#' @importFrom S4Vectors splitAsList
.merge_row_ranges <- function(x, f, archetype = 1L){
  # input check
  if(length(levels(f)) == length(x)){
    return(x)
  }
  .check_f(length(x), f)
  .check_archetype(f, archetype)
  #
  archetype <- .norm_archetype(f, archetype)
  x_split <- S4Vectors::splitAsList(x,f)
  x_split <- x_split[IRanges::IntegerList(as.list(archetype))]
  unlist(x_split,use.names = FALSE)
}

#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment colData rowRanges rowData rowData<-
#' @importFrom SingleCellExperiment int_colData int_colData<-
#'   int_elementMetadata int_elementMetadata<-
#'   int_metadata int_metadata<-
#' @importFrom TreeSummarizedExperiment rowTree colTree transNode rowLinks
#'   changeTree
#' @importFrom ape keep.tip
.merge_rows <- function(x, f, archetype = 1L, mergeTree = FALSE){
  # input check
  if(length(levels(f)) == nrow(x)){
    return(x)
  }
  .check_f(nrow(x), f)
  .check_archetype(f, archetype)
  if(!.is_a_bool(mergeTree)){
    stop("'mergeTree' must be TRUE or FALSE.", call. = FALSE)
  }
  # merge row data
  row_data <- .merge_row_or_col_data(rowData(x), f, archetype)
  # merge assays
  assays <- assays(x)
  assays <- S4Vectors::SimpleList(lapply(assays, .merge_assay, f, "row"))
  names(assays) <- names(assays(x))
  # merge to result
  class_name <- class(x)
  class_FUN <- class_name
  if(grepl("^Ranged",class_FUN)){
    class_FUN <- gsub("^Ranged","",class_FUN)
  }
  args <- list(assays = assays,
               colData = colData(x),
               metadata = metadata(x))
  ####################################
  # class specific pre-processing
  if(is(x,"RangedSummarizedExperiment")){
    # merge row ranges
    row_ranges <- .merge_row_ranges(rowRanges(x), f, archetype)
    args <- c(args,
              list(rowRanges = row_ranges))
  } else {
    args <- c(args,
              list(rowData = row_data))
  }
  if(is(x,"TreeSummarizedExperiment")){
    args <- c(args,
              list(rowTree = rowTree(x),
                   colTree = colTree(x)))
  }
  ####################################
  ans <- do.call(match.fun(class_FUN), args)
  ####################################
  # class specific post-processing
  if(is(x,"RangedSummarizedExperiment")){
    rowData(ans) <- row_data
  }
  if(is(x,"SingleCellExperiment")){
      int_row_data <- .merge_row_or_col_data(int_elementMetadata(x), f, archetype)
      int_elementMetadata(ans) <- int_row_data
      int_colData(ans) <- int_colData(x)
      int_metadata(ans) <- int_metadata(x)
  }
  if(is(x,"TreeSummarizedExperiment")){
    # optionally merge tree
    row_tree <- rowTree(ans)
    if(!is.null(row_tree) && mergeTree){
      row_leaf <- transNode(tree = row_tree, node = rowLinks(ans)$nodeNum)
      row_tree <- ape::keep.tip(phy = row_tree, tip = row_leaf)
      ans <- changeTree(ans, rowTree = row_tree)
    }
  }
  if(is(x,"MicrobiomeExperiment") && !is.null(referenceSeq(x))){
      referenceSeq <- .merge_row_ranges(referenceSeq(x), f, archetype)
      referenceSeq(ans) <- referenceSeq
  }
  ####################################
  ans
}

#' @importFrom S4Vectors SimpleList metadata
#' @importFrom SummarizedExperiment colData rowRanges rowData rowData<-
#' @importFrom SingleCellExperiment int_colData int_colData<-
#'   int_elementMetadata int_elementMetadata<-
#'   int_metadata int_metadata<-
#' @importFrom TreeSummarizedExperiment rowTree colTree transNode colLinks
#'   changeTree
#' @importFrom ape keep.tip
.merge_cols <- function(x, f, archetype = 1L, mergeTree = FALSE){
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
  col_data <- .merge_row_or_col_data(colData(x), f, archetype)
  # merge assays
  assays <- assays(x)
  assays <- S4Vectors::SimpleList(lapply(assays, .merge_assay, f, "col"))
  names(assays) <- names(assays(x))
  # merge to result
  class_name <- class(x)
  class_FUN <- class_name
  if(grepl("^Ranged",class_FUN)){
    class_FUN <- gsub("^Ranged","",class_FUN)
  }
  args <- list(assays = assays,
               colData = col_data,
               metadata = metadata(x))
  ####################################
  # class specific pre-processing
  if(is(x,"RangedSummarizedExperiment")){
    args <- c(args,
              list(rowRanges = rowRanges(x)))
  } else {
    args <- c(args,
              list(rowData = rowData(x)))
  }
  if(is(x,"TreeSummarizedExperiment")){
    args <- c(args,
              list(rowTree = rowTree(x),
                   colTree = colTree(x)))
  }
  ####################################
  ans <- do.call(match.fun(class_FUN), args)
  ####################################
  # class specific post-processing
  if(is(x,"RangedSummarizedExperiment")){
    rowData(ans) <- rowData(x)
  }
  if(is(x,"SingleCellExperiment")){
    int_col_data <- .merge_row_or_col_data(int_colData(x), f, archetype)
    int_elementMetadata(ans) <- int_elementMetadata(x)
    int_colData(ans) <- int_col_data
    int_metadata(ans) <- int_metadata(x)
  }
  if(is(x,"TreeSummarizedExperiment")){
    # optionally merge tree
    col_tree <- colTree(ans)
    if(!is.null(col_tree) && mergeTree){
      col_leaf <- transNode(tree = col_tree, node = colLinks(ans)$nodeNum)
      col_tree <- ape::keep.tip(phy = col_tree, tip = col_leaf)
      ans <- changeTree(ans, colTree = col_tree)
    }
  }
  ####################################
  ans
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
  function(x, f, archetype = 1L, mergeTree = FALSE){
    .merge_rows(x, f, archetype = archetype, mergeTree = mergeTree)
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
  function(x, f, archetype = 1L, mergeTree = FALSE){
    .merge_cols(x, f, archetype, mergeTree = mergeTree)
  }
)
