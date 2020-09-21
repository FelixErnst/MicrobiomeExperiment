#' @name taxonomy-methods
#'
#' @title Taxonomy related functions
#'
#' @description
#' These function work on optional data present in \code{rowData}.
#'
#' \code{taxonomyRanks} returns, which columns of \code{rowData(x)} are regarded
#' as columns containing taxonomic information.
#'
#' \code{taxonomyRankEmpty} checks, if a selected rank is empty of information.
#'
#' @param x \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @return
#' \itemize{
#'   \item{\code{taxonomyRanks}:} {a \code{character} vector with all the
#'     taxonomic ranks found in \code{colnames(rowData(x))}}
#'   \item{\code{taxonomyRankEmpty}:} {a \code{logical} value}
#' }
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

setGeneric("taxonomyRanks", signature = c("x"),
           function(x)
             standardGeneric("taxonomyRanks"))

#' @rdname taxonomy-methods
#' @aliases taxonomicRanks
#' @export
setMethod("taxonomyRanks", signature = c(x = "SummarizedExperiment"),
          function(x){
            ranks <- colnames(rowData(x))
            ranks[.get_tax_cols(ranks)]
          }
)

setGeneric("taxonomyRankEmpty",
           signature = "x",
           function(x, rank = taxonomyRanks(x)[1L],
                    empty.fields = c(NA, "", " ", "\t", "-"))
             standardGeneric("taxonomyRankEmpty"))

#' @rdname taxonomy-methods
#' @aliases taxonomyRankEmpty
#' @export
setMethod("taxonomyRankEmpty", signature = c(x = "SummarizedExperiment"),
          function(x, rank = taxonomyRanks(x)[1],
                   empty.fields = c(NA, "", " ", "\t", "-")){
            # input check
            if(!.is_non_empty_string(rank)){
              stop("'rank' must be an non empty single character value.", call. = FALSE)
            }
            if(ncol(rowData(x)) == 0L){
              stop("taxonomyData needs to be populated.", call. = FALSE)
            }
            .check_taxonomic_rank(rank, x)
            .check_for_taxonomic_data_order(x)
            #
            rowData(x)[,rank] %in% empty.fields
          }
)

################################################################################
# helper functions

.get_tax_cols_logical <- function(x){
  tolower(x) %in% TAXONOMY_RANKS
}

.get_tax_cols <- function(x){
  which(.get_tax_cols_logical(x))
}

.get_tax_cols_from_se <- function(x){
  .get_tax_cols(colnames(rowData(x)))
}

.get_tax_groups <- function(x, col, onRankOnly = FALSE){
  tax_cols <- .get_tax_cols_from_se(x)
  tax_col_n <- seq_along(tax_cols)
  if(length(tax_col_n) < col){
    stop(".")
  }
  if(onRankOnly){
    groups <- rowData(x)[,tax_cols[tax_col_n == col],drop=TRUE]
  } else {
    groups <- rowData(x)[,tax_cols[tax_col_n <= col],drop=FALSE]
    groups <- apply(groups,1L,paste,collapse="_")
  }
  factor(groups, unique(groups))
}

#' @importFrom IRanges CharacterList LogicalList
.get_taxonomic_label <- function(x, empty.fields = c(NA, "", " ", "\t", "-"),
                                 with_type = FALSE){
  rd <- rowData(x)
  tax_cols <- .get_tax_cols_from_se(x)
  tax_ranks_non_empty <- !is.na(CharacterList(t(rd[,tax_cols]))) &
    !LogicalList(lapply(CharacterList(t(rd[,tax_cols])),"%in%",empty.fields))
  tax_ranks_non_empty <- t(as(tax_ranks_non_empty,"matrix"))
  tax_ranks_selected <- apply(tax_ranks_non_empty,1L,which)
  if(any(lengths(tax_ranks_selected) == 0L)){
    stop("Only empty taxonomic information detected. Some rows contain only ",
         "entries selected by 'empty.fields'. Cannot generated labels.",
         call. = FALSE)
  }
  if(is.matrix(tax_ranks_selected)){
    tax_ranks_selected <- apply(tax_ranks_selected,2L,max)
  } else if(is.list(tax_ranks_selected)) {
    tax_ranks_selected <- lapply(tax_ranks_selected,max)
    tax_ranks_selected <- unlist(tax_ranks_selected)
  } else if(is.vector(tax_ranks_selected)){
    tax_ranks_selected <- max(tax_ranks_selected)
  } else {
    stop(".")
  }
  tax_cols_selected <- tax_cols[tax_ranks_selected]
  all_same_rank <- length(unique(tax_cols_selected)) == 1L
  ans <- mapply("[",
                as.data.frame(t(as.data.frame(rd))),
                tax_cols_selected,
                SIMPLIFY = FALSE)
  if(with_type || !all_same_rank){
    TR <- toupper(colnames(rd)[tax_cols])
    ans <- paste0(colnames(rd)[unlist(tax_cols_selected)],
                  "_",
                  unlist(ans, use.names = FALSE))
  } else {
    ans <- unlist(ans, use.names = FALSE)
  }
  # last resort - this happens, if annotation data contains ambiguous data
  # sometimes labeled as "circles"
  if(anyDuplicated(ans)){
    dup <- which(ans %in% ans[which(duplicated(ans))])
    ans[dup] <- make.unique(ans[dup], sep = "_")
  }
  ans
}

