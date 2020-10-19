#' Agglomerate taxa of the same type.
#'
#' \code{agglomerateByRank} can be used to sum up data based on the association
#' to certain taxonomic ranks given as \code{rowData}. Only available
#' \code{\link{taxonomicRanks}} can be used.
#'
#' @param x \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'
#' @param rank a single character defining a taxonomic rank. Must be a value of
#'   \code{taxonomicRanks()} function.
#'
#' @param ranks a character vector defining taxonomic ranks. All values must be
#'   values of \code{taxonomicRanks()} function.
#'
#' @param onRankOnly \code{TRUE} or \code{FALSE}: Should information only from
#'   the specified rank used or from ranks equal and above?.
#'   (default: \code{onRankOnly = FALSE})
#'
#' @param na.rm \code{TRUE} or \code{FALSE}: Should taxa with an empty rank be
#'   removed? Use it with caution, since result with NA on the selected rank
#'   will be dropped. This setting can be tweaked by defining
#'   \code{empty.fields} to your needs. (default: \code{na.rm = TRUE})
#'
#' @param empty.fields a \code{character} value defining, which values should be
#'   regarded as empty. (Default: \code{c(NA, "", " ", "\t")}). They will be
#'   removed if \code{na.rm = TRUE} before agglomeration.
#'
#' @param agglomerateTree \code{TRUE} or \code{FALSE}: should to
#'   \code{rowTree()} also be agglomerated? (Default:
#'   \code{agglomerateTree = FALSE})
#'
#'
#' @return A taxonomically-agglomerated, optionally-pruned object of the same
#'   class \code{x}.
#'
#' @seealso
#'   \code{\link[=merge-methods]{mergeRows}}
#'
#' @name agglomerate-methods
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' colnames(rowData(GlobalPatterns))
#' taxonomyRanks(GlobalPatterns)
#'
#' # agglomerate at the Family taxonomic rank
#' x1 <- agglomerateByRank(GlobalPatterns, rank="Family")
#' ## How many taxa before/after agglomeration?
#' nrow(GlobalPatterns)
#' nrow(x1)
#'
#' # with agglomeration of the tree
#' x2 <- agglomerateByRank(GlobalPatterns, rank="Family",
#'                         agglomerateTree = TRUE)
#' nrow(x2) # same number of rows, but
#' rowTree(x1) # ... different
#' rowTree(x2) # ... tree
#'
#' ## Look at enterotype dataset...
#' data(enterotype)
#' ## print the available taxonomic ranks. Shows only 1 rank available
#' ## not useful for agglomerateByRank
#' taxonomyRanks(enterotype)
NULL

setGeneric("agglomerateByRank",
           signature = "x",
           function(x, rank = taxonomyRanks(x)[1L], onRankOnly = FALSE,
                    na.rm = FALSE, empty.fields = c(NA, "", " ", "\t", "-"),
                    agglomerateTree = FALSE)
               standardGeneric("agglomerateByRank"))


#' @rdname agglomerate-methods
#' @aliases agglomerateByRank
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("agglomerateByRank", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], onRankOnly = FALSE, na.rm = FALSE,
       empty.fields = c(NA, "", " ", "\t", "-"), agglomerateTree = FALSE){
        # input check
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be an non empty single character value.",
                 call. = FALSE)
        }
        if(!.is_a_bool(onRankOnly)){
            stop("'onRankOnly' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(na.rm)){
            stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
        }
        if(ncol(rowData(x)) == 0L){
            stop("taxonomyData needs to be populated.", call. = FALSE)
        }
        .check_taxonomic_rank(rank, x)
        if(!.is_a_bool(agglomerateTree)){
            stop("'agglomerateTree' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_for_taxonomic_data_order(x)
        #

        # Make a vector from the taxonomic data.
        col <- which( taxonomyRanks(x) %in% rank )
        tax_cols <- .get_tax_cols_from_se(x)

        # if na.rm is TRUE, remove the empty, white-space, NA values from
        # tree will be pruned later, if agglomerateTree = TRUE
        if( na.rm ){
            tax <- as.character(rowData(x)[,tax_cols[col]])
            f <- !(tax %in% empty.fields)
            x <- x[f, , drop=FALSE]
        }

        # get groups of taxonomy entries
        tax_factors <- .get_tax_groups(x, col = col, onRankOnly = onRankOnly)

        # merge taxa
        x <- mergeRows(x, f = tax_factors, mergeTree = agglomerateTree)

        # "Empty" the values to the right of the rank, using NA_character_.
        if( col < length(taxonomyRanks(x)) ){
            badcolumns <- tax_cols[seq_along(tax_cols) > col]
            if(length(badcolumns) > 0L){
                row_data <- rowData(x)
                row_data[, badcolumns] <- NA_character_
                rowData(x) <- row_data
            }
        }
        # adjust rownames
        rownames(x) <- .get_taxonomic_label(x, empty.fields)
        x
    }
)

setGeneric("getAgglomerateData",
           signature = "x",
           function(x, ranks = taxonomyRanks(x)[1:2],
                    empty.fields = c(NA, "", " ", "\t", "-"),
                    rankThresholds = 0)
               standardGeneric("getAgglomerateData"))

#' @rdname agglomerate-methods
#' @aliases getAgglomerateData
#'
#' @param rankThresholds a single numeric value to define threshold for grouping
#'   values in "Others".
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr %>% contains
#' @importFrom tidyr pivot_longer
#'
#' @export
setMethod("getAgglomerateData", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ranks = taxonomyRanks(x)[1:2],
       empty.fields = c(NA, "", " ", "\t", "-"),
       rankThresholds = 0){
        # input checks
        if(!.is_non_empty_character(ranks)){
            stop("'ranks' must be a character vector with no empty values.",
                 call. = FALSE)
        }
        if(ncol(rowData(x)) == 0L){
            stop("taxonomyData needs to be populated.", call. = FALSE)
        }
        .check_taxonomic_ranks(ranks, x)
        .check_for_taxonomic_data_order(x)
        if(!all(.is_numeric_string(rankThresholds))){
            stop("'rankThresholds' must be numeric vector.", call. = FALSE)
        }
        if(length(rankThresholds) != 1L &&
         length(rankThresholds) != length(ranks)){
            stop("'rankThresholds' must have the length == 1L or the same ",
                 "length as 'ranks'.", call. = FALSE)
        }
        if(length(rankThresholds) == 1L){
            rankThresholds <- rep(rankThresholds,length(ranks))
        }
        #
        data <- .get_agglomerate_data(x, ranks)
        data <- .blank_empty_data(data, ranks, empty.fields, rankThresholds)
        ans <- data %>%
          pivot_longer(cols = !contains(ranks), names_to = "Samples",
                       values_to = "RelAbundance")
        ans
    }
)

#' @importFrom SummarizedExperiment colData assays
#' @importFrom dplyr %>% everything mutate
#' @importFrom tibble as_tibble
#' @importFrom rlang :=
.get_agglomerate_data <- function(x, ranks){
    ans <- assays(x)$relabundance %>%
        as_tibble()
    for(i in seq_along(ranks)){
        rank <- rev(ranks)[i]
        ans <- ans %>%
            mutate(!!rank := factor(rowData(x)[,rank]),
                   .before = everything())
    }
    ans
}

#' @importFrom dplyr %>% select contains
.check_for_rank_threshold <- function(data, .rank, ranks, rankThreshold){
    mean(colSums(data %>% select(!contains(ranks)))) >= rankThreshold
}

#' @importFrom dplyr %>% group_by group_map across group_keys pull
.blank_empty_data <- function(data, ranks, empty.fields, rankThresholds){
    data[,ranks] <- lapply(data[,ranks],
                           function(d){
                               d[d %in% empty.fields] <- NA
                               factor(as.character(d),unique(d))
                           })
    na_values <- mapply(
        function(i, rankThreshold){
            rank <- ranks[i]
            previous_ranks <- ranks[seq_along(ranks) < i]
            data_grouped <- data %>%
                group_by(across(c(previous_ranks,rank)))
            threshold_passed <- data_grouped %>%
                group_map(.f = .check_for_rank_threshold, ranks, rankThreshold)
            threshold_passed <- unlist(threshold_passed)
            keys <- data_grouped %>%
                group_keys() %>%
                pull(!!rank)
            unique(as.character(keys[!threshold_passed]))
        },
        seq_along(ranks),
        rankThresholds)
    data[,ranks] <- mapply(
        function(d, na_val){
            d[d %in% na_val] <- NA
            d
        },
        data[,ranks],
        na_values)
    data
}
