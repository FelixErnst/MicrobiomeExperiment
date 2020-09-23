################################################################################
# testing

.is_a_bool <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
  is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
  .is_non_empty_character(x) && length(x) == 1L
}

.is_a_string <- function(x){
  is.character(x) && length(x) == 1L
}

.are_whole_numbers <- function(x){
  tol <- 100 * .Machine$double.eps
  abs(x - round(x)) <= tol && !is.infinite(x)
}

.is_numeric_string <- function(x){
  x <- as.character(x)
  suppressWarnings({x <- as.numeric(x)})
  !is.na(x)
}

.is_function <- function(x){
  typeof(x) == "closure" && is(x, "function")
}

.all_are_existing_files <- function(x){
  all(file.exists(x))
}

.get_name_in_parent <- function(x) {
  .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
  paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}

################################################################################
# checks

.check_altVal <- function(x, abund_values, altVal,
                          abund_valuesName = .get_name_in_parent(abund_values),
                          altValName = .get_name_in_parent(altVal)){
  if(ncol(x) != ncol(altVal)){
    stop("Number of columns of 'x' and '",altVal,"' does not match.",
         call. = FALSE)
  }
  if(is(altVal,"SummarizedExperiment")){
    .check_abund_values(abund_values, altVal)
  }
}

#' @importFrom SummarizedExperiment assays
.check_abund_values <- function(abund_values, x,
                                name = .get_name_in_parent(abund_values)){
  if(!.is_non_empty_string(abund_values)){
    stop("'",name,"' must be a single non-empty character value.",
         call. = FALSE)
  }
  if(!(abund_values %in% names(assays(x)))){
    stop("'",name,"' must be a valid name of assays(x)", call. = FALSE)
  }
}

.check_relative_abundance <- function(data){
  if(sum(colSums(data)) != ncol(data)){
    warning("The selected assay probably does not contain relative abundance ",
            "data.", call. = FALSE)
  }
}

#' @importFrom SummarizedExperiment colData
.check_var <- function(var, x, name = .get_name_in_parent(var)){
  if(missing(var)
     || !.is_non_empty_character(var)
     || !all(var %in% colnames(colData(x)))){
    stop("'",name,"' must contain valid colnames of colData(x)", call. = FALSE)
  }
}

#' @importFrom SummarizedExperiment colData
.check_single_var <- function(var, x, name = .get_name_in_parent(var)){
  if(missing(var)
     || !.is_non_empty_string(var)
     || !(var %in% colnames(colData(x)))){
    stop("'",name,"' must be a single colname of colData(x)", call. = FALSE)
  }
}

.check_var_data <- function(data, name){
  if(is.character(data) && !is.factor(data)){
    warning("variable data '",name,"' is a character and not a factor.
            Converting ...\n",
            call. = FALSE)
    data <- as.factor(data)
  }
  if(is.factor(data) && length(levels(data)) < 2L || length(unique(data)) < 2L){
    warning("variable data '",name,"' with only one level (all values are the ",
            "same). Removing from correlation calculation.",
            call. = FALSE)
    return(NULL)
  }
  data
}

.check_base_var <- function(base_var, var_data){
  if(length(base_var) > 1L ||
     (!is.integer(base_var)
      && (!is.numeric(base_var) || base_var != as.integer(base_var))
      && !(base_var %in% levels(var_data)))){
    stop("'base_var' must be a single integer or valid name for a level of ",
         "the select 'var'", call. = FALSE)
  }
}

.check_taxonomic_rank <- function(rank, x){
  if( !(rank %in% taxonomyRanks(x) ) ){
    stop("'rank' must be a value from 'taxonomyRanks()'")
  }
}
.check_taxonomic_ranks <- function(ranks, x){
  if( !all(ranks %in% taxonomyRanks(x) ) ){
    stop("'ranks' must contain values from 'taxonomyRanks()'")
  }
}

#' @importFrom SummarizedExperiment rowData
.check_for_taxonomic_data_order <- function(x){
  ranks <- colnames(rowData(x))
  f <- tolower(ranks) %in% TAXONOMY_RANKS
  if(!any(f)){
    stop("no taxonomic ranks detected in rowData(). Colnames named on of the ",
         "following values can be used: '",
         paste(TAXONOMY_RANKS, collapse = "', '"), "'", call. = FALSE)
  }
  m <- match(TAXONOMY_RANKS, tolower(ranks[f]))
  m <- m[!is.na(m)]
  # check that taxonomic ranks are in order. If they are all value in check
  # should be 1 or 0
  check <- unique(c(m[-1], m[length(m)]) - m )
  if(!all(check %in% c(1L,0L))){
    stop("Taxonomic ranks are not in order. Please reorder columns, which ",
         "correspond to taxonomic ranks like this:\n'",
         paste(TAXONOMY_RANKS, collapse = "', '"), "'.",
         call. = FALSE)
  }
}

################################################################################
# low-level accessors

.get_correlation_row_data_id <- function(abund_values,method,var,base_var){
  paste0(abund_values,"_",method,"_",var,"_",base_var)
}

.get_base_var <- function(var_data, base_var){
  if(is.numeric(base_var)){
    base_var <- levels(var_data)[[base_var]]
  } else {
    base_var <- levels(var_data)[levels(var_data) %in% base_var]
  }
  base_var
}
