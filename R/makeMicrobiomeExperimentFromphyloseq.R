#' Coerce phyloseq object
#'
#' @param obj a \code{phyloseq} object
#'
#' @return An object of class MicrobiomeExperiment
#'
#' @importFrom phyloseq phyloseq
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#' @importClassesFrom phyloseq phyloseq
#'
#' @export
#'
#' @examples
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     as(GlobalPatterns, "MicrobiomeExperiment")
#'     data(enterotype, package="phyloseq")
#'     as(enterotype, "MicrobiomeExperiment")
#'     data(esophagus, package="phyloseq")
#'     as(esophagus, "MicrobiomeExperiment")
#' }
makeMicrobiomeExperimentFromphyloseq <- function(obj) {
    if(!is(obj,"phyloseq")){
        stop("'obj' must be a 'phyloseq' object")
    }
    assays <- SimpleList(counts = obj@otu_table@.Data)
    rowData <- S4Vectors:::make_zero_col_DataFrame(nrow(assays$counts))
    colData <- S4Vectors:::make_zero_col_DataFrame(ncol(assays$counts))
    if(!is.null(obj@tax_table@.Data)){
        rowData <- DataFrame(data.frame(obj@tax_table@.Data))
    }
    if(!is.null(obj@sam_data)){
        colData <- DataFrame(data.frame(obj@sam_data))
    }
    if(!is.null(obj@phy_tree)){
        rowTree <- obj@phy_tree
    } else {
        rowTree <- NULL
    }
    if (!is.null(obj@refseq)) {
        referenceSeq <- obj@refseq
    } else {
        referenceSeq <- NULL
    }
    MicrobiomeExperiment(assays = assays,
                         rowData = obj@tax_table@.Data,
                         colData = colData,
                         rowTree = rowTree,
                         referenceSeq = referenceSeq)
}

setAs("phyloseq", "MicrobiomeExperiment", function(from)
{
    makeMicrobiomeExperimentFromphyloseq(from)
})
