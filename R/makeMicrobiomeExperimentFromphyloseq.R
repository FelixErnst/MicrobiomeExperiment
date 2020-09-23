#' Coerce phyloseq object
#'
#' @param obj a \code{phyloseq} object
#'
#' @return An object of class MicrobiomeExperiment
#'
#' @importFrom phyloseq phyloseq
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment colData
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
    otu <- obj@otu_table@.Data
    taxa <- obj@tax_table@.Data
    if(is.null(taxa)){
        taxa <- matrix(nrow = nrow(otu), ncol=0)
        rownames(taxa) <- rownames(otu)
    }
    if(!is.null(obj@phy_tree)){
        tree <- obj@phy_tree
    } else {
        tree <- NULL
    }
    mf <- MicrobiomeFeatures(taxa = taxa, tree = tree)
    if (!is.null(obj@refseq)) {
        mf@refDbSeq <- obj@refseq
    }
    output <- MicrobiomeExperiment(
        assays = SimpleList(counts = obj@otu_table@.Data),
        rowData = mf
    )
    if(!is.null(obj@sam_data)){
        colData(output) <- DataFrame(data.frame(obj@sam_data))
    }
    return(output)
}

setAs("phyloseq", "MicrobiomeExperiment", function(from)
{
    makeMicrobiomeExperimentFromphyloseq(from)
})
