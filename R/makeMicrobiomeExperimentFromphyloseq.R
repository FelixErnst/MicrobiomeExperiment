#' Coerce phyloseq object
#'
#' @param obj a \code{phyloseq} object
#'
#' @return An object of class MicrobiomeExperiment
#'
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @examples
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     makeMicrobiomeExperimentFromphyloseq(esophagus)
#' }
makeMicrobiomeExperimentFromphyloseq <- function(obj) {
    # input check
    if(!requireNamespace("phyloseq")){
        stop("'phyloseq' package not found. Please install it to use this ",
             "function.",
             call. = FALSE)
    }
    if(!is(obj,"phyloseq")){
        stop("'obj' must be a 'phyloseq' object")
    }
    #
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
        rowTree = tree,
        microbiomeData = mf
    )
    if(!is.null(obj@sam_data)){
        colData(output) <- DataFrame(data.frame(obj@sam_data))
    }
    return(output)
}
