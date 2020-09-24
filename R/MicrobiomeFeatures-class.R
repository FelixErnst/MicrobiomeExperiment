## Class constructor
.MicrobiomeFeatures <- setClass("MicrobiomeFeatures",
    contains = "mgFeatures")

#' The MicrobiomeFeatures class
#'
#' Alias for mgFeatures class from metagenomeFeatures package
#'
#' @param taxa taxonomic information as a \code{DataFrame}.
#'   See \code{\link[metagenomeFeatures:mgFeatures-class]{mgFeatures}} for more
#'   details.
#' @param tree phylogenetic tree. See
#'   \code{\link[metagenomeFeatures:mgFeatures-class]{mgFeatures}} for more
#'   details.
#' @param seq reference sequence. See
#'   \code{\link[metagenomeFeatures:mgFeatures-class]{mgFeatures}} for more
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importClassesFrom metagenomeFeatures mgFeatures
#'
#' @name MicrobiomeFeatures-class
#'
#' @aliases MicrobiomeFeatures
#'
#' @export
MicrobiomeFeatures <- function(taxa = DataFrame(),
                               tree = NULL,
                               seq = NULL) {
    if(!is(taxa,"DataFrame")){
        taxa <- as(taxa,"DataFrame")
    }
    mgF <- metagenomeFeatures::mgFeatures(taxa, tree, seq, metadata = list())
    .MicrobiomeFeatures(mgF)
}

setAs("DataFrame", "MicrobiomeFeatures", function(from) {
    MicrobiomeFeatures(from)
})
setAs("list", "MicrobiomeFeatures", function(from) {
    MicrobiomeFeatures(from)
})
