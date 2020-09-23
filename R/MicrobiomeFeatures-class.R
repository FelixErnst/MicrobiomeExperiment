## Class constructor
.MicrobiomeFeatures <- setClass("MicrobiomeFeatures",
    contains = "mgFeatures")

#' Alias for mgFeatures class from metagenomeFeatures package
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importClassesFrom metagenomeFeatures mgFeatures
#'
#' @aliases MicrobiomeFeatures-class
#' @export
MicrobiomeFeatures <- function(taxa = DataFrame(),
                               tree = NULL,
                               seq = NULL,
                               metadata = list()) {
    mgF <- metagenomeFeatures::mgFeatures(taxa, tree, seq, metadata)
    .MicrobiomeFeatures(mgF)
}

setAs("DataFrame", "MicrobiomeFeatures", function(from) {
    MicrobiomeFeatures(from)
})
setAs("list", "MicrobiomeFeatures", function(from) {
    MicrobiomeFeatures(from)
})
