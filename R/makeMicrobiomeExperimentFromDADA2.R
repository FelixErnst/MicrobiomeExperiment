#' Coerce DADA2 results to \code{MicrobiomeExperiment}
#'
#' \code{makeMicrobiomeExperimentFromDADA2} is a wrapper for the
#' \code{mergePairs} function from the \code{dada2 package}.
#'
#' @param ... See \code{mergePairs} function for
#'   more details.
#'
#' @return An object of class \code{MicrobiomeExperiment}
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom Biostrings DNAStringSet
#' @importFrom stringr str_pad
#'
#' @export
makeMicrobiomeExperimentFromDADA2 <- function(...) {
    # input checks
    if(!requireNamespace("dada2")){
        stop("'dada2' package not found. Please install it to use this ",
             "function.",
             call. = FALSE)
    }
    #
    mergers <- dada2::mergePairs(...)
    seqtab <- dada2::makeSequenceTable(mergers)
    seqtab <- t(seqtab)
    # generate row and col names
    rName <- paste0("ASV",
                    stringr::str_pad(seq.int(1L,nrow(seqtab)),
                                     nchar(nrow(seqtab)) + 1L,
                                     pad="0"))
    cName <- colnames(seqtab)
    # retrieve count data and reference sequence
    assays <- S4Vectors::SimpleList(counts = unname(seqtab))
    refseq <- Biostrings::DNAStringSet(rownames(seqtab))
    # construct ME an name rows and cols
    output <- MicrobiomeExperiment(assays = assays,
                                   refSeq = refseq)
    colnames(output) <- cName
    rownames(output) <- rName
    output
}
