#' @include aaa.R
#' @include generics.R
#' @include pgFull.R
#' @include pgLM.R
#' @include pgFullLoc.R
#' @include pgLMLoc.R
NULL

#' Construct a pangenome from fasta files
#' 
#' This function constructs an initial pangenome object from a set of fasta 
#' files. Note that the actual pangenome is not calculated here. As such this
#' function mainly sets everything up before beginning the more lengthly
#' pangenome calculation.
#' 
#' @param paths A character vector with location of fasta files
#' 
#' @param translated A boolean indicating if the fasta files contain amino acid
#' sequences
#' 
#' @param geneLocation A function, string or dataframe. If it is a data.frame it
#' should contain the columns 'contig', 'start', 'end' and 'strand' with a row
#' for each gene. If it is a function it should take the name (fasta 
#' description) for each gene and output a data.frame similar to described 
#' above. If it is a string it should specify the format of the gene names. 
#' Currently only 'prodigal' is supported.
#' 
#' @param lowMem Boolean. Should FindMyFriends avoid storing sequences in 
#' memory.
#' 
#' @param ... Additional defaults to set on the object
#' 
#' @return A pgVirtual subclass object depending on geneLocation and lowMem.
#' \tabular{lll}{
#'  \bold{geneLocation} \tab \bold{lowMem} \tab \bold{Resulting class}         \cr
#'  NULL                \tab FALSE         \tab \code{\linkS4class{pgFull}}    \cr
#'  NULL                \tab TRUE          \tab \code{\linkS4class{pgLM}}      \cr
#'  !NULL               \tab FALSE         \tab \code{\linkS4class{pgFullLoc}} \cr
#'  !NULL               \tab TRUE          \tab \code{\linkS4class{pgLMLoc}}   \cr
#' }
#' 
#' @examples 
#' location <- tempdir()
#' unzip(system.file('extdata', 'Mycoplasma.zip', package='FindMyFriends'),
#'       exdir=location)
#' genomeFiles <- list.files(location, full.names=TRUE, pattern='*.fasta')
#' 
#' # Create pgFull
#' pangenome(genomeFiles, TRUE)
#' 
#' # Create pgFullLoc
#' pangenome(genomeFiles, TRUE, geneLocation='prodigal')
#' 
#' # Create pgLM
#' pangenome(genomeFiles, TRUE, lowMem=TRUE)
#' 
#' # Create pgLMLoc
#' pangenome(genomeFiles, TRUE, geneLocation='prodigal', lowMem=TRUE)
#' 
#' @export
#' 
#' @importFrom Biostrings readAAStringSet readDNAStringSet fasta.index
#' @importFrom tools file_path_sans_ext
#' 
pangenome <- function(paths, translated, geneLocation = NULL, lowMem = FALSE, 
                      ...) {
    settings <- .pkg_variables$defaults
    sOverwrite <- list(...)
    for (i in names(sOverwrite)) {
        settings[[i]] <- sOverwrite[[i]]
    }
    settings$translated <- translated
    settings$lowMem <- lowMem
    args <- list(.settings = settings)
    args$Class <- if (lowMem) {
        if (is.null(geneLocation)) {
            'pgLM'
        } else {
            'pgLMLoc'
        }
    } else {
        if (is.null(geneLocation)) {
            'pgFull'
        } else {
            'pgFullLoc'
        }
    }
    sequenceFileLength <- as.integer(nSeqs(paths))
    if (lowMem) {
        args$seqIndex <- fasta.index(paths, 
                                     seqtype = if (translated) 'AA' else 'DNA')
    } else {
        args$sequences <- if (translated) {
            readAAStringSet(paths)
        } else {
            readDNAStringSet(paths)
        }
    }
    if (any(sequenceFileLength == 0)) {
        warning('The following files contained no sequences\n\n', 
                paste(paths[sequenceFileLength == 0], collapse = '\n'))
    }
    args$seqToOrg <- as.integer(unlist(sapply(1:length(paths), function(x) {
        rep(x, sequenceFileLength[x])
    })))
    
    orgNames <- file_path_sans_ext(basename(paths))
    args$orgInfo <- data.frame(
        nGenes = sequenceFileLength, 
        row.names = orgNames, 
        check.names = FALSE, 
        stringsAsFactors = FALSE)
    args$matrix <- matrix(nrow = 0, ncol = length(paths), 
                          dimnames = list(NULL, orgNames))
    if (!is.null(geneLocation)) {
        geneNames <- if (lowMem) {
            args$seqIndex$desc
        } else {
            names(args$sequences)
        }
        args$geneLocation <- getSeqInfo(geneLocation, geneNames)
    }
    do.call(new, args)
}

.pkg_variables$defaults <- list(
    groupPrefix = 'OG',
    nextGroup = 1,
    kmerSize = 4,
    lowerLimit = 0.5,
    algorithm = 'infomap',
    flankSize = 4,
    minFlank = 0,
    forceParalogues = TRUE,
    rescale = TRUE,
    transform = FALSE,
    maxLengthDif = 0.1
)