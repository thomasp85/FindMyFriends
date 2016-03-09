#' @include generics.R
#' @include aaa.R
#' @include pgInMem.R
NULL

#' Class for reference based pangenome data
#' 
#' This class handles pangenome information where gene sequences are kept on 
#' disc instead of stored in memory. As long as the original fasta files are
#' not modified, this class will take care of indexing the genes correctly. This
#' class has a substantially lower memory footprint than the 
#' \code{\linkS4class{pgFull}} class at the expense of longer sequence lookup
#' times. For massive pangenomes containing Gb of sequence data there is no 
#' alternative though.
#' 
#' @slot seqIndex A data.frame as produced by 
#' \code{\link[Biostrings]{fasta.index}} with random access information for each
#' gene.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgLM',
    contains = 'pgInMem',
    slots = list(
        seqIndex = 'data.frame'
    ),
    validity = function(object) {
        if (!all(names(object@seqIndex) %in% c('recno', 
                                               'fileno', 
                                               'offset', 
                                               'desc', 
                                               'seqlength', 
                                               'filepath'))) {
            return('seqIndex must be a valid indexing data.frame for 
                   Biostrings')
        }
        if (length(object@seqToOrg) != nrow(object@seqIndex)) {
            return('Gene indexes of different length')
        }
        return(TRUE)
    },
    prototype = list(
        seqIndex = data.frame(
            recno = integer(), 
            fileno = integer(), 
            offset = integer(), 
            desc = character(), 
            seqlength = integer(), 
            filepath = character()
        )
    )
)

### UTILITY FUNCTIONS

#' @describeIn genes Gene access for pgLM and subclasses
#' 
setMethod(
    'genes', c('pgLM', 'missing'),
    function(object, split, subset) {
        if (missing(subset)) {
            if (translated(object)) {
                safeAAread(object@seqIndex)
            } else {
                safeDNAread(object@seqIndex)
            }
        } else {
            if (translated(object)) {
                safeAAread(object@seqIndex[subset,]) 
            } else {
                safeDNAread(object@seqIndex[subset,])
            }
        }
    }
)
#' @describeIn genes Gene access for pgFull and subclasses with group splitting
#' 
setMethod(
    'genes', c('pgLM', 'character'),
    function(object, split, subset) {
        if (!split %in% c('organism', 'group', 'paralogue', 'paralog')) {
            stop('Can only split by organism, gene group or paralogue link')
        }
        if (split == 'organism') {
            if (missing(subset)) subset <- 1:nOrganisms(object)
            if (inherits(subset, 'character')) {
                subset <- match(subset, orgNames(object))
            }
            seqSubset <- which(object@seqToOrg %in% subset)
            ans <- genes(object, subset=seqSubset)
            ans <- splitStringSet(ans, object@seqToOrg[seqSubset])
            names(ans) <- orgNames(object)[as.integer(names(ans))]
        } else if (split == 'group') {
            if (!hasGeneGroups(object)) {
                stop('No gene groups created')
            }
            if (missing(subset)) subset <- 1:nGeneGroups(object)
            if (inherits(subset, 'character')) {
                subset <- match(subset, groupNames(object))
            }
            seqSubset <- which(object@seqToGeneGroup %in% subset)
            ans <- genes(object, subset=seqSubset)
            ans <- splitStringSet(ans, object@seqToGeneGroup[seqSubset])
            names(ans) <- groupNames(object)[as.integer(names(ans))]
        } else if (split %in% c('paralogue', 'paralog')) {
            if (!hasParalogueLinks(object)) {
                stop('No paralogue links created')
            }
            seqToPar <- paralogueInd(object@seqToGeneGroup, 
                                     groupInfo(object)$paralogue)
            if (missing(subset)) subset <- 1:max(seqToPar)
            seqSubset <- which(seqToPar %in% subset)
            ans <- genes(object, subset=seqSubset)
            ans <- splitStringSet(ans, seqToPar[seqSubset])
            names(ans) <- seq_along(ans)
        }
        ans
    }
)
#' @describeIn geneNames Get genenames for pgLM and subclasses
#' 
setMethod(
    'geneNames', 'pgLM',
    function(object) {
        object@seqIndex$desc
    }
)
#' @describeIn geneNames Set genenames for pgLM and subclasses
#' 
setMethod(
    'geneNames<-', 'pgLM',
    function(object, value) {
        object@seqIndex$desc <- value
        object
    }
)
#' @describeIn geneWidth Get gene width for pgLM and subclasses
#' 
setMethod(
    'geneWidth', 'pgLM',
    function(object) {
        object@seqIndex$seqlength
    }
)
#' @rdname internalMergePangenomes
#' 
setMethod(
    'mergePangenomes', c('pgLM', 'pgLM'),
    function(pg1, pg2, geneGrouping, groupInfo, ...) {
        if (class(pg1) != class(pg2)) {
            stop('pangenomes must be instances of the same class')
        }
        callNextMethod(pg1, pg2, geneGrouping, groupInfo, 
                       seqIndex = rbind(pg1@seqIndex, pg2@seqIndex), ...)
    }
)

## HELPERS

#' @importFrom Biostrings readAAStringSet
#' 
safeAAread <- function(index) {
    if (length(unique(index$fileno)) > 2000) {
        nIndex <- index[order(index$fileno), ]
        files <- unique(nIndex$fileno)
        splits <- rep(seq_len(ceiling(length(files)/2000)), 
                      each = 2000, length.out = length(files))
        seqs <- lapply(split(files, splits), function(f) {
            readAAStringSet(nIndex[nIndex$fileno %in% f, ])
        })
        seqs <- Reduce(c, seqs)
        seqs[match(index$recno, nIndex$recno)]
    } else {
        readAAStringSet(index)
    }
}

#' @importFrom Biostrings readDNAStringSet
#' 
safeDNAread <- function(index) {
    if (length(unique(index$fileno)) > 2000) {
        nIndex <- index[order(index$fileno), ]
        files <- unique(nIndex$fileno)
        splits <- rep(seq_len(ceiling(length(files)/2000)), 
                      each = 2000, length.out = length(files))
        seqs <- lapply(split(files, splits), function(f) {
            readDNAStringSet(nIndex[nIndex$fileno %in% f, ])
        })
        seqs <- Reduce(c, seqs)
        seqs[match(index$recno, nIndex$recno)]
    } else {
        readDNAStringSet(index)
    }
}