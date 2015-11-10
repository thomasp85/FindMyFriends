################################################################################
# TODO: Add matchMembers and matchSim for matching gene groups across Pangenome
#       objects.

#' @include generics.R
#' @include aaa.R
#' @include pgInMem.R
NULL

#' Class for in memory pangenome data
#' 
#' This class handles pangenome data without gene location information and with
#' all sequences stored in memory. This makes sequence lookup much faster but
#' also increases the memory footprint of the object thus making it a bad choice
#' for very large pangenome with millions of genes.
#' 
#' @slot sequences Either an AAStringSet or DNAStringSet containing all 
#' sequences in the pangenome.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
#' @importFrom Biostrings DNAStringSet
#' 
setClass(
    'pgFull',
    contains = 'pgInMem',
    slots = list(
        sequences = 'XStringSet'
    ),
    validity = function(object) {
        if (length(object@sequences) != length(object@seqToOrg)) {
            return('Sequence and index length differ')
        }
        return(TRUE)
    },
    prototype = list(
        sequences = DNAStringSet()
    )
)

### UTILITY FUNCTIONS

#' @describeIn genes Gene access for pgFull and subclasses
#' 
setMethod(
    'genes', c('pgFull', 'missing'),
    function(object, split, subset) {
        if (missing(subset)) {
            object@sequences
        } else {
            object@sequences[subset]
        }
    }
)
#' @describeIn genes Gene access for pgFull and subclasses with group splitting
#' 
setMethod(
    'genes', c('pgFull', 'character'),
    function(object, split, subset) {
        if (!split %in% c('organism', 'group', 'paralogue', 'paralog')) {
            stop('Can only split by organism, gene group or paralogue link')
        }
        if (split == 'organism') {
            ans <- splitStringSet(object@sequences, object@seqToOrg)
            names(ans) <- orgNames(object)[as.integer(names(ans))]
        } else if (split == 'group') {
            if (!hasGeneGroups(object)) {
                stop('No gene groups created')
            }
            ans <- splitStringSet(object@sequences, object@seqToGeneGroup)
            names(ans) <- groupNames(object)[as.integer(names(ans))]
        } else if (split %in% c('paralogue', 'paralog')) {
            if (!hasParalogueLinks(object)) {
                stop('No paralogue links created')
            }
            ans <- splitStringSet(object@sequences, 
                                  paralogueInd(object@seqToGeneGroup, 
                                               groupInfo(object)$paralogue))
            names(ans) <- seq_along(ans)
        }
        if (missing(subset)) {
            ans
        } else {
            ans[subset]
        }
    }
)
#' @describeIn geneNames Get genenames for pgFull and subclasses
#' 
setMethod(
    'geneNames', 'pgFull',
    function(object) {
        names(object@sequences)
    }
)
#' @describeIn geneNames Set genenames for pgFull and subclasses
#' 
setMethod(
    'geneNames<-', 'pgFull',
    function(object, value) {
        names(object@sequences) <- value
        object
    }
)
#' @describeIn geneWidth Get gene widths for pgFull and subclasses
#' 
#' @importFrom Biostrings width
#' 
setMethod(
    'geneWidth', 'pgFull',
    function(object) {
        width(object@sequences)
    }
)
#' @rdname internalMergePangenomes
#' 
setMethod(
    'mergePangenomes', c('pgFull', 'pgFull'),
    function(pg1, pg2, geneGrouping, groupInfo, ...) {
        if (class(pg1) != class(pg2)) {
            stop('pangenomes must be instances of the same class')
        }
        callNextMethod(pg1, pg2, geneGrouping, groupInfo, 
                       sequences = c(pg1@sequences, pg2@sequences), ...)
    }
)