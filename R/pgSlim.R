#' @include generics.R
#' @include aaa.R
#' @include pgInMem.R
NULL

#' Class for pangenome data with no reference to genes
#' 
#' This class is a slim version of pgLM and pgFull that does not store any
#' information pertaining to the actual genes. This means that this class cannot
#' be the basis for the creation of a pangenome but that pgLM or pgFull objects
#' can be coerced down to this representation after the pangenome has been 
#' created to make it less burdensome to work with, while still keeping a lot
#' of the functionality of the FindMyFriends framework.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgSlim',
    contains = c('pgInMem')
)

### UTILITY FUNCTIONS
setAs(
    'pgInMem', 'pgSlim',
    function(from) {
        new('pgSlim', 
            seqToOrg = seqToOrg(from),
            seqToGeneGroup = seqToGeneGroup(from),
            groupInfo = groupInfo(from),
            orgInfo = orgInfo(from),
            .settings = from@.settings)
    }
)

#' @describeIn genes Throws error for pgSlim
#' 
setMethod(
    'genes', c('pgSlim', 'missing'),
    function(object, split, subset) {
        stop('pgSlim objects don\'t support gene querying')
    }
)
#' @describeIn genes Throws error for pgSlim
#' 
setMethod(
    'genes', c('pgSlim', 'character'),
    function(object, split, subset) {
        stop('pgSlim objects don\'t support gene querying')
    }
)
#' @describeIn geneNames Throws error for pgSlim
#' 
setMethod(
    'geneNames', 'pgSlim',
    function(object) {
        stop('pgSlim objects don\'t support gene querying')
    }
)
#' @describeIn geneNames Throws error for pgSlim
#' 
setMethod(
    'geneNames<-', 'pgSlim',
    function(object, value) {
        stop('pgSlim objects don\'t support gene querying')
    }
)
#' @describeIn geneWidth Throws error for pgSlim
#' 
setMethod(
    'geneWidth', 'pgSlim',
    function(object) {
        stop('pgSlim objects don\'t support gene querying')
    }
)
