#' @include pgSlim.R
#' @include pgInMemLoc.R
#' @include aaa.R
#' @include generics.R
NULL

#' Class for pangenome data with no reference to genes
#' 
#' This class extends \code{\linkS4class{pgSlim}} by subclassing 
#' \code{\linkS4class{pgInMemLoc}} and thus adding gene location information
#' to each gene. See the respective superclasses for more information.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgSlimLoc',
    contains = c('pgSlim', 'pgInMemLoc')
)

setAs(
    'pgInMemLoc', 'pgSlimLoc',
    function(from) {
        new('pgSlimLoc', 
            seqToOrg = seqToOrg(from),
            seqToGeneGroup = seqToGeneGroup(from),
            groupInfo = groupInfo(from),
            orgInfo = orgInfo(from),
            geneLocation = geneLocation(from),
            .settings = from@.settings)
    }
)

#' @rdname internalMergePangenomes
#' 
setMethod(
    'mergePangenomes', c('pgSlimLoc', 'pgSlimLoc'),
    function(pg1, pg2, geneGrouping, groupInfo) {
        if (class(pg1) != class(pg2)) {
            stop('pangenomes must be instances of the same class')
        }
        callNextMethod(pg1, pg2, geneGrouping, groupInfo, 
                       geneLocation = rbind(pg1@geneLocation, 
                                            pg2@geneLocation))
    }
)