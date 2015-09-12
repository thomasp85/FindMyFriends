#' @include pgFull.R
#' @include pgInMemLoc.R
#' @include aaa.R
#' @include generics.R
NULL

#' Class for in memory pangenome data with location information
#' 
#' This class extends \code{\linkS4class{pgFull}} by subclassing 
#' \code{\linkS4class{pgInMemLoc}} and thus adding gene location information
#' to each gene. See the respective superclasses for more information.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgFullLoc',
    contains = c('pgFull', 'pgInMemLoc')
)
#' @rdname internalMergePangenomes
#' 
setMethod(
    'mergePangenomes', c('pgFullLoc', 'pgFullLoc'),
    function(pg1, pg2, geneGrouping, groupInfo) {
        if (class(pg1) != class(pg2)) {
            stop('pangenomes must be instances of the same class')
        }
        callNextMethod(pg1, pg2, geneGrouping, groupInfo, 
                       geneLocation = rbind(pg1@geneLocation, 
                                            pg2@geneLocation))
    }
)