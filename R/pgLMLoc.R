#' @include pgLM.R
#' @include pgInMemLoc.R
#' @include aaa.R
#' @include generics.R
NULL

#' Class for reference based pangenome data with location information
#' 
#' This class extends \code{\linkS4class{pgLM}} by subclassing 
#' \code{\linkS4class{pgVirtualLoc}} and thus adding gene location information
#' to each gene. See the respective superclasses for more information.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgLMLoc',
    contains=c('pgLM', 'pgInMemLoc')
)
#' @rdname internalMergePangenomes
#' 
setMethod(
    'mergePangenomes', c('pgLMLoc', 'pgLMLoc'),
    function(pg1, pg2, geneGrouping, groupInfo) {
        if(class(pg1) != class(pg2)) stop('pangenomes must be instances of the same class')
        pg <- callNextMethod()
        new(
            class(pg1),
            pg,
            sequenceInfo = rbind(pg1@sequenceInfo, pg2@sequenceInfo)
        )
    }
)