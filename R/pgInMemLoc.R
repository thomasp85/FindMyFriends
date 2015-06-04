#' @include pgVirtualLoc.R
NULL

#' Superclass for gene location aware pangenome
#' 
#' This virtual class is the superclass for all standard, location aware, 
#' pangenome classes in FindMyFriends. It stores all chromosomal information in
#' a data.frame.
#' 
#' @slot geneLocation A data.frame containing the columns 'contig', 'start', 
#' 'end' and 'strand' and a row for each gene in the pangenome.
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgInMemLoc',
    contains=c('VIRTUAL', 'pgVirtualLoc'),
    slots=list(
        geneLocation='data.frame'
    ),
    prototype = list(
        sequenceInfo = data.frame(contig=character(), start=integer(), end=integer(), strand=integer())
    )
)
#' @describeIn geneLocation Get gene location for pgInMemLoc subclasses
#' 
setMethod(
    'geneLocation', 'pgInMemLoc',
    function(object) {
        object@geneLocation
    }
)