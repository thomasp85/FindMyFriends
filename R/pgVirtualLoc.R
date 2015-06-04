#' @include generics.R
#' @include aaa.R
#' @include pgVirtual.R
NULL

#' Superclass for gene location aware pangenome
#' 
#' This virtual class should be subclassed by all classes that include 
#' chromosomal position of the genes (along with subclassing pgVirtual). The
#' class itself is an empty shell that only takes care of dispatching and 
#' checking the promises of subclasses are held.
#' 
#' Subclasses of pgVirtualLoc must implement the following methods:
#' 
#' \describe{
#'  \item{geneLocation(object)}{Return a data.frame with a row for each gene,
#'  describing the chromosomal position of the gene. The data.frame must contain
#'  the columns 'contig', 'start', 'end' and 'strand'. Contig is 
#'  self-explanatory, start and end is the respective start and end positions on
#'  the contig (start must be lower than end) and strand defines the coding 
#'  direction as -1 or 1.}
#' }
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgVirtualLoc',
    contains=c('VIRTUAL'),
    validity = function(object) {
        if(!hasMethod('geneLocation', class(object))) {
            return('The method "geneLocation" must be implemented for ', class(object))
        }
        gL <- geneLocation(object)
        if(nrow(gL) != nGenes(object)) {
            return('Sequence and sequenceInfo length differ')
        }
        if(!all(c('contig', 'start', 'end', 'strand') %in% names(gL))) {
            return('Missing columns in sequenceInfo')
        }
        if(!all(unique(gL$strand)) %in% c(-1, 1)) {
            return('Strand must be coded with -1 and 1')
        }
        return(TRUE)
    }
)
