#' @include aaa.R
#' @include generics.R
#' @include grouping.R
#' @include pgVirtual.R
NULL

#' @describeIn kmerLink Linking for pgVirtual subclasses
#' 
#' @param lowMem logical. Should low memory footprint be ensured over 
#' computation speed
#' 
#' @param kmerSize The size of kmers to use for similarity calculations.
#' 
#' @param lowerLimit The lower threshold for similarity below which it is set to
#' 0
#' 
#' @param rescale Should Similarities be normalised between lowerLimit and 1
#' 
#' @param transform Transformation function to apply to similarities
#' 
#' @param pParam An optional BiocParallelParam object that defines the workers 
#' used for parallelisation.
#'   
#' @param nSplits The number of jobs to split the computations into
#' 
#' @param algorithm The name of the community detection algorithm from igraph to
#' use for gene grouping. See \code{\link[igraph]{communities}} for an overview.
#' The trailing '.community' can be omitted from the name. Default is 'infomap',
#' which is also the recommended.
#' 
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' @importFrom BiocParallel SerialParam
#' 
setMethod(
    'kmerLink', 'pgVirtual',
    function(object, lowMem, kmerSize, lowerLimit, rescale, transform, pParam, nSplits, algorithm, ...) {
        .fillDefaults(defaults(object))
        
        gRep <- getRep(object, 'random')
        er <- getExRep(gRep, spectrumKernel(kmerSize), sparse = TRUE)
        if(missing(pParam)) {
            sim <- linearKernel(er, sparse=TRUE, diag=FALSE, lowerLimit=lowerLimit)
        } else {
            sim <- lkParallel(er, pParam, nSplits, lowerLimit=lowerLimit)
        }
        sim <- transformSim(sim, lowerLimit, rescale, transform)
        members <- igGroup(sim, algorithm, ...)
        groupInfo(object)$paralogue <- members
        object
    }
)