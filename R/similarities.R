#' @include aaa.R
#' @include pgFull.R
#' @include pgLM.R
NULL

#' @describeIn kmerSimilarity Kmer based similarities for pgVirtual subclasses
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
#' @importFrom kebabs spectrumKernel getExRep linearKernel
#' @importFrom BiocParallel SerialParam bpworkers
#' 
setMethod(
    'kmerSimilarity', 'pgVirtual',
    function(object, lowMem, kmerSize, lowerLimit, rescale, transform, pParam) {
        .fillDefaults(defaults(object))
        
        if (lowMem) {
            if (missing(pParam)) {
                pParam <- SerialParam()
                nSplits <- nOrganisms(object)
            } else {
                nSplits <- bpworkers(pParam)
            }
            res <- lkParallelLM(object, kmerSize, pParam, nSplits, 
                                lowerLimit = lowerLimit)
        } else {
            kernel <- spectrumKernel(kmerSize)
            er <- getExRep(genes(object), kernel, sparse = TRUE)
            if (missing(pParam)) {
                res <- linearKernel(er, sparse = TRUE, diag = FALSE, 
                                    lowerLimit = lowerLimit)
            } else {
                res <- lkParallel(er, pParam, bpworkers(pParam), 
                                  lowerLimit = lowerLimit)
            }
        }
        transformSim(res, lowerLimit, rescale, transform)
    }
)
