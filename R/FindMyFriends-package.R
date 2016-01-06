#' FindMyFriends: Comparative microbial genomics in R
#' 
#' @details
#' This package has two objectives: Define a framework for working with
#' pangenomic data in R and provide speed and memory effecient algorithms that
#' makes it possible to create huge pangenomes in a reasonable amount of time.
#' While providing novel algorithms itself it also makes it possible to import
#' results from other algorithms into the framework thus facilitating doing
#' post-processing of results from other tools that only provides an initial
#' grouping of genes.
#' 
#' In order to balance speed and memory consumption FindMyFriends provides two
#' different sequence storage modes - either in-memory or as a reference to the 
#' original fasta file. The former excels in lookup speed but can end up too 
#' unwieldy for big pangenomes with Gb of sequence data. The latter in contrast
#' can handle extremely huge sets of genes but can in turn slow down 
#' calculations due to longer sequence lookups.
#' 
#' The novelty of the FindMyFriends algorithms lie primarily in the fact that 
#' they utilise allignment-free sequence comparisons based on cosine similarity
#' of kmer feature vectors. This is substantially faster than BLAST while 
#' retaining the needed resolution. Another novelty is the introduction of 
#' Guided Pairwise Comparison - a different approach than standard all-vs-all
#' comparisons.
#' 
#' @docType package
#' @name FindMyFriends-package
#' @author Thomas Lin Pedersen
#' 
#' @import methods
#' @importFrom Rcpp evalCpp
#' @useDynLib FindMyFriends
#' 
NULL