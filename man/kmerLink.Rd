% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/linking.R
\docType{methods}
\name{kmerLink}
\alias{kmerLink}
\alias{kmerLink,pgVirtual-method}
\title{Link gene groups by homology}
\usage{
kmerLink(object, ...)

\S4method{kmerLink}{pgVirtual}(object, lowMem, kmerSize, lowerLimit,
  rescale, transform, pParam, algorithm, ...)
}
\arguments{
\item{object}{A pgVirtual subclass}

\item{...}{parameters passed on to the community detection algorithm.}

\item{lowMem}{logical. Should low memory footprint be ensured over 
computation speed}

\item{kmerSize}{The size of kmers to use for similarity calculations.}

\item{lowerLimit}{The lower threshold for similarity below which it is set to
0}

\item{rescale}{Should Similarities be normalised between lowerLimit and 1}

\item{transform}{Transformation function to apply to similarities}

\item{pParam}{An optional BiocParallelParam object that defines the workers 
used for parallelisation.}

\item{algorithm}{The name of the community detection algorithm from igraph to
use for gene grouping. See \code{\link[igraph]{communities}} for an overview.
The trailing '.community' can be omitted from the name. Default is 'infomap',
which is also the recommended.}
}
\value{
An object with the same class as object with linking between gene 
groups.
}
\description{
This method allows the user to define a secondary grouping of genes be 
linking gene groups based on sequence similarity (paralogues). A 
representative for each gene group is used for the calculations and the 
similarity is assessed using the kmer based cosine similarity.
}
\section{Methods (by class)}{
\itemize{
\item \code{pgVirtual}: Linking for pgVirtual subclasses
}}

\examples{
testPG <- .loadPgExample(withGroups=TRUE)

# No paralogue links
hasParalogueLinks(testPG)

# Create the links
testPG <- kmerLink(testPG)

}
