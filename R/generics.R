#' Access default values for a pgVirtual subclass object
#' 
#' This method lets the user view and set the default values used for the 
#' different algorithms in FindMyFriends. Many of the parameters are reoccuring
#' and it can become laborious to type them in at each step. These 
#' functionalities makes it easy to set defaults on a per-pangenome basis.
#' 
#' Currently the following methods support reading defaults from a pgVirtual
#' object. Note that only directly named arguments are supported - arguments
#' passed on through the \code{...}-mechanism are not supported unless they are
#' passed to a function that support it.
#' 
#' \itemize{
#'   \item \code{\link{graphGrouping}}
#'   \item \code{\link{gpcGrouping}}
#'   \item \code{\link{variableRegions}}
#'   \item \code{\link{plotGroup}}
#'   \item \code{\link{kmerLink}}
#'   \item \code{\link{plotSimilarity}}
#'   \item \code{\link{plotTree}}
#'   \item \code{\link{kmerSimilarity}}
#' }
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A named list of default values
#' 
#' @examples 
#' # Get all object defaults
#' testPG <- .loadPgExample()
#' defaults(testPG)
#' 
#' # Set a new default
#' defaults(testPG)$minFlank <- 2
#' 
#' @export
#' 
setGeneric("defaults", function(object) {
    standardGeneric("defaults")
})
#' @rdname defaults
#' 
#' @param value The new values to set
#' 
#' @export
#' 
setGeneric("defaults<-", function(object, value) {
    standardGeneric("defaults<-")
})
#' Get the number of organisms represented in a pangenome
#' 
#' This method returns the current number of organisms in a pgVirtual 
#' subclass. This is also the result of calling \code{length()} on the object.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer giving the number of organisms
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' nOrganisms(testPG)
#' 
#' @export
#' 
setGeneric( 'nOrganisms', def = function(object) {
    standardGeneric('nOrganisms')
})
#' Get the total number of genes in a pangenome
#' 
#' This method returns the total number of genes in a pangenome (i.e. the sum
#' of genes in each organism in the pangenome)
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer giving the number of genes in the object
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' nGenes(testPG)
#' 
#' @export
#' 
setGeneric('nGenes', def = function(object) {
    standardGeneric('nGenes')
})
#' Get the number of gene groups in a pangenome
#' 
#' This method gives the number of different gene groups in the object.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer giving the number of gene groups
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' nGeneGroups(testPG)
#' 
#' @export
#' 
setGeneric('nGeneGroups', def = function(object) {
    standardGeneric('nGeneGroups')
})
#' Check whether gene groups are defined
#' 
#' This method checks whether any grouping of genes has been done on the 
#' object and returns TRUE if that is the case.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A boolean indicating whether gene groups have been defined (TRUE) or
#' not (FALSE)
#' 
#' @examples 
#' # Empty pangenome
#' testPG <- .loadPgExample()
#' hasGeneGroups(testPG)
#' 
#' # With gene groups
#' testPG <- .loadPgExample(withGroups=TRUE)
#' hasGeneGroups(testPG)
#' 
#' @export
#' 
setGeneric('hasGeneGroups', def = function(object) {
    standardGeneric('hasGeneGroups')
})
#' Checks whether linking of paralogues has been done
#' 
#' This method checks for the existance of paralogue links in the object.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A boolean indicating whether paralogue links have been defined (TRUE)
#' or not (FALSE)
#' 
#' @examples 
#' # No paralogues
#' testPG <- .loadPgExample(withGroups=TRUE)
#' hasParalogueLinks(testPG)
#' 
#' # With paralogues
#' testPG <- .loadPgExample(withGroups=TRUE, withParalogues=TRUE)
#' hasParalogueLinks(testPG)
#' 
#' @export
#' 
setGeneric('hasParalogueLinks', def = function(object) {
    standardGeneric('hasParalogueLinks')
})
#' Checks for existance of gene location information
#' 
#' This method checks whether gene location information is present in the 
#' object i.e. if the object inherits from pgVirtualLoc
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A boolean indicating whether gene location information is present 
#' (TRUE) or not (FALSE)
#' 
#' @examples 
#' # Exclusive pgVirtual subclasses
#' testPG <- .loadPgExample()
#' hasGeneInfo(testPG)
#' 
#' # pgVirtualLoc subclasses
#' testPG <- .loadPgExample(geneLoc=TRUE)
#' hasGeneInfo(testPG)
#' 
#' @export
#' 
setGeneric('hasGeneInfo', def = function(object) {
    standardGeneric('hasGeneInfo')
})
#' Get gene location for all genes
#' 
#' This method returns the gene location of all genes as a data.frame with each
#' row corresponding to a gene in the pangenome. The data.frame will have the 
#' columns 'start', 'end', 'contig' and 'strand' (order of columns not ensured)
#' with start and end giving the start and end position of the gene on the 
#' contig/chromosome given in the contig column. Strand gives the direction of
#' translation, 1 is from start to end and -1 is from end to start (thus start 
#' should always be lower than end no matter the direction of translation)
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A data.frame as described above
#' 
#' @note Required for subclasses of pgVirtualLoc in order to extend the class
#' system of FindMyFriends
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE)
#' head(geneLocation(testPG))
#' 
#' @export
#' 
setGeneric('geneLocation', def = function(object) {
    standardGeneric('geneLocation')
})
#' Check the sequence type of the pangenome
#' 
#' This method checks whether the genes in the pangenome are on translated
#' form (amino acid sequences) or not. A return value of FALSE only indicates
#' that the storage mode for the genes is not an AAStringSet. While this leaves
#' room for both RNA-, DNA- and BStringSet, only DNAStringSet makes much sense
#' and is therefore assumed
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A boolean indicating whether genes are translated (TRUE) or not 
#' (FALSE)
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Genes are translated
#' translated(testPG)
#' 
#' # ... and therefore returned as AAStringSet instead of DNAStringSet
#' class(genes(testPG, subset=1))
#' 
#' @export
#' 
setGeneric('translated', def = function(object) {
    standardGeneric('translated')
})
#' Extract gene sequences from a pangenome
#' 
#' This method is used to extract the genomic sequences that is the basis for
#' the pangenome. Genes can be split and subsetted upfront based on other 
#' information in the pangenome, such as gene groups and organisms. For some
#' pgVirtual subclasses the subset parameter is mandatory in order to avoid 
#' reading all genes into memory at once.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param split A string giving the optional splitting type. Either 'organism', 
#' 'group' or 'paralogue'.
#' 
#' @param subset A subsetting of the result equal to using '[]' on the result.
#' It is generally recommended to use this instead of subsetting the result, as
#' it avoids unneeded memory allocation.
#' 
#' @return An XStringSet if split is missing or an XStringSetList if it is 
#' present
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE, withParalogues=TRUE)
#' # Direct gene access
#' genes(testPG)
#' 
#' # Early subsetting
#' genes(testPG, subset=1:10)
#' 
#' # Split by membership
#' genes(testPG, split='organism')
#' genes(testPG, split='group')
#' genes(testPG, split='paralogue')
#' 
#' # Split and subset - get genes from the first organism
#' genes(testPG, split='organism', subset=1)
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('genes', def = function(object, split, subset) {
    standardGeneric('genes')
})
#' Get a representative sequence for each gene group
#' 
#' This method returns a representative sequence for each of the gene groups
#' defined in the pangenome. Currently the methods defined for selecting 
#' sequences are 'random', 'shortest', and 'longest. In case of tie for the two 
#' latter the first occurence gets returned. Consensus sequence might be added
#' at a latter stage.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param method The method to use to get a representative. Either 'random',
#' 'shortest' or 'longest'.
#' 
#' @return An XStringSet
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Get a random sequence from each group
#' getRep(testPG, 'random')
#' 
#' @export
#' 
setGeneric('getRep', def = function(object, method) {
    standardGeneric('getRep')
})
#' Get and set the names of the genes in the pangenome
#' 
#' These methods lets you query and change the naming of genes in your 
#' pangenome. Take note that even though sequences are not in memory for pgLM
#' objects, the names are. This means that changes to the description header in
#' the underlying fasta files have no effect on the naming in your pangenome
#' 
#' @param object A pgVirtual subclass
#' 
#' @return In case of the getter, a character vector containing the names of 
#' each gene.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' head(geneNames(testPG))
#' 
#' geneNames(testPG)[10] <- 'Gene number 10'
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('geneNames', def = function(object) {
    standardGeneric('geneNames')
})
#' @rdname geneNames
#' 
#' @param value A character vector with new names
#' 
#' @export
#' 
setGeneric('geneNames<-', def = function(object, value) {
    standardGeneric('geneNames<-')
})
#' Get the sequence length of each gene
#' 
#' This method extracts the width (i.e. number of residues) of each gene in
#' the pangenome.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer vector with the length of each sequence
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' head(geneWidth(testPG))
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('geneWidth', def = function(object) {
    standardGeneric('geneWidth')
})
#' Get and set the names of organisms in the pangenome
#' 
#' These methods lets you manipulate the naming of organisms in the pangenome.
#' By default organisms are named after the fasta file they are defined by, but
#' this can be changed at will.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return In case of the getter a character vector with names
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' orgNames(testPG)
#' 
#' orgNames(testPG)[3] <- 'Organism 3'
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('orgNames', def = function(object) {
    standardGeneric('orgNames')
})
#' @rdname orgNames
#' 
#' @param value A vector with new names - will be coerced to characters
#' 
#' @export
#' 
setGeneric('orgNames<-', def = function(object, value) {
    standardGeneric('orgNames<-')
})
#' Get and set the names of gene groups in the pangenome
#' 
#' These methods lets you manipulate the naming of gene groups in the 
#' pangenome. By default organisms are numbered consecutively but this can be
#' changed at will. New gene groups will be numbered though despite what naming
#' scheme has been introduced before.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return In case of the getter a character vector with names
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' head(groupNames(testPG))
#' 
#' groupNames(testPG)[20] <- 'Gene group 20'
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('groupNames', def = function(object) {
    standardGeneric('groupNames')
})
#' @rdname groupNames
#' 
#' @param value A vector with new names - will be coerced to characters
#' 
#' @export
#' 
setGeneric('groupNames<-', def = function(object, value) {
    standardGeneric('groupNames<-')
})
#' Get and set information about organisms
#' 
#' These methods lets you access the information stored about each organism and
#' add to it or modify it. The only information present up front is the number 
#' of genes present in each organism. While possible, this information should 
#' not be changed manually but through the \code{\link{removeGene}} functions.
#' 
#' @param object A pgVirtual subclass
#'   
#' @return In case of the getter a data.frame with organism information.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' orgInfo(testPG)
#' 
#' orgInfo(testPG)$Genus <- 'Mycoplasma'
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @family Metadata
#' 
#' @export
#' 
setGeneric('orgInfo', def = function(object) {
    standardGeneric('orgInfo')
})
#' @rdname orgInfo
#' 
#' @param value A data.frame with a row for each organism
#' 
#' @export
#' 
setGeneric('orgInfo<-', def = function(object, value) {
    standardGeneric('orgInfo<-')
})
#' Get and set information about gene group
#' 
#' These methods lets you access the information stored about each gene group 
#' and add to it or modify it. Upfront the following columns are present: 
#' 'description', 'group', 'paralogue', 'GO', 'EC', 'nOrg' and 'nGenes'. All 
#' except 'group', 'nOrg' and 'nGenes' are filled with NA as default. The latter
#' are prefilled with information derived from the grouping itself and should 
#' not be modified manually. 'description' is meant to contain a human readable
#' description of the functionality of the gene group, 'GO' should contain GO 
#' terms (stored in a list of character vectors) and EC should contain enzyme
#' numbers (again stored as a list of character vectors). There is no check for 
#' the validity of the content so it is up to the user to ensure that the terms
#' added are valid. Additional columns can be added at will.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return In case of the getter a data.frame with organism information.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' head(groupInfo(testPG))
#' 
#' groupInfo(testPG)$description[1] <- 'transposase'
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @family Metadata
#' 
#' @export
#' 
setGeneric('groupInfo', def = function(object) {
    standardGeneric('groupInfo')
})
#' @rdname groupInfo
#' 
#' @param value A data.frame with a row for each group
#' 
#' @export
#' 
setGeneric('groupInfo<-', def = function(object, value) {
    standardGeneric('groupInfo<-')
})
#' Get the pangenome matrix
#' 
#' This method lets you extract the pangenome matrix of the pangenome. It is not
#' possible to directly change the pangenome matrix as it not necessary stored 
#' in the object but might be calculated on request. Either way the pangenome
#' matrix is a function of the gene grouping and should be changed by changing 
#' the gene grouping instead of being manipulated downstream.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A matrix with organisms as columns and gene groups as rows
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' head(pgMatrix(testPG))
#' 
#' @export
#' 
setGeneric('pgMatrix', def = function(object) {
    standardGeneric('pgMatrix')
})
#' Use igraph to create gene grouping from a similarity matrix
#' 
#' This method takes a similarity matrix based on all genes in the pangenome,
#' converts it to a graph representation and uses one of igraphs community 
#' detection algorithms to split all genes into groups. Within the FindMyFriends
#' framework the similarity matrix would usually come from 
#' \code{\link{kmerSimilarity}}, but it can just as well be defined in other 
#' ways e.g. be blast derived.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... parameters to be passed on to the community detection algorithm
#' 
#' @return An object of the same class as 'object'.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Too heavy to include
#' \dontrun{
#' # Generate similarity matrix
#' simMat <- kmerSimilarity(testPG, lowerLimit=0.75)
#' 
#' # Group genes
#' testPG <- graphGrouping(testPG, simMat)
#' }
#' 
#' @family grouping algorithms
#' 
#' @export
#' 
setGeneric('graphGrouping', def = function(object, ...) {
    standardGeneric('graphGrouping')
})
#' Guided Pairwise Comparison grouping of genes
#' 
#' This algorithm recursively builds up a pangenome by merging subpangenomes. 
#' The recursion follows either a supplied hierarchical clustering or one 
#' created using kmer comparison for the full organism. At each step a 
#' representative for each gene group is selected randomly as a representative 
#' and gets compared to all other representatives. Gene groups are then merged 
#' based on the pangenome created for the representatives. Due to the sampling 
#' of representatives at each step there is a certain randomness to the 
#' algorithm. Results should be fairly stable though, as gene groups are 
#' compared multiple times.
#' 
#' @param object A pgVirtual subclass
#'   
#' @param ... parameters passed on.
#'   
#' @return An object of the same class as 'object'.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Too heavy to include
#' \dontrun{
#' testPG <- gpcGrouping(testPG)
#' }
#' 
#' @family grouping algorithms
#'   
#' @export
#' 
setGeneric('gpcGrouping', def = function(object, ...) {
    standardGeneric('gpcGrouping')
})
#' Gene grouping by preclustering with CD-HIT
#' 
#' This grouping algorithm partly mimicks the approach used by Roary, but
#' instead of using BLAST in the second pass it uses cosine similarity of kmer
#' feature vectors, thus providing an even greater speedup. The algorithm uses
#' the CD-HIT algorithm to precluster highly similar sequences and then groups
#' these clusters by extracting a representative and clustering these using the
#' standard FindMyFriends kmer cosine similarity.
#' 
#' @param object A pgVirtual subclass
#'   
#' @param ... parameters passed on.
#'   
#' @return An object of the same class as 'object'.
#' 
#' @references 
#' Page, A. J., Cummins, C. A., Hunt, M., Wong, V. K., Reuter, S., Holden, M. T. 
#' G., et al. (2015). Roary: rapid large-scale prokaryote pan genome analysis. 
#' \emph{Bioinformatics}, btv421.
#' 
#' Fu, L., Niu, B., Zhu, Z., Wu, S., Li, W. (2012). CD-HIT: 
#' accelerated for clustering the next generation sequencing data. 
#' \emph{Bioinformatics}, \bold{28} (23), 3150--3152.
#' 
#' Li, W. and Godzik, A. (2006) Cd-hit: a fast program for clustering and 
#' comparing large sets of protein or nucleotide sequences. 
#' \emph{Bioinformatics}, \bold{22}, 1658--9.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' testPG <- cdhitGrouping(testPG)
#' 
#' @family grouping algorithms
#'   
#' @export
#' 
setGeneric('cdhitGrouping', def = function(object, ...) {
    standardGeneric('cdhitGrouping')
})
setGeneric('precluster', def = function(object, ...) {
    standardGeneric('precluster')
})
#' Define gene grouping manually
#' 
#' In cases where results from other algorithms are wished to be imported into
#' the FindMyFriends framework, this method ensures that the proper formatting 
#' is done. The grouping can  be defined as an integer vector with an element 
#' for each gene. The value of each element is then used as the gene group 
#' classifier. Alternatively groups can be defined by a list of integer vectors.
#' Each element of the list defines a group and the content of each element 
#' refers to gene indexes.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param groups Either a list or integer vector defining the grouping
#' 
#' @return An object of the same class as 'object'.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Load grouping data
#' groups <- system.file('extdata', 'examplePG', 'groupsWG.txt', 
#'     package='FindMyFriends'
#' )
#' groups <- scan(groups, what=integer(), quiet=TRUE)
#' 
#' # Do the grouping
#' testPG <- manualGrouping(testPG, groups)
#' 
#' @family grouping algorithms
#' 
#' @export
#' 
setGeneric('manualGrouping', def = function(object, groups) {
    standardGeneric('manualGrouping')
})
#' Calculate a similarity matrix based on kmers
#' 
#' This method takes a pangenome and calculate a similarity matrix based on 
#' cosine similarity of kmer feature vectors in an all-vs-all fashion. The 
#' result can subsequently be used to group genes either using 
#' \code{\link{graphGrouping}} or homebrewed grouping scheme. In case of the 
#' latter \code{\link{manualGrouping}} should be used to add the grouping back 
#' to the pangenome.
#' 
#' @param object A pgVirtual subclass
#'   
#' @param ... parameters passed on.
#'  
#' @return A matrix (sparse or normal) with cosine similarity for each gene pair 
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Too heavy to include
#' \dontrun{
#' kmerSim <- kmerSimilarity(testPG, lowerLimit=0.75)
#' }
#' 
#' @export
#' 
setGeneric('kmerSimilarity', def = function(object, ...) {
    standardGeneric('kmerSimilarity')
})
#' Split gene groups by neighborhood synteny
#' 
#' This function evaluates already created gene groups and splits the members 
#' into new groups based on the synteny of the flanking genes and the similarity
#' of the sequences. In general the splitting is based on multiple stages that
#' all gene pairs must pass in order to remain in the same group. First the link
#' between the genes is removed if they are part of the same organism. Then the 
#' synteny of the flanking genes are assessed and if it doesn't passes the 
#' defined threshold the link between the gene pair is removed. Then the kmer
#' similarity of the two sequences are compared and if below a certain threshold
#' the link is removed. Lastly the length of the two sequences are compared and 
#' if below a certain threshold the link is removed. Based on this new graph 
#' cliques are detected and sorted based on the lowest within-clique sequence 
#' similarity and neighborhood synteny. The cliques are then added as new groups
#' if the members are not already members of a new group until all members are
#' part of a new group. This approach ensures that all members of the new 
#' groupings passes certain conditions when compared to all other members of the
#' same group. After the splitting a refinement step is done where gene groups 
#' with high similarity and sharing a neighbor either up- or downstream are 
#' merged together to avoid spurius errors resulting from the initial grouping.
#' 
#' @param object A pgVirtualLoc subclass
#'   
#' @param ... parameters passed on.
#' 
#' @return An object with the same class as object containing the new grouping.
#' 
#' @family group-splitting
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE, withGroups=TRUE)
#' 
#' # Too heavy to run
#' \dontrun{
#' testPG <- neighborhoodSplit(testPG, lowerLimit=0.75)
#' }
#' 
#' @export
#' 
setGeneric('neighborhoodSplit', def = function(object, ...) {
    standardGeneric('neighborhoodSplit')
})
#' Link gene groups by homology
#' 
#' This method allows the user to define a secondary grouping of genes be 
#' linking gene groups based on sequence similarity (paralogues). A 
#' representative for each gene group is used for the calculations and the 
#' similarity is assessed using the kmer based cosine similarity.
#' 
#' @param object A pgVirtual subclass
#'   
#' @param ... parameters passed on to the community detection algorithm.
#' 
#' @return An object with the same class as object with linking between gene 
#' groups.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # No paralogue links
#' hasParalogueLinks(testPG)
#' 
#' # Create the links
#' testPG <- kmerLink(testPG)
#' 
#' @export
#' 
setGeneric('kmerLink', def = function(object, ...) {
    standardGeneric('kmerLink')
})
#' Add new organisms to an existing pangenome
#' 
#' This method allows new genomes to be added to an already processed pangenome,
#' preserving existing grouping and adding new genes to their relevant groups.
#' This makes it possible to gradually grow the pangenome as new sequences 
#' becomes available without redoing the grouping at each time, loosing the gene
#' group metadata.
#' 
#' @param object A pgVirtual subclass to merge the new genomes into
#' 
#' @param newSet An object of the same class as object containing the new 
#' organisms to add. Grouping of the genes contained in this object can already
#' exist, if not it will be done automatically.
#'   
#' @param ... parameters passed on.
#' 
#' @return An object of the same class as object containing the new organisms
#' from newSet and possible new gene groups from genes with no orthologues in 
#' the original pangenome.
#' 
#' @examples 
#' # Get base pangenome
#' pg <- .loadPgExample(geneLoc = TRUE, withGroups = TRUE, 
#'                      withNeighborhoodSplit = TRUE)
#' # Get some additional genomes
#' location <- tempdir()
#' unzip(system.file('extdata', 'Mycoplasma.zip', package = 'FindMyFriends'),
#'       exdir = location)
#' genomeFiles <- list.files(location, full.names = TRUE, pattern = '*.fasta')[6:10]
#' pg2 <- pangenome(genomeFiles, translated = TRUE, geneLocation = 'prodigal')
#' 
#' # Combine the two (too computational heavy to include)
#' \dontrun{
#' pg3 <- addGenomes(pg, pg2, nsParam = list(lowerLimit = 0.8))
#' }
#' 
#' @export
#' 
setGeneric('addGenomes', def = function(object, newSet, ...) {
    standardGeneric('addGenomes')
})
#' Remove genes from a pangenome
#' 
#' This method makes it possible to safely remove genes from a pangenome using a
#' variaty of selection mechanisms depending on the supplied parameters. The 
#' name parameter refers to the gene name, organism refers to either organism 
#' name or index, group refers to either gene group name or index and ind refers
#' to the gene index. See examples for details of the different possibilities.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param name A character vector of names of genes to remove
#' 
#' @param organism Either an integer or character vector of orgnanisms to remove
#' genes from. If neither name nor ind is given all genes in the organisms are 
#' removed.
#' 
#' @param group Either an integer or character vector of gene groups to remove
#' genes from. If ind is not given all genes in the groups are removed.
#' 
#' @param ind Indexes of the selections to remove. If both name, organism and 
#' group is not given, it indexes into the raw gene index, otherwise it indexes 
#' into the element defined by organism or group.
#' 
#' @param ... parameters passed on (currently ignored).
#' 
#' @return An object of the same class as object without the genes that should
#' be removed.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' nGenes(testPG)
#' 
#' # Remove gene number 6
#' removeGene(testPG, ind=5)
#' 
#' # Remove all genes from organism 'AE017244'
#' removeGene(testPG, organism='AE017244')
#' 
#' # Remove first gene in gene group 10
#' removeGene(testPG, group=10, ind=1)
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @export
#' 
setGeneric('removeGene', def = function(object, name, organism, group, ind, 
                                        ...) {
    standardGeneric('removeGene')
})
#' Merge paralogue gene groups into new gene groups
#' 
#' This method allows for merging of paralogue gene groups defined using 
#' \code{\link{kmerLink}} into new, bigger, gene groups.
#' 
#' @param object A pgVirtual subclass
#'   
#' @param ... parameters passed on to metadata collapse function. For 
#' combineInfo='merge' sep specifies the separator - sep='none' collapses 
#' information into list elements instead of strings. For combineInfo='largest'
#' no addition arguments are given.
#' 
#' @return An object of the same class as object with the new grouping.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE, withParalogues=TRUE)
#' 
#' # Number of gene groups before collapse
#' nGeneGroups(testPG)
#' 
#' # Number of gene groups after collapse
#' testPG <- collapseParalogues(testPG, combineInfo='largest')
#' nGeneGroups(testPG)
#' 
#' @export
#' 
setGeneric('collapseParalogues', def = function(object, ...) {
    standardGeneric('collapseParalogues')
})
#' Safely add group info
#' 
#' This method allows for adding of group metadata by specifying the name of the
#' metadata and the gene groups it should be added to. It protects the user from
#' overwriting information that is derived from the data, and ensures the proper
#' formatting. Should be prefered to \code{\link{groupInfo<-}} for all but the
#' simplest cases.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... parameters passed on.
#' 
#' @return An object of the same class as object with the new gene group 
#' information.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Create some info
#' info <- data.frame(nickname=c('Tessie', 'Johnny'), index=c(4, 500))
#' 
#' # Add it to the object
#' testPG <- addGroupInfo(testPG, info=info, key='index')
#' 
#' @family Metadata
#' 
#' @export
#' 
setGeneric('addGroupInfo', def = function(object, ...) {
    standardGeneric('addGroupInfo')
})
#' Safely add organisms info
#' 
#' This method allows for adding of organism metadata by specifying the name of
#' the metadata and the organisms it should be added to. It protects the user 
#' from overwriting information that is derived from the data and ensures proper
#' formatting. Should be prefered to \code{\link{orgInfo<-}} for all but the 
#' simplest cases.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... parameters passed on.
#' 
#' @return An object of the same class as object with the added organism 
#' information.
#' 
#' @examples 
#' testPG <- .loadPgExample()
#' 
#' # Create some information
#' info <- data.frame(location=c('Copenhagen', 'Paris', 'London'), 
#'     name=c('AE017243', 'AP012303', 'AE017244')
#' )
#' 
#' # Add the information
#' testPG <- addOrgInfo(testPG, info=info, key='name')
#' 
#' @family Metadata
#' 
#' @export
#' 
setGeneric('addOrgInfo', def = function(object, ...) {
    standardGeneric('addOrgInfo')
})
#' Calculate statistics about each gene group
#' 
#' This method calculates a range of statistics and positional information
#' about each gene group. The information returned are. Maximum number of genes
#' from the same organism (paralogues), shortest sequence length, longest 
#' sequence length, standard deviation of sequence lengths, index of genes in 
#' group, downstream and upstream gene groups.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... parameters passed on.
#' 
#' @return A list with an element for each gene group, each with the following 
#' elements.
#' \describe{
#'  \item{maxOrg}{The highest number of distinct genes from the same organism 
#'  present in the group. A number above 1 indicate the presence of paralogues.}
#'  \item{minLength}{The length of the shortest sequence in the group.}
#'  \item{maxLength}{The length of the longest sequence in the group.}
#'  \item{sdLength}{The standard deviation of lengths in the group.}
#'  \item{genes}{The index for the genes present in the group.}
#'  \item{backward}{A character vector with gene groups separated by ';' that 
#'  lies downstream of the gene group. The number of gene groups for each gene
#'  is controlled by the flankSize argument. If the contig stops before the 
#'  required number of flanking genes have been reached, NA will be added. 
#'  Downstream is defined in relation to the strand of the contig/chromosome, 
#'  and not the translational direction of the gene in question.}
#'  \item{forward}{As above in the other direction.}
#' }
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' grStats <- groupStat(testPG)
#' 
#' @export
#' 
setGeneric('groupStat', def = function(object, ...) {
    standardGeneric('groupStat')
})
#' Calculate statistics about each organism
#' 
#' This method, much like {code{\link{groupStat}}} calculates different 
#' statistics for each organism in the pangenome. Depending on the parameters
#' the statistics are: number of genes, minimum length of gene, maximum length
#' of gene standard deviation of gene lengths, residue frequency, number of gene
#' groups and number of paralogues.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... parameters passed on.
#' 
#' @return A data.frame with a row per organism, with each statistic in a column
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' orgStats <- orgStat(testPG)
#' 
#' @export
#' 
setGeneric('orgStat', def = function(object, ...) {
    standardGeneric('orgStat')
})
#' Calculate the panchromosome graph
#' 
#' This method creates a graph representation of the panchromosome - The 
#' complete set of gene groups linked together by chromosomal position. Each 
#' vertice in the graph represent a gene group and each edge represent a 
#' positional relation between two gene groups (neighboring each other). 
#' Vertices are annotated with number of genes, organism names and strand while
#' edges are annotated with numer of genes (as weight), and organism names.
#' 
#' @param object A pgVirtualLoc subclass
#' 
#' @param ... parameters passed on
#' 
#' @return An igraph object
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE)
#' 
#' panchromosome <- pcGraph(testPG)
#' 
#' @export
#' 
setGeneric('pcGraph', def = function(object, ...) {
    standardGeneric('pcGraph')
})
#' Detect regions of high variability in the panchromosome
#' 
#' This method analyses the panchromosome and detects regions of local 
#' non-linearity. These regions often corresponds to areas with 
#' insertion/deletions, frameshifts or general high plasticity. It works by 
#' examining each vertice of the panchromosome with an out degree above 2 and
#' detect cycles within the neighborhood of these vertices. Adjacent cycles are
#' then joined together to form bigger groups of high variability.
#' 
#' @param object A pgVirtualLoc subclass
#' 
#' @param ... parameters to pass on
#' 
#' @return A list of variable regions. Each element contains the following
#' elements: 
#' \describe{
#'  \item{type}{Either 'ins/den', 'frameshift', 'hub', 'plastic' or 'end'. 
#'  ins/del are regions where the two outgoing vertices are directly connected.
#'  frameshift are regions where the two outgoing vertices are connected through
#'  two different routes, but not directly. hub are regions with more than two 
#'  outgoing vertices. plastic are regions where the two outgoing vertices are
#'  connected through multiple different paths. end are regions with only one
#'  outgoing vertice.}
#'  \item{members}{The gene groups being part of the region.}
#'  \item{flank}{The outgoing vertices connecting the region to the rest of the
#'  panchromosome.}
#'  \item{connectsTo}{The gene group(s) each flank connects to outside of the 
#'  region}
#'  \item{graph}{The subgraph of the panchromosome representing the region}
#' }
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE)
#' 
#' # Too heavy to include
#' \dontrun{
#' regions <- variableRegions(testPG)
#' 
#' # Have a look at the first region
#' regions[[1]]
#' }
#' 
#' @export
#' 
setGeneric('variableRegions', def = function(object, ...) {
    standardGeneric('variableRegions')
})
#' Extract a graph representation of a gene group neighborhood
#' 
#' This method creates a graph representation of the imidiate neighborhood of
#' a gene group. It is different from creating a subgraph of the panchromosome
#' in that only vertices and edges directly reachable from the gene group is 
#' included. The vertices will be annotated with a centerGroup property 
#' indicating whether or not the node is the queried gene group.
#' 
#' @param object A pgVirtualLoc subclass
#' 
#' @param ... Parameters passed on.
#' 
#' @return An igraph object with gene groups as vertices and positional 
#' connections as edges. The edges is weighted according to the number of genes 
#' sharing the connection. All vertices have a centerGroup attribute, which is 
#' FALSE for all but the center group.
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE)
#' 
#' # Look at the surroundings of group 10
#' neighborhood <- getNeighborhood(testPG, group=10)
#' 
#' @seealso \code{\link{plotNeighborhood}} for nice plotting of the neighborhood
#' 
#' @export
#' 
setGeneric('getNeighborhood', def = function(object, ...) {
    standardGeneric('getNeighborhood')
})
#' Plot (very) basic statistics on the pangenome
#' 
#' This method plots the number of genes in each organism and, if gene groups
#' have been defined, the number of singleton, accessory and core gene groups.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters passed on to color scale.
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Should make a nice little plot
#' plotStat(testPG)
#' 
#' @export
#' 
setGeneric('plotStat', def = function(object, ...) {
    standardGeneric('plotStat')
})
#' Plot the neighborhood of a gene group
#' 
#' This method plots the neighborhood extracted using 
#' \code{\link{getNeighborhood}} in a visually pleasing way. It is mainly a 
#' wrapper around \code{\link[igraph]{plot.igraph}} to ensure the proper 
#' information is visualised.
#' 
#' @param object A pgVirtualLoc subclass
#' 
#' @param ... Parameter passed on to igraph's plot method.
#' 
#' @return Called for the side effect of creating a plot. Invisibly returns an 
#' igraph object with all visual parameters set as node and edge attributes.
#' 
#' @examples 
#' testPG <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE)
#' 
#' # Nice little overview of the neighborhood of gene group 30
#' plotNeighborhood(testPG, 30)
#' 
#' @export
#' 
setGeneric('plotNeighborhood', def = function(object, ...) {
    standardGeneric('plotNeighborhood')
})
#' Plot the similarities of genes within a group
#' 
#' This method plots a gene group with genes as vertices and cosine similarities
#' as weighted edges. Mildly informative at best :-)
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters to be passed on to igraphs plotting method
#' 
#' @return Called for the side effect of creating a plot. Invisibly returns an 
#' igraph object with all visual parameters set as node and edge attributes.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' plotGroup(testPG, 10, lowerLimit=0.25)
#' 
#' @export
#' 
setGeneric('plotGroup', def = function(object, ...) {
    standardGeneric('plotGroup')
})
#' Plot the evolution in gene groups
#' 
#' This method constructs a plot showing how the number of singleton, accessory
#' and core gene groups evolve as the size of the pangenome increases. Different
#' ways of increasing the size of the pangenome is available.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters to be passed on
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Standard type - organisms ordered by their index in the pangenome
#' plotEvolution(testPG, ordering='none')
#' 
#' # Bootstrapped with confidence intervals
#' plotEvolution(testPG, ordering='bootstrap')
#' 
#' @export
#' 
setGeneric('plotEvolution', def = function(object, ...) {
    standardGeneric('plotEvolution')
})
#' Create a heatplot with similarities between all organisms
#' 
#' This method creates a heatplot showing the similarity between all organisms
#' in the pangenome. The similarity can either be derived from the pangenome
#' matrix or from kmer calculations of the genes themselves.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters to be passed on.
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Use kmers
#' plotSimilarity(testPG, type='kmer')
#' 
#' # Use pangenome matrix
#' plotSimilarity(testPG, type='pangenome')
#' 
#' @seealso \code{\link{plotTree}} for a dendrogram plot of the same data.
#' 
#' @export
#' 
setGeneric('plotSimilarity', def = function(object, ...) {
    standardGeneric('plotSimilarity')
})
#' Plot a dendrogram of the organisms in a pangenome
#' 
#' This method plots a dendrogram of the relationship between the organisms in
#' the pangenome. It does not tries to by phylogenetic in any way but merely 
#' shows the relationship in data. As with \code{\link{plotSimilarity}} it can
#' be based on either the pangenome matrix or kmer feature vectors.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters to be passed on.
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' plotTree(testPG, type='pangenome', dist='binary', clust='ward.D2')
#' 
#' # And now in a circle (type defaults to 'pangenome')
#' plotTree(testPG, circular=TRUE, dist='binary', clust='ward.D2')
#' 
#' @seealso \code{\link{plotSimilarity}} for a heatmap plot of the same data.
#' 
#' @export
#' 
setGeneric('plotTree', def = function(object, ...) {
    standardGeneric('plotTree')
})
#' Add gene grouping to pangenome
#' 
#' This is an internal function, not meant to be called directly. For adding
#' gene grouping manually see \code{\link{manualGrouping}}. This method is a
#' requirement for classes inheriting from pgVirtual and is not relevant for
#' everyday users.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Parameters to be passed on
#' 
#' @return An object with the same class as object
#' 
#' @note For internal use only. Required for extending the class system of
#' FindMyFriends
#' 
#' @note Required for subclasses of pgVirtual in order to extend the class
#' system of FindMyFriends
#' 
#' @rdname internalGroupGenes
#' @name internal-groupGenes
#' @aliases groupGenes
#' @keywords internal
#' 
#' @export
#' 
setGeneric('groupGenes', def = function(object, ...) {
    standardGeneric('groupGenes')
})
#' Get gene-to-organism relationship
#' 
#' This method returns the organism membership for each gene in the pangenome as
#' a vector of indices. Element 1 corresponds to gene 1 and the value is the
#' index of the corresponding organism.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer vector with an element for each gene in the pangenome.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Stored sequentially so the first will belong to organism 1
#' head(seqToOrg(testPG))
#' 
#' @note Required for extending the class system of FindMyFriends
#' 
#' @seealso \code{\link{seqToGeneGroup}} for gene-to-genegroup relationship
#' 
#' @export
#' 
setGeneric('seqToOrg', def = function(object) {
    standardGeneric('seqToOrg')
})
#' Get gene-to-genegroup relationship
#' 
#' This method returns the group membership for each gene in the pangenome as
#' a vector of indices. Element 1 corresponds to gene 1 and the value is the
#' index of the corresponding gene group. If gene groups have yet to be defined
#' it returns a vector of length 0.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return An integer vector with an element for each gene in the pangenome.
#' 
#' @examples 
#' testPG <- .loadPgExample(withGroups=TRUE)
#' 
#' # Have a look at what the first six genes belongs to
#' head(seqToGeneGroup(testPG))
#' 
#' @note Required for extending the class system of FindMyFriends
#' 
#' @seealso \code{\link{seqToOrg}} for gene-to-organism relationship
#' 
#' @export
#' 
setGeneric('seqToGeneGroup', def = function(object) {
    standardGeneric('seqToGeneGroup')
})
#' Add metadata to the pangenome
#' 
#' These methods are only for internal use and not relevant for regular users.
#' They are required for subclasses of pgVirtual and allows FindMyFriends to add
#' and change metadata as part of the pipelines.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param name The name of the metadata to set
#' 
#' @param info A vector of metadata
#' 
#' @param key The indexes the metadata in info pertains to
#' 
#' @param ... Parameters to be passed on
#' 
#' @return An object of the same class as object
#' 
#' @note For internal use only. Use \code{\link{addGroupInfo}} and 
#' \code{\link{addOrgInfo}} instead. Required for extending the class system of
#' FindMyFriends.
#' 
#' @rdname internalMetadata
#' @name internal-metadata
#' @aliases setGroupInfo
#' @keywords internal
#' 
#' @export
#' 
setGeneric('setGroupInfo', def = function(object, ...) {
    standardGeneric('setGroupInfo')
})
#' @rdname internalMetadata
#' 
#' @export
#' 
setGeneric('setOrgInfo', def = function(object, ...) {
    standardGeneric('setOrgInfo')
})
#' Merge information from two pangenomes
#' 
#' This method is for internal use only and should not be called directly. Use
#' \code{\link{addGenomes}} instead. It is required for subclasses of pgVirtual.
#' 
#' @param pg1 A pgVirtual subclass
#' 
#' @param pg2 An object of the same class as pg1
#' 
#' @param geneGrouping The grouping of the genees in the merged pangenome.
#' Equivalent to calling seqToGeneGroup on the new object
#' 
#' @param groupInfo The metadata on the gene groups in the merged pangenome. 
#' Equivalent to calling groupInfo on the new object
#' 
#' @param ... Parameters to be passed on
#' 
#' @return An object of the same class as pg1
#' 
#' @note For internal use only. Required for extending the class system of 
#' FindMyFriends.
#' 
#' @rdname internalMergePangenomes
#' @name internal-mergePangenomes
#' @aliases mergePangenomes
#' @keywords internal
#' 
#' @export
#' 
setGeneric('mergePangenomes', def = function(pg1, pg2, ...) {
    standardGeneric('mergePangenomes')
})
#' Split gene groups based on similarity
#' 
#' This function splits up gene groups based on cosine similarity of kmer 
#' feature vectors. It uses hard splitting based on a similarity cutoff where
#' unconnected components constitutes new groups. Unlike 
#' \code{\link{neighborhoodSplit}}, paralogues cannot be forced into separate 
#' groups as information needed for this is not present.
#' 
#' @param object A pgVirtual subclass
#' 
#' @param ... Arguments passed on
#' 
#' @return A new pgVirtual subclass object of the same class as 'object'
#' 
#' @family group-splitting
#' 
#' @examples 
#' # Get a grouped pangenome
#' pg <- .loadPgExample(withGroups = TRUE)
#' 
#' \dontrun{
#' # Split groups by similarity (Too heavy to include)
#' pg <- kmerSplit(pg, lowerLimit = 0.8)
#' }
#' 
#' @export
#' 
setGeneric('kmerSplit', def = function(object, ...) {
    standardGeneric('kmerSplit')
})