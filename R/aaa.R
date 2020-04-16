globalVariables(
    names = c(
        '.'
    )
)
.pkg_variables <- new.env()
#' Load an example pangenome
#' 
#' This function loads an example pangenome at various stages of calculation,
#' useful for examples and tests.
#' 
#' @param lowMem logical. Should the returned object inherit from pgLM
#' 
#' @param geneLoc logical. Should the returned object inherit from pgVirtualLoc
#' 
#' @param withGroups logical. Should gene groups be defined
#' 
#' @param withNeighborhoodSplit logical. Should neighborhoodsplitting have been
#' performed
#' 
#' @param withParalogues logical. Should paralogue linking have been performed
#' 
#' @return A pgVirtual subclass object to the specifications defined
#' 
#' @examples 
#' # Load standard (pgFull)
#' .loadPgExample()
#' 
#' # Use pgLM
#' .loadPgExample(lowMem=TRUE)
#' 
#' # Create with pgVirtualLoc subclass (here pgFullLoc)
#' .loadPgExample(geneLoc=TRUE)
#' 
#' # Create with grouping information
#' .loadPgExample(withGroups=TRUE)
#' 
#' # Create with gene groups split by neighborhood (pgVirtualLoc implied)
#' .loadPgExample(withNeighborhoodSplit=TRUE)
#' 
#' # Create with paralogue links
#' .loadPgExample(withGroups=TRUE, withParalogues=TRUE)
#' 
#' @importFrom utils unzip
#' @export
#' 
#' @rdname loadPgExample
#' 
.loadPgExample <- function(lowMem = FALSE, geneLoc = FALSE, withGroups = FALSE, 
                           withNeighborhoodSplit = FALSE, 
                           withParalogues = FALSE) {
    location <- tempdir()
    unzip(system.file('extdata', 'Mycoplasma.zip', package = 'FindMyFriends'),
          exdir = location)
    genomeFiles <- list.files(location, full.names = TRUE, pattern = '*.fasta')
    args <- list(paths = genomeFiles[1:5], translated = TRUE, lowMem = lowMem)
    if (geneLoc || withNeighborhoodSplit) {
        args$geneLocation <- 'prodigal'
    }
    obj <- do.call(pangenome, args)
    if (!(withGroups || withNeighborhoodSplit)) {
        return(obj) # Nothing more to do
    }
    if (withNeighborhoodSplit) {
        grFile <- 'groupsNS.txt'
    } else {
        grFile <- 'groupsWG.txt'
    }
    grFile <- system.file('extdata', 'examplePG', grFile, 
                          package = 'FindMyFriends')
    obj <- manualGrouping(obj, scan(grFile, what = integer(), quiet = TRUE))
    if (withParalogues) {
        if (withNeighborhoodSplit) {
            pFile <- 'paraNS.txt'
        } else {
            pFile <- 'paraWG.txt'
        }
        pFile <- system.file('extdata', 'examplePG', pFile, 
                             package = 'FindMyFriends')
        groupInfo(obj)$paralogue <- scan(pFile, what = integer(), quiet = TRUE)
    }
    obj
}
#' Assign object defaults to missing values
#' 
#' This function takes care of investigating the enclosing functions arguments
#' and identifying the missing ones. If they are missing and a default is given
#' this value is assigned to the enclosing functions environment
#' 
#' @param def A named list of default values
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' # Should only be called within methods/functions
#' 
#' # This will obviously fail
#' \dontrun{
#'   t <- function(x) {
#'     x+1
#'   }
#'   t()
#' }
#' 
#' # Using .fillDefaults
#' t <- function(x, defs) {
#'   .fillDefaults(defs)
#'   x+1
#' }
#' 
#' # With defaults
#' t(defs=list(x=5))
#' 
#' # Direct setting takes precedence
#' t(x=2, defs=list(x=5))
#' 
#' # Still fails if defs doesn't contain the needed parameter
#' \dontrun{
#'   t(defs=list(y='no no'))
#' }
#' 
#' # Usually defs are derived from the object in a method:
#' \dontrun{
#'   setMethod('fillDefExample', 'pgFull',
#'     function(object, x, y) {
#'       .fillDefaults(defaults(object))
#'       x+y
#'     }
#'   )
#' }
#' 
#' @seealso Set and get pangenome defaults with \code{\link{defaults}}
#' 
#' @export
#' 
#' @rdname fillDefaults
#' 
.fillDefaults <- function(def) {
    frame <- sys.frame(-1)
    args <- as.list(frame)
    for (i in names(args)) {
        if (identical(args[[i]], quote(expr = )) && !is.null(def[[i]])) {
            if (!is.null(def$verbose) && def$verbose) {
                message('Using object default ', i, '=', def[[i]])
            }
            assign(i, def[[i]], envir = sys.frame(-1))
        }
    }
}
#' Parallel version of linearKernel from kebabs
#' 
#' This function splits up the computations of linearKernel into chunks that can
#' be distributed and computed in parallel. Each chunk consists of a square of 
#' the lower triangle.
#' 
#' @param x An ExplicitRepresentation
#' 
#' @param pParam The parallelisation parameters to be passed to bplapply
#' 
#' @param nSplits How many groups should x be split up in. The number of chunks
#' will be (nSplits+1)*nSplits/2
#' 
#' @param diag Should the diagonal be retained
#' 
#' @param lowerLimit Threshold for setting values to 0, thus decreasing the 
#' memory requierement for the resulting sparse matrix
#' 
#' @return A lower triangular sparse matrix with dimension nrow(x)*nrow(x)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom kebabs linearKernel
#' 
#' @noRd
#' 
lkParallel <- function(x, pParam, nSplits, diag = FALSE, lowerLimit = 0) {
    chunks <- getChunks(nrow(x), nSplits)
    
    res <- bplapply(
        1:nrow(chunks$combs), 
        function(i, x, chunks, combs, diag, lowerLimit) {
            comb <- combs[i,]
            if (comb$row == comb$col) {
                interval <- chunks[comb$row, 1]:chunks[comb$row, 2]
                linearKernel(x[interval,], sparse = TRUE, diag = diag, 
                             lowerLimit = lowerLimit)
            } else {
                intervalRow <- chunks[comb$row, 1]:chunks[comb$row, 2]
                intervalCol <- chunks[comb$col, 1]:chunks[comb$col, 2]
                linearKernel(x[intervalRow,], x[intervalCol, ], sparse = TRUE, 
                             lowerLimit = lowerLimit)
            }
        }, 
        x = x, 
        chunks = chunks$chunks, 
        combs = chunks$combs, 
        diag = diag, 
        lowerLimit = lowerLimit, 
        BPPARAM = pParam)
    
    weaveChunks(res, chunks)
}
#' Parallel version of linearKernel with dynamic kmer coputations
#' 
#' This version of lkParallel does not accept a precalculated explicit 
#' representation but calculate it in each chunk to put less strain on memory
#' consumption.
#' 
#' @param pangenome A pgVirtual subclass
#' 
#' @param kmerSize The size of kmer
#' 
#' @param pParam The parallelisation parameters to be passed to bplapply
#' 
#' @param nSplits How many groups should x be split up in. The number of chunks
#' will be (nSplits+1)*nSplits/2
#' 
#' @param diag Should the diagonal be retained
#' 
#' @param lowerLimit Threshold for setting values to 0, thus decreasing the 
#' memory requierement for the resulting sparse matrix
#' 
#' @return A lower triangular sparse matrix with dimension nrow(x)*nrow(x)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' 
#' @noRd
#' 
lkParallelLM <- function(pangenome, kmerSize, pParam, nSplits, diag = FALSE, 
                         lowerLimit = 0) {
    chunks <- getChunks(nGenes(pangenome), nSplits)
    
    res <- bplapply(
        1:nrow(chunks$combs), 
        function(i, pangenome, kmerSize, chunks, combs, diag, lowerLimit) {
            comb <- combs[i,]
            if (comb$row == comb$col) {
                interval <- chunks[comb$row, 1]:chunks[comb$row, 2]
                x <- getExRep(genes(pangenome, subset = interval), 
                              spectrumKernel(kmerSize))
                linearKernel(x, sparse = TRUE, diag = diag, 
                             lowerLimit = lowerLimit)
            } else {
                intervalRow <- chunks[comb$row, 1]:chunks[comb$row, 2]
                intervalCol <- chunks[comb$col, 1]:chunks[comb$col, 2]
                x <- getExRep(genes(pangenome, 
                                    subset = c(intervalRow, intervalCol)), 
                              spectrumKernel(kmerSize))
                linearKernel(x[seq_along(intervalRow),], 
                             x[seq_along(intervalCol) + length(intervalRow), ], 
                             sparse = TRUE, 
                             lowerLimit = lowerLimit)
            }
        }, 
        pangenome = pangenome, 
        kmerSize = kmerSize, 
        chunks = chunks$chunks, 
        combs = chunks$combs, 
        diag = diag, 
        lowerLimit = lowerLimit, 
        BPPARAM = pParam)
    
    weaveChunks(res, chunks)
}
#' Split a set of a certain size into chunks of a certain size
#' 
#' This function calculates the number of chunks based on the size of the set 
#' and the size of the chunks. Furthermore it creates a dataframe with chunk 
#' combinations to ensure that all chunks are compared with each others
#' 
#' @param size The number of elements to split
#' 
#' @param nSplits The number of chunks to create
#' 
#' @return A list with two elements: 'combs' contains all combinations of chunks 
#' and 'chunks' contain the start end end indices of the elements in each chunk.
#' 
#' @importFrom utils combn
#' 
#' @noRd
#' 
getChunks <- function(size, nSplits) {
    if (nSplits == 1) {
        return(list(
            combs = data.frame(
                col = 1,
                row = 1
            ),
            chunks = matrix(c(1, size), nrow = 1)
        ))
    }
    chunkSize <- ceiling(size/nSplits)
    from <- seq(from = 1, to = size, by = chunkSize)
    to <- c(from[-1] - 1, size)
    chunks <- cbind(from, to)
    combs <- data.frame(rbind(t(combn(1:nrow(chunks), 2)), 
                              matrix(rep(1:nrow(chunks), 2), ncol = 2)))
    names(combs) <- c('col', 'row')
    combs$origInd <- 1:nrow(combs)
    list(combs = combs, chunks = chunks)
}
#' Combine result from chunk vs chunk computations
#' 
#' This function combines results from doing chunk against chunk computations 
#' into a full matrix - equal to having done all-vs-all on the full set instead
#' of in chunks.
#' 
#' @param squares A list of sparse matrices with results of the chunk vs chunk 
#' computation.
#' 
#' @param split The output of getChunks
#' 
#' @return A sparse matrix with combined results
#' 
#' @importFrom dplyr %>% arrange_ group_by_ do
#' @importFrom Matrix Matrix
#' 
#' @noRd
#' 
weaveChunks <- function(squares, split) {
    oldOptions <- options(dplyr.show_progress = FALSE)
    on.exit({
        options(oldOptions)
    })
    if (length(squares) == 1) {
        return(squares[[1]])
    }
    res <- split$combs %>%
        arrange_(~col) %>%
        group_by_(~col) %>%
        arrange_(~row) %>%
        do(cols = {
            nCol <- ncol(squares[[.$origInd[1]]])
            missingRows <- split$chunks[.$row[1], 1] - 1
            if (missingRows != 0) {
                mat <- list(
                    Matrix(0, ncol = nCol, nrow = missingRows, sparse = TRUE)
                )
            } else {
                mat <- list()
            }
            mat <- append(mat, squares[.$origInd])
            do.call(rbind, mat)
        })
    
    do.call(cbind, res$cols)
}

#' Recursively calculate and merge pangenomes
#' 
#' This function takes care of calculating the pangenome at each branch point
#' of a tree by getting representatives of gene groups from each subtree and 
#' grouping the representatives. The algorithm can cache the result at each
#' branch point so the computations can be taken up if they are interupted.
#' 
#' @param pangenome A pgVirtual subclass
#' 
#' @param tree A dendrogram with a hierarchcal clustering of the genomes in the 
#' pangenome.
#' 
#' @param er An explicit representation of the genes in the pangenome. Either 
#' this or kmerSize must be supplied.
#' 
#' @param kmerSize The size of the kmers to use for calculating the explicit 
#' representation. Either this or er must be supplied.
#' 
#' @param lowerLimit The lower limit of the cosine similarity in order for it to
#' be reported.
#' 
#' @param cacheDB A subclass of filehash to use for result cahcing
#' 
#' @return A list grouping gene indexes in pangenome
#' 
#' @importFrom digest digest
#' @importFrom filehash dbExists dbFetch dbInsert
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom Biostrings order
#' @importFrom stats is.leaf runif
#' 
#' @noRd
#' 
recurseCompare <- function(pangenome, tree, er, clusters, kmerSize, lowerLimit, 
                           cacheDB, pParam) {
    args <- mget(ls())
    args$tree <- NULL
    
    if ('pgGroups' %in% names(attributes(tree))) {
        return(attr(tree, 'pgGroups'))
    }
    if (!missing(cacheDB)) {
        key <- digest(tree)
        if (dbExists(cacheDB, key)) {
            return(dbFetch(cacheDB, key))
        }
    }
    if (is.leaf(tree)) {
        lab <- attr(tree, 'label')
        org <- which(orgNames(pangenome) == lab)
        if (length(org) == 0) {
            org <- suppressWarnings(as.integer(lab))
            if (org > length(pangenome) || is.na(org)) {
                stop('Invalid leaf label - no match')
            }
        } else if (length(org) > 1) {
            stop('Invalid leaf label - multiple matches')
        }
        groups <- as.list(which(seqToOrg(pangenome) == org))
    } else {
        group1 <- do.call(recurseCompare, append(args, list(tree = tree[[1]])))
        group2 <- do.call(recurseCompare, append(args, list(tree = tree[[2]])))
        groups <- c(group1, group2)
    }
    groups <- mergeGroups(groups, clusters)
    represent <- sapply(groups, function(x) {
        x[sample.int(length(x), size = 1L)]
    })
    
    if (length(unlist(groups)) > 1e6 || runif(1) < 0.01) 
        gc() # Ensures always gc() when ngenes reaches 1mill
    
    gen <- genes(pangenome, subset = represent)
    if (missing(er)) {
        erRep <- getExRep(gen,  spectrumKernel(kmerSize))
    } else {
        erRep <- er[represent,]
    }
    members <- lkFMF(erRep, order = order(gen), lowerLimit = lowerLimit, 
                 upperLimit = lowerLimit)
    rm(gen)
    newGroups <- lapply(split(groups, members), unlist)
    if (!missing(cacheDB)) {
        dbInsert(cacheDB, key, newGroups)
    }
    newGroups
}
#' @importFrom igraph make_undirected_graph
mergeGroups <- function(groups, clusters) {
    if (is.na(clusters[1])) {
        return(groups)
    }
    inds <- unlist(groups)
    if (length(unique(clusters[inds])) == length(inds)) {
        return(groups)
    }
    clusters <- rbind(rep(seq_along(groups), lengths(groups)),
                      clusters[unlist(groups)] + length(groups))
    edges <- as.integer(clusters)
    gr <- make_undirected_graph(edges)
    groupClusters <- components(gr)$membership[seq_along(groups)]
    lapply(split(groups, groupClusters), unlist)
}
#' Parallel version of recurseCompare
#' 
#' This function is a parallel version of recurseCompare. It works by splitting 
#' the tree into subtrees and process these in parallel. This will not give a 
#' perfect paralellization as the tree is seldom symmetrical. Furthermore the
#' paralellization can only be used for the lower part of the tree.
#' 
#' @param pangenome A pgVirtual subclass
#' 
#' @param tree A dendrogram with a hierarchcal clustering of the genomes in the 
#' pangenome.
#' 
#' @param er An explicit representation of the genes in the pangenome. Either 
#' this or kmerSize must be supplied.
#' 
#' @param kmerSize The size of the kmers to use for calculating the explicit 
#' representation. Either this or er must be supplied.
#' 
#' @param lowerLimit The lower limit of the cosine similarity in order for it to
#' be reported.
#' 
#' @param pParam A BiocParallelParam subclass
#' 
#' @param Approximate number of subtrees to create
#' 
#' @param cacheDB A subclass of filehash to use for result cahcing
#' 
#' @return A list grouping gene indexes in pangenome
#' 
#' @importFrom BiocParallel bplapply
#' 
#' @noRd
#' 
recurseCompPar <- function(pangenome, tree, er, kmerSize, lowerLimit, pParam, 
                           nSplits, cacheDB) {
    args <- list(pangenome = pangenome, lowerLimit = lowerLimit)
    if (missing(er)) {
        args$kmerSize <- kmerSize
    } else {
        args$er <- er
    }
    if (!missing(cacheDB)) args$cacheDB <- cacheDB
    
    nTrees <- nSplits
    if (attributes(tree)$members/4 > nTrees) {
        nTrees <- floor(attributes(tree)$members/4)
    }
    
    while (nTrees > 4) {
        subtrees <- cutK(tree, nTrees)
        subRes <- bplapply(subtrees$lower, function(subtree, args) {
            args$tree <- subtree
            do.call(recurseCompare, args)
        }, args = args, BPPARAM = pParam)
        tree <- fillTree(subtrees$upper, subRes)
        nTrees <- floor(attributes(tree)$members/4)
    }
    args$tree <- tree
    do.call(recurseCompare, args)
}
#' Add results of subtrees to the branch points of the upper tree
#' 
#' This function adds the result of recurseCompare to the branch point of the
#' subtrees as the pgGroups attribute
#' 
#' @param upper The upper tree that needs to be filled out
#' 
#' @param lowerRes A list with the result of recurseCompare on the subtrees.
#' 
#' @return upper with the leafs filled with lowerRes
#' 
#' @importFrom stats is.leaf
#' 
#' @noRd
#' 
fillTree <- function(upper, lowerRes) {
    if (is.leaf(upper)) {
        leafInd <- as.integer(sub('Branch ', '', attributes(upper)$label))
        attributes(upper)$pgGroups <- lowerRes[[leafInd]]
    } else {
        upper[[1]] <- fillTree(upper[[1]], lowerRes)
        upper[[2]] <- fillTree(upper[[2]], lowerRes)
    }
    return(upper)
}

#' Cluster organisms in pangenome
#' 
#' This function clusters the genomes in a pangenome, either based on the 
#' pangenome matrix or the cosine similarity of the combined genome (all genes
#' concatenated).
#' 
#' @param pangenome A subclass of pgVirtual
#' 
#' @param type Either 'pangenome' or 'kmer'. Should the clustering be based on 
#' the pangenome matrix or kmer counts of all genes in the genomes
#' 
#' @param kmerSize If type='kmer' the kmerSize to use for kmer based similarity.
#' 
#' @param dist The distance function to use. All possible values of method in 
#' dist() are allowed as well as 'cosine' for type='kmer'
#' 
#' @param clust The clustering function to use. Passed on to method in hclust()
#' 
#' @param chunkSize The number of organisms to process at once in parallel
#' 
#' @param pParam A BiocParallelParam subclass
#' 
#' @return A dendrogram object
#' 
#' @importFrom stats hclust as.dendrogram
#' 
#' @noRd
#' 
orgTree <- function(pangenome, type, kmerSize, dist, clust, chunkSize = 100, 
                    pParam) {
    distances <- switch(
        type,
        pangenome = pgDist(pangenome, method = dist),
        kmer = kmerDist(pangenome, kmerSize, chunkSize = chunkSize, pParam, 
                        method = dist)
    )
    clusters <- hclust(distances, method = clust)
    as.dendrogram(clusters)
}
#' Get cosine similarity of organisms
#' 
#' This function calculate the cosine similarity of organisms based on the total
#' kmer count of all genes in the organism.
#' 
#' @param pangenome A subclass of pgVirtual
#' 
#' @param kmerSize If type='kmer' the kmerSize to use for kmer based similarity.
#' 
#' @param chunkSize The number of organisms to process at once in parallel
#' 
#' @param pParam A BiocParallelParam subclass
#' 
#' @return A matrix with cosine similarities in the lower triangle
#' 
#' @importFrom kebabs linearKernel
#' 
#' @noRd
#' 
kmerSim <- function(pangenome, kmerSize, chunkSize = 100, pParam) {
    er <- orgExRep(pangenome, kmerSize, chunkSize = chunkSize, pParam)
    sim <- linearKernel(er, sparse = FALSE)@.Data
    diag(sim) <- 1
    sim[sim < 0] <- 0
    sim
}
#' Get the distance between organisms based on kmers
#' 
#' This function calculates the distance between all organisms based on their
#' kmer feature vector. If method='cosine' The cosine similarity is calculated 
#' and the distance is defined as sqrt(1-cosSim). Otherwise the distance is 
#' calculated directly on the kmer feature vector using dist().
#' 
#' @param pangenome A subclass of pgVirtual
#' 
#' @param kmerSize If type='kmer' the kmerSize to use for kmer based similarity.
#' 
#' @param chunkSize The number of organisms to process at once in parallel
#' 
#' @param pParam A BiocParallelParam subclass
#' 
#' @param method The method to use for distance calculation. Either 'cosine' or
#' a valid value for the method parameter of dist().
#' 
#' @return A distance matrix
#' 
#' @importFrom stats as.dist dist
#' 
#' @noRd
#' 
kmerDist <- function(pangenome, kmerSize, chunkSize = 100, pParam, 
                     method = 'cosine') {
    if (method == 'cosine') {
        sim <- kmerSim(pangenome, kmerSize, chunkSize = chunkSize, pParam)
        as.dist(sqrt(1 - sim))
    } else {
        er <- orgExRep(pangenome, kmerSize, chunkSize = chunkSize, pParam)
        dist(er@.Data, method = method)
    }
}
#' Calculate explicit representations of genomes
#' 
#' This function concatenates all genes in each genome and calculate the 
#' explicit representation of the result. The '-' separator is used to ensure 
#' that kmers do not 'walk over' joined genes. Because of this, this approach
#' should be equal to calculating the explicit representation of each gene and
#' summing up.
#' 
#' @param pangenome A subclass of pgVirtual
#' 
#' @param kmerSize If type='kmer' the kmerSize to use for kmer based similarity.
#' 
#' @param chunkSize The number of organisms to process at once in parallel
#' 
#' @param pParam A BiocParallelParam subclass
#' 
#' @return An ExplicitRepresentationDense object.
#' 
#' @import S4Vectors
#' @importFrom Biostrings AAStringSet DNAStringSet
#' @importFrom kebabs getExRep spectrumKernel
#' 
#' @noRd
#' 
orgExRep <- function(pangenome, kmerSize, chunkSize = 100, pParam) {
    nChunks <- ceiling(length(pangenome)/chunkSize)
    chunks <- rep(1:nChunks, each = chunkSize, length.out = length(pangenome))
    if (missing(pParam)) {
        erList <- lapply(1:nChunks, function(i) {
            genes <- lapply(as.list(genes(pangenome, 'organism', 
                                          subset = which(chunks == i))), 
                            as.character)
            genomes <- unstrsplit(genes, sep = '-')
            if (translated(pangenome)) {
                genomes <- AAStringSet(genomes)
            } else {
                genomes <- DNAStringSet(genomes)
            }
            er <- getExRep(genomes, spectrumKernel(kmerSize), sparse = FALSE)
            er@.Data
        })
    } else {
        erList <- bplapply(1:nChunks, function(i, pangenome, chunks, kmerSize) {
            genes <- lapply(as.list(genes(pangenome, 'organism', 
                                          subset = which(chunks == i))), 
                            as.character)
            genomes <- unstrsplit(genes, sep = '-')
            if (translated(pangenome)) {
                genomes <- AAStringSet(genomes)
            } else {
                genomes <- DNAStringSet(genomes)
            }
            er <- getExRep(genomes, spectrumKernel(kmerSize), sparse = FALSE)
            er@.Data
        }, 
        pangenome = pangenome, 
        chunks = chunks, 
        kmerSize = kmerSize, 
        BPPARAM = pParam)
    }
    kmerMatrix <- rbindMat(erList, fill = 0)
    new('ExplicitRepresentationDense', 
        .Data = kmerMatrix, 
        usedKernel = spectrumKernel(kmerSize), 
        quadratic = FALSE)
}
#' Rbind matrices based on their colnames
#' 
#' This function is much like dplyr's bind_rows, but works with matrices
#' 
#' @param x A matrix or a list of matrices
#' 
#' @param ... Additional matrices if x is a matrix
#' 
#' @param fill The value to fill into non-existing cells
#' 
#' @return A matrix
#' 
#' @noRd
#' 
rbindMat <- function(x, ..., fill = NA) {
    if (inherits(x, 'matrix')) x <- list(x, ...)
    if (length(x) == 1) return(x[[1]])
    cols <- unique(unlist(lapply(x, colnames)))
    rows <- seq(sum(sapply(x, nrow)))
    ans <- matrix(fill, ncol = length(cols), nrow = length(rows), 
                  dimnames = list(NULL, cols))
    currentEnd <- 1
    rowN <- FALSE
    for (i in seq_along(x)) {
        rowInd <- seq(currentEnd, length.out = nrow(x[[i]]))
        ans[rowInd, colnames(x[[i]])] <- x[[i]]
        if (!is.null(rownames(x[[i]]))) {
            if (!rowN) {
                rowN <- TRUE
                rownames(ans) <- 1:nrow(ans)
            }
            rownames(ans)[rowInd] <- rownames(x[[i]])
        }
        currentEnd <- currentEnd + nrow(x[[i]])
    }
    ans
}
#' Extracts location information from prodigal generated files
#' 
#' This function parses prodigal generated sequence names and extract the 
#' contig, start, stop and strand of each protein.
#' 
#' @param desc A character vector of sequence descriptions
#' 
#' @return A data.frame with the columns contig, start, end and strand
#' 
#' @noRd
#' 
prodigalParse <- function(desc) {
    info <- do.call(rbind, strsplit(desc, ' # '))[,-5]
    info[, 1] <- sub('_\\d+$', '', info[,1])
    info <- data.frame(
        contig = info[, 1], 
        start = as.integer(info[,2]), 
        end = as.integer(info[,3]), 
        strand = as.integer(info[,4]), 
        stringsAsFactors = FALSE)
    info
}
#' Extracts and check validity of location information
#' 
#' This function is used to handle the multiple ways gene location can be 
#' supplied in. Furthermore it checks the length and columns of the output.
#' 
#' @param format A data.frame, a list of data.frames, a function taking gene
#' names as the sole input, or 'prodigal'. If it is a data.frame it will be 
#' passed on as-is. If it is a list of data.frames they will be rbinded 
#' together. If it is a function the function will be called with the gene names
#' as input. If it is 'prodigal' the function prodigalParse will be called with
#' the gene names as input. Other named formats (e.g. 'prodigal') might be added
#' in the future.
#' 
#' @param desc The gene names of the genes to get gene location from
#' 
#' @return A data.frame with number of rows equal to length of desc and the 
#' columns 'start', 'end', 'contig' and 'strand'.
#' 
#' @noRd
#' 
getSeqInfo <- function(format, desc) {
    if (is.null(format)) {
        return(data.frame(
            contig = character(), 
            start = integer(), 
            end = integer(), 
            strand = integer())
            )
    }
    
    if (inherits(format, 'data.frame')) {
        seqInfo <- format
    } else if (inherits(format, 'list')) {
        seqInfo <- do.call(rbind, format)
    } else if (inherits(format, 'function')) {
        seqInfo <- do.call(format, list(desc))
    } else if (inherits(format, 'character')) {
        seqInfo <- switch(
            format,
            prodigal = prodigalParse(desc),
            stop('Unknown format: ', format)
        )
    } else {
        stop('Wrong input. Must be NULL, data.frame, list, function or 
             character')
    }
    row.names(seqInfo) <- NULL
    
    if (nrow(seqInfo) != length(desc) || 
        !all(c('contig', 'start', 'end', 'strand') %in% names(seqInfo))) {
        stop('Bad sequenceInfo formatting. Rows must match number of genes and 
             columns must be at least \'contig\', \'start\', \'end\' and 
             \'strand\'')
    }
    seqInfo
}

#' Cut dendrogram into K groups
#' 
#' This function replicates the hclust cutree function with the k parameter for 
#' dendrogram object, which officially only supports cutting at a specified 
#' height.
#' 
#' @param x A dendrogram
#' 
#' @param k The number of groups to cut at
#' 
#' @return A dendrogram
#' 
#' @importFrom stats is.leaf
#' 
#' @noRd
#' 
cutK <- function(x, k) {
    getHeights <- function(x) {
        heights <- attr(x, 'height')
        if (!is.leaf(x)) {
            heights <- c(heights, getHeights(x[[1]]), getHeights(x[[2]]))
        }
        heights
    }
    heights <- sort(getHeights(x), decreasing = TRUE)
    middleHeight <- heights + c(diff(heights), 0)/2
    i <- 1
    while (TRUE) {
        dendro <- cut(x, middleHeight[i])
        if (length(dendro$lower) == k) break
        i <- i + 1
    }
    dendro
}

#' Split an XStringSet by a vector of factors
#' 
#' This function is equivalent to split but works on XStringSet and returns an
#' XStringSetList
#' 
#' @param x An XStringSet
#' 
#' @param f A vector of same length as x with grouping information
#' 
#' @return An XStringSetList subclass equivalent to the XStringSet subclass of x
#' 
#' @noRd
#' 
#' @import IRanges
#' @importClassesFrom Biostrings AAStringSetList DNAStringSetList RNAStringSetList BStringSetList
#' 
splitStringSet <- function(x, f) {
    groups <- split(seq_along(x), f)
    partitioning <- PartitioningByEnd(groups)
    data <- x[unlist(groups)]
    type <- as.character(class(x))
    new(paste0(type, 'List'), 
        unlistData = data, 
        partitioning = partitioning, 
        elementType = type, 
        elementMetadata = NULL, 
        metadata = list())
}
#' Group genes based on grouping of gene groups
#' 
#' This function expands a grouping of gene groups (such as a paralogue link) 
#' into a grouping of the underlying genes
#' 
#' @param groupInd A grouping of genes as returned by seqToGeneGroup()
#' 
#' @param links A grouping of gene groups as a vector of factors, as recorded in
#' groupInfo(object)$paralogue
#' 
#' @return An integer vector with the same length as groupInd giving the 
#' grouping of genes according to group grouping
#' 
#' @noRd
#' 
paralogueInd <- function(groupInd, links) {
    linkInd <- split(seq_along(links), links)
    groupLookup <- rep(seq_along(linkInd), lengths(linkInd))
    groupLookup[match(groupInd, unlist(linkInd))]
}
#' Replacement for the gtable rbind
#' 
#' Used to allow calling of rbind_gtable that can set the size of columns to
#' min and max in addition to first and last.
#' 
#' @param ... gtable objects
#' 
#' @param size A string specifying how the size of each column should be 
#' calculated. If first it is taken from the first gtable, if last, from the 
#' last. If max it is taken as the maximum of all, if min the minimum.
#' 
#' @param z The stack order of the gtables if any
#' 
#' @return A gtable object
#' 
#' @noRd
#' 
rbindGtable <- function(..., size = "max", z = NULL) {
    gtables <- list(...)
    if (!is.null(z)) {
        gtables <- gtable:::z_arrange_gtables(gtables, z)
    }
    Reduce(function(x, y) rbind_gtable(x, y, size = size), gtables)
}
#' Replacement for gtables rbind.gtable
#' 
#' Binds two gtables together, allowing the user to specifiy min and max, in 
#' addition to first and last for column width calculation
#' 
#' @param x A gtable object
#' 
#' @param y A gtable object
#' 
#' @param size A string specifying how the size of each column should be 
#' calculated. If first it is taken from the first gtable, if last, from the 
#' last. If max it is taken as the maximum of all, if min the minimum.
#' 
#' @return A gtable object
#' 
#' @importFrom gtable gtable_add_cols
#' @importFrom grid unit unit.c
#' @noRd
#' 
rbind_gtable <- function(x, y, size = "max") {
    if (nrow(x) == 0) return(y)
    if (nrow(y) == 0) return(x)
    if (ncol(x) > ncol(y)) {
        y <- gtable_add_cols(y, rep(unit(1e-6, 'mm'), ncol(x) - ncol(y)))
        background <- grep('background', y$layout$name)
        y$layout$r[background] <- ncol(y)
    }
    if (ncol(x) < ncol(y)) {
        x <- gtable_add_cols(x, rep(unit(1e-6, 'mm'), ncol(y) - ncol(x)))
        background <- grep('background', x$layout$name)
        x$layout$r[background] <- ncol(x)
    }
    x_row <- length(x$heights)
    y_row <- length(y$heights)
    if (x_row == 0) return(y)
    if (y_row == 0) return(x)
    
    lay_x <- unclass(x$layout)
    lay_y <- unclass(y$layout)
    
    x$layout <- data.frame(
        t = c(lay_x$t, lay_y$t + x_row),
        l = c(lay_x$l, lay_y$l),
        b = c(lay_x$b, lay_y$b + x_row),
        r = c(lay_x$r, lay_y$r),
        z = c(lay_x$z, lay_y$z),
        clip = c(lay_x$clip, lay_y$clip),
        name = c(lay_x$name, lay_y$name)
    )
    
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    
    size <- match.arg(size, c("first", "last", "max", "min"))
    x$widths <- switch(size,
                       first = x$widths,
                       last = y$widths,
                       min = gtable:::compare_unit(x$widths, y$widths, base::pmin),
                       max = gtable:::compare_unit(x$widths, y$widths, base::pmax)
    )
    
    x$grobs <- append(x$grobs, y$grobs)
    
    x
}
#' Import annotation from an .annot file
#' 
#' This function makes it easy to import annotation create in Blast2GO or other
#' programs supporting .annot exporting of results.
#' 
#' @param file The .annot file to import
#' 
#' @return A data.frame ready to merge with a pangenome object using 
#' \code{\link{addGroupInfo}} with the \code{key} argument set to 'name'.
#' 
#' @examples 
#' # Get path to file
#' annot <- system.file('extdata', 'examplePG', 'example.annot', package='FindMyFriends')
#' 
#' # Parse the file
#' readAnnot(annot)
#' 
#' @importFrom dplyr %>% mutate_ group_by_ summarise_
#' @importFrom utils read.table
#' @export
#' 
readAnnot <- function(file) {
    data <- read.table(file, header = FALSE, sep = '\t', fill = TRUE, 
                       stringsAsFactors = FALSE)
    names(data) <- c('name', 'annot', 'desc')
    data <- data %>% 
        mutate_(ontology = ~grepl('GO:', annot)) %>%
        group_by_(~name) %>%
        summarise_(description = ~desc[1], GO = ~I(list(annot[ontology])), 
                  EC = ~I(list(annot[!ontology])))
    data <- as.data.frame(data)
    data$GO <- unclass(data$GO)
    data$EC <- unclass(data$EC)
    data
}
#' Transform similarities
#' 
#' This function is used to transform the similarity matrix created by 
#' linearKernel or the parallel variants, based on whether to normalize and do
#' some kind of transformation on them
#' 
#' @param similarity A matrix (sparse or normal) to transform
#' 
#' @param low The lower limit of the new range (upper is fixed to 1)
#' 
#' @param rescale logical. Should the values be normalized to a new range
#' 
#' @param transform A function taking a matrix (sparse or normal), transforms
#' the values in it, and returns a new matrix of the same dimensions.
#' 
#' @return A matrix of the same type and dimensions as similarity
#' 
#' @importFrom Matrix Matrix
#' 
#' @noRd
#' 
transformSim <- function(similarity, low, rescale, transform) {
    if (!inherits(similarity, 'sparseMatrix')) {
        similarity <- Matrix(similarity, sparse = TRUE)
    }
    if (rescale) {
        similarity@x <- (similarity@x - low)/(1 - low)
    }
    if (inherits(transform, 'function')) {
        transform(similarity)
    } else {
        similarity
    }
}
#' Calculate a pangenome matrix
#' 
#' This function calculates a pangenome matrix based on the basic information
#' available in a pgVirtual subclass. Depending on the class implementation,
#' better approaches might be available, but this approach is ensured to be 
#' possible.
#' 
#' @param object A pgVirtual subclass
#' 
#' @return A dgCMatrix with organisms in columns and gene groups in rows
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @noRd
#' 
getPgMatrix <- function(object) {
    mat <- createPanMatrix(seqToOrg(object), seqToGeneGroup(object))
    sparseMatrix(i = mat$i, p = mat$p, x = mat$x, 
                 dims = c(nGeneGroups(object), nOrganisms(object)),
                 dimnames = list(groupNames(object), orgNames(object)))
}
#' Convert between list and vector type indexing
#' 
#' This function converts back and forth between storing grouping as a list of
#' integer vectors with members of each groups and one integer vector with
#' group membership for each element.
#' 
#' @param groups Integer vector or list of integer vector
#' 
#' @return Depending on the input. The reverse of the input
#' 
#' @noRd
#' 
convertGrouping <- function(groups) {
    if (is.list(groups)) {
        members <- rep(seq_along(groups), lengths(groups))
        members[unlist(groups)] <- as.integer(members)
    } else {
        members <- split(seq_along(groups), groups)
    }
    members
}
formatSeconds <- function(sec) {
    names <- c('day', 'hour', 'minute', 'second')
    time <- c(0, 0, 0, 0)
    if (sec >= 86400) {
        time[1] <- floor(sec / 86400)
        sec <- sec %% 86400
    }
    if (sec >= 3600) {
        time[2] <- floor(sec / 3600)
        sec <- sec %% 3600
    }
    if (sec >= 60) {
        time[3] <- floor(sec / 60)
        sec <- sec %% 60
    }
    time[4] <- round(sec, 3)
    names <- paste0(names, ifelse(time == 1, '', 's'))
    names <- names[time != 0]
    time <- time[time != 0]
    combTime <- paste(time, names)
    if (length(combTime) > 1) {
        paste(paste(combTime[-length(combTime)], collapse = ', '), 
              tail(combTime, 1), 
              sep = ' and ')
    } else {
        combTime
    }
}
