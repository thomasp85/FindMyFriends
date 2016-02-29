#' @include aaa.R
#' @include pgVirtualLoc.R
#' @include progress.R
NULL

#' @describeIn neighborhoodSplit Neighborhood-based gene group splitting for
#' pgVirtualLoc subclasses
#' 
#' @param flankSize The number of flanking genes on each side of the gene to use
#' for comparison.
#' 
#' @param forceParalogues Force similarity of paralogue genes to 0
#' 
#' @param kmerSize The length of kmers used for sequence similarity
#' 
#' @param lowerLimit The lower limit of sequence similarity below which it will
#' be set to 0
#' 
#' @param maxLengthDif The maximum deviation in sequence length to allow. 
#' Between 0 and 1 it describes a percentage. Above 1 it describes a fixed 
#' length
#' 
#' @param guideGroups An integer vector with prior grouping that, all else being
#' equal, should be prioritized. Used internally.
#' 
#' @importFrom Matrix sparseMatrix
#' @importFrom igraph graph_from_adjacency_matrix components
#' 
setMethod(
    'neighborhoodSplit', 'pgVirtualLoc',
    function(object, flankSize, forceParalogues, kmerSize, lowerLimit, 
             maxLengthDif, guideGroups = NULL) {
        .fillDefaults(defaults(object))
        
        if (is.null(guideGroups)) {
            guideGroups <- rep(1, nGenes(object))
        }
        
        # Presplit by length
        widths <- geneWidth(object)
        grouping <- widthSim(split(seq_len(nGenes(object)), seqToGeneGroup(object)), widths, maxLengthDif, 'Presplitting')
        object <- manualGrouping(object, as.integer(grouping))
        rm(grouping)
        cat('\nPresplitting resulted in ', nGeneGroups(object), ' gene groups\n', sep = '')
        
        gLoc <- getNeighbors(object)
        startGrouping <- seqToGeneGroup(object)
        startGroupingSplit <- split(seq_len(nGenes(object)), startGrouping)
        finalGrouping <- startGrouping
        org <- seqToOrg(object)
        containsParalogues <- groupHasParalogues(startGroupingSplit, org)
        
        prog <- makeProgress(length(startGroupingSplit), 'Splitting   ', 1000, 12)
        easySplits <- lapply(
            which(!containsParalogues), 
            neighborSplitting,
            object = object, seqToOrg = org, neighbors = gLoc, 
            grouping = finalGrouping, widths = widths, 
            maxLengthDif = maxLengthDif, forceParalogues = forceParalogues,
            flankSize = flankSize, kmerSize = kmerSize, lowerLimit = lowerLimit,
            guide = guideGroups,
            prog = prog
        )
        easySplits <- unlist(easySplits, recursive = FALSE)
        finalGrouping[unlist(easySplits)] <- rep(seq_along(easySplits) + max(finalGrouping), lengths(easySplits))
        
        pending <- containsParalogues
        lastCall = 0
        while (any(pending)) {
            if (lastCall < 5) {
                currentRound <- getPotentials(gLoc$down, gLoc$up, pending, gLoc$reverse, startGroupingSplit, startGrouping)
                if (length(currentRound) == 0) {
                    lastCall <- lastCall + 1
                    chosen <- order(lengths(startGroupingSplit[pending]))[seq_len(min(sum(pending), 10))]
                    currentRound <- which(pending)[chosen]
                } else {
                    lastCall <- 0
                }
            } else {
                currentRound <- which(pending)
            }
            
            splits <- lapply(
                currentRound, 
                neighborSplitting,
                object = object, seqToOrg = org, neighbors = gLoc, 
                grouping = finalGrouping, widths = widths, 
                maxLengthDif = maxLengthDif, forceParalogues = forceParalogues,
                flankSize = flankSize, kmerSize = kmerSize, lowerLimit = lowerLimit,
                guide = guideGroups,
                prog = prog
            )
            splits <- unlist(splits, recursive = FALSE)
            finalGrouping[unlist(splits)] <- rep(seq_along(splits) + max(finalGrouping), lengths(splits))
            pending[currentRound] <- FALSE
        }
        object <- manualGrouping(object, split(seq_len(nGenes(object)), finalGrouping))
        cat('\nSplitting resulted in ', nGeneGroups(object), ' gene groups\n', sep = '')
        object <- neighborhoodMerge(object, maxLengthDif)
        cat('\nMerging resulted in ', nGeneGroups(object), ' gene groups\n', sep = '')
        object
    }
)
#' @describeIn kmerSplit Kmer similarity based group splitting for pgVirtual 
#' subclasses
#' 
#' @param kmerSize The length of kmers used for sequence similarity
#' 
#' @param lowerLimit The lower limit of sequence similarity below which it will
#' be set to 0
#' 
#' @param maxLengthDif The maximum deviation in sequence length to allow. 
#' Between 0 and 1 it describes a percentage. Above 1 it describes a fixed 
#' length
#' 
#' @param pParam An optional BiocParallelParam object that defines the workers 
#' used for parallelisation.
#' 
#' @importFrom BiocParallel SerialParam bplapply
#' 
setMethod(
    'kmerSplit', 'pgVirtual',
    function(object, kmerSize, lowerLimit, maxLengthDif, pParam) {
        .fillDefaults(defaults(object))
        
        if (missing(pParam)) pParam <- SerialParam()
        
        groups <- split(1:nGenes(object), seqToGeneGroup(object))
        newGroups <- bplapply(groups, kmerSplitting, pangenome = object, 
                              kmerSize = kmerSize, lowerLimit = lowerLimit, 
                              maxLengthDif = maxLengthDif, BPPARAM = pParam)
        newGroups <- unlist(newGroups, recursive = FALSE)
        manualGrouping(object, newGroups)
    }
)

### SPLITTING HELPER FUNCTIONS

#' Split a group of genes based on kmer similarity
#' 
#' This function splits a set of group into subgroups based on their sequence 
#' similarities and optionally lengths
#' 
#' @param i The indexes of the genes in the group
#' 
#' @param kmerSize The word size used for cosine similarity
#' 
#' @param lowerLimit The lower limit below which sequences are deemed unsimilar
#' 
#' @param maxLegthDif The maximum deviation in sequence length to allow. 
#' Between 0 and 1 it describes a percentage. Above 1 it describes a fixed 
#' length
#' 
#' @return A list of integer vectors with members of the new groups
#' 
#' @importFrom kebabs getExRep linearKernel spectrumKernel
#' @importFrom Biostrings width
#' @importFrom igraph graph_from_adjacency_matrix components V
#' 
#' @noRd
#' 
kmerSplitting <- function(i, pangenome, kmerSize, lowerLimit, maxLengthDif) {
    geneSeqs <- genes(pangenome, subset = i)
    er <- getExRep(geneSeqs, spectrumKernel(kmerSize))
    lkMat <- kebabs::linearKernel(er, sparse = TRUE, diag = FALSE, 
                          lowerLimit = lowerLimit)
    if (!is.null(maxLengthDif)) {
        if (maxLengthDif < 1) {
            lMat <- outer(width(geneSeqs), width(geneSeqs), 
                          function(a, b) {
                              abs(a - b)/pmax(a, b)
                          }) < maxLengthDif
        } else {
            lMat <- outer(width(geneSeqs), width(geneSeqs), 
                          function(a, b) {abs(a - b)}) < maxLengthDif
        }
        lkMat[!lMat] <- 0
    }
    dimnames(lkMat) <- list(i, i)
    
    gr <- graph_from_adjacency_matrix(lkMat, mode = 'lower', diag = FALSE, 
                                      weighted = TRUE)
    
    members <- components(gr)$membership
    
    split(as.integer(names(members)), members)
}
#' Split gene groups based on neighborhood and similarity
#' 
#' This function takes a gene group and splits it into subgroups based on the
#' neighborhood similarity as calculated by neighborhoodSimilarity(), sequence
#' similarity based on kmers and sequence length difference. The splitting is
#' done by creating a graph where gene groups are only connected if they pass
#' all 3 tests (neighborhood, sequence and length). Based on this graphs cliques
#' are extracted, favouring cliques with high internal sequence and neighborhood 
#' similarity, until all genes have been assigned to a group.
#' 
#' @param group The index of the group to split
#' 
#' @param object The pgVirtual subclass object
#' 
#' @param seqToOrg The value returned when calling seqToOrg() on object
#' 
#' @param neighbors A data.frame as returned by calling getNeighbors
#' 
#' @param grouping The current grouping of the genes, similar in format to 
#' calling seqToGeneGroup on object, but possibly modified in advance
#' 
#' @param widths The sequence length of each sequence
#' 
#' @param maxLengthDif The maximal deviation in sequence length allowed when
#' grouping two sequences
#' 
#' @param forceParalogues Logical. Should genes from the same genome be forced
#' apart?
#' 
#' @param flankSize The length to traverse the neighbors to either side
#' 
#' @param kmerSize The size of the kmers to use for sequence comparison
#' 
#' @param lowerLimit The lower threshold in cosine similarity that sequences
#' must have in order to get grouped
#' 
#' @param guideGroups An a priori grouping that will, all else being equal, be
#' prefered when the grouping is performed
#' 
#' @return A list containing the new groups as integer vectors holding the gene
#' index.
#' 
#' @importFrom kebabs linearKernel getExRep spectrumKernel
#' 
#' @noRd
#' 
neighborSplitting <- function(group, object, seqToOrg, neighbors, grouping, 
                              widths, maxLengthDif, forceParalogues, flankSize, 
                              kmerSize, lowerLimit, guideGroups, prog) {
    members <- which(grouping == group)
    if (length(members) == 1) {
        if (!missing(prog)) {
            progress(prog)
        }
        return(list(members))
    }
    nSim <- neighborhoodSim(members - 1, grouping, seqToOrg, flankSize, 
                            neighbors$down, neighbors$up, neighbors$reverse, 
                            widths, maxLengthDif, forceParalogues)
    seqs <- genes(object, subset=members)
    sSim <- lkFMFmat(getExRep(seqs, spectrumKernel(kmerSize)), order = order(seqs),
                  lowerLimit = lowerLimit, upperLimit = 1)
    edges <- mergeSims(nSim$i, nSim$p, nSim$x, sSim@i, sSim@p, sSim@x, 
                       guideGroups)
    if (nrow(edges) == 0) {
        if (!missing(prog)) {
            progress(prog)
        }
        return(split(members, seq_along(members)))
    }
    res <- split(members, extractCliques(edges, length(members)))
    if (!missing(prog)) {
        progress(prog)
    }
    res
}
#' Extract the \emph{best} clique from a graph
#' 
#' This function first reduce the density of the graph by removing the worst 
#' edges until the k-core of the graph equals the maximum possible size of 
#' cliques in the graph. After this reduction all largest cliques are detected
#' and sorted by their minimum sWeight and nWeight edges. The clique with the
#' highest minimum weights are returned.
#' 
#' @param edges A data.frame with the edges in the graph. The start and end of
#' the edgesis recorded in the "to" and "from" columns. The data.frame must
#' additionally contain the columns "nSim", "sSim" and "gSim", with the 
#' similarity of the neighborhood, sequence and guidegroup respectively.
#' 
#' @param nNodes The number of nodes in the graph. If the last node is 
#' unconnected it will not be recorded in the edges data.frame, so this 
#' information is necessary.
#' 
#' @return An integer vector with the membership of each node
#' 
#' @importFrom igraph make_undirected_graph gorder degree V<-
#' 
#' @noRd
#' 
extractCliques <- function(edges, nNodes) {
    edges <- edges[order(edges$nSim, edges$sSim, edges$gSim, decreasing = TRUE),]
    gr <- make_undirected_graph(rbind(edges$from, edges$to), n = nNodes)
    if (all(degree(gr) == nNodes - 1)) {
        return(list(rep(1, nNodes)))
    }
    V(gr)$ID <- seq_len(nNodes) - 1
    getCliques(gr)
}
#' Get the adjacent genes for for each gene in a pangenome
#' 
#' This function is used to extract neighbor information for each gene in a 
#' pangenome. The returned information is 0-indexed as it is mainly used to pass
#' into C++ functions.
#' 
#' @param pg A pgVirtualLoc subclass object
#' 
#' @return A data.frame with column "id", "down", "up" and "reverse". Id gives 
#' the id of the group, down the neighbor downwards and up the neighbor upward.
#' Reverse tells if the gene is located on the reverse strand. If a gene is 
#' located at the beginning or end of a DNA string, -1 will be used to indicate
#' absence of neighbors in up and down
#' 
#' @importFrom dplyr %>% group_by arrange transmute n ungroup arrange
#' 
#' @noRd
#' 
getNeighbors <- function(pg) {
    oldOptions <- options(dplyr.show_progress = FALSE)
    on.exit({
        options(oldOptions)
    })
    gLoc <- geneLocation(pg)
    gLoc$id <- seq_len(nGenes(pg))
    gLoc$org <- seqToOrg(pg)
    
    gLoc <- gLoc %>% group_by(org, contig) %>% 
        arrange(start, end) %>% 
        transmute(id = id-1, 
                  down = c(-1, id[-n()]), 
                  up = c(id[-1], -1), 
                  reverse = strand == -1) %>% 
        ungroup() %>% 
        arrange(id)
    
    as.data.frame(gLoc)[, -(1:2)]
}
#' Determine which gene groups contains paralogues
#' 
#' This function simply investigates whether or not each group contain multiple 
#' genes from the same organism
#' 
#' @param pangenome A pgVirtual subclass object
#' 
#' @return A logical vector indicating if each gene group in the pangenome 
#' contains paralogues
#' 
#' @importFrom dplyr %>% group_by summarise
#' 
#' @noRd
#' 
anyParalogues <- function(pangenome) {
    groups <- data.frame(organism = seqToOrg(pangenome), 
                         geneGroup = seqToGeneGroup(pangenome))
    groups <- groups %>%
        group_by(geneGroup) %>%
        summarise(paralogues = anyDuplicated(organism) != 0)
    groups$paralogues
}
#' @importFrom igraph degree 
neighborhoodMerge <- function(pangenome, maxLengthDif) {
    cdhitOpts <- list(l = 5, n = 5)
    if (maxLengthDif < 1) {
        cdhitOpts$s <- 1 - maxLengthDif
    } else {
        cdhitOpts$S <- maxLengthDif
    }
    cdhitOpts <- lapply(cdhitOpts, as.character)
    first <- TRUE
    while (TRUE) {
        pc <- pcGraph(pangenome)
        knots <- which(degree(pc) > 2)
        if (length(knots) == 0) break
        knots <- match(V(pc)$name[knots], groupNames(pangenome))
        geneInd <- which(seqToGeneGroup(pangenome) %in% knots)
        neighbors <- trailGroups(geneInd, pangenome, 1)
        neighbors <- split(neighbors, seqToGeneGroup(pangenome)[geneInd])
        neighbors <- lapply(seq_along(neighbors), function(i) {
            groupInd <- as.integer(names(neighbors)[i])
            n <- neighbors[[i]]
            nNeighbors <- lengths(n)
            groupPos <- sapply(n, function(nn) which(nn == groupInd))
            down <- groupPos - 1
            up <- groupPos + 1
            down <- unlist(Map(`[`, n[down != 0], down[down != 0]))
            up <- unlist(Map(`[`, n[up <= nNeighbors], up[up <= nNeighbors]))
            list(sort(unique(down)), sort(unique(up)))
        })
        neighbors <- unlist(neighbors, recursive = FALSE)
        neighbors <- unique(neighbors[lengths(neighbors) > 1])
        neighborlookup <- data.frame(
            OG = unlist(neighbors), 
            NG = rep(seq_along(neighbors), lengths(neighbors))
        )
        GOI <- unique(unlist(neighbors))
        repGOI <- getRep(pangenome, 'longest')[GOI]
        if (first) {
            first <- FALSE
        } else {
            cat('\n')
        }
        equals <- cdhit(repGOI, cdhitOpts, 'Merging     ')
        GOI <- split(GOI, equals)
        GOI <- GOI[lengths(GOI) > 1]
        pairs <- lapply(GOI, function(groups) {
            pairs <- combn(groups, 2)
            matched <- sapply(seq_len(ncol(pairs)), function(i) {
                any(neighborlookup$NG[neighborlookup$OG == pairs[1, i]] %in%
                        neighborlookup$NG[neighborlookup$OG == pairs[2, i]])
            })
            pairs[, matched]
        })
        pairs <- do.call(cbind, pairs)
        
        if (ncol(pairs) == 0) break
        
        dupPairs <- apply(matrix(duplicated(as.vector(pairs)), nrow = 2), 2, any)
        pairs <- pairs[, !dupPairs, drop = FALSE]
        currentGroups <- seqToGeneGroup(pangenome)
        toChange <- lapply(seq_len(ncol(pairs)), function(i) {
            which(currentGroups %in% pairs[,i])
        })
        currentGroups[unlist(toChange)] <- rep(seq.int(max(currentGroups) + 1, 
                                                       length.out = length(toChange)),
                                               lengths(toChange))
        pangenome <- manualGrouping(pangenome, split(seq_len(nGenes(pangenome)), currentGroups))
    }
    pangenome
}
#' igraph functions to access through Rcpp
#' 
#' This list allows one to access igraph functions without putting igraph in 
#' Depends, thus polluting the seqrch path
#' 
#' @keywords Internal
#' 
#' @export
#' 
#' @noRd
#' 
.igraphFunctions <- list(
    neighbors = igraph::neighbors,
    gorder = igraph::gorder,
    gsize = igraph::gsize,
    ends = igraph::ends,
    vertex_attr = igraph::vertex_attr,
    delete_vertices = igraph::delete_vertices
)