#' @include aaa.R
#' @include pgVirtualLoc.R
NULL

#' @describeIn neighborhoodSplit Neighborhood-based gene group splitting for
#' pgVirtualLoc subclasses
#' 
#' @param flankSize The number of flanking genes on each side of the gene to use
#' for comparison.
#' 
#' @param minFlank The lowest number of matching genes on either side to allow
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
#' @importFrom dplyr %>% arrange group_by do mutate ungroup
#' 
setMethod(
    'neighborhoodSplit', 'pgVirtualLoc',
    function(object, flankSize, minFlank, forceParalogues, kmerSize, lowerLimit, maxLengthDif=0.1) {
        .fillDefaults(defaults(object))
        
        pending <- rep(TRUE, nGeneGroups(object))
        containsParalogues <- anyParalogues(object)
        currentGrouping <- seqToGeneGroup(object)
        lastCall <- FALSE
        while(any(pending)) {
            if(!lastCall) {
                hadChanges <- FALSE
                thisRound <- which(pending)
                neighbors <- trailGroups2(thisRound, currentGrouping, object, flankSize)
            }
            for(i in seq_along(thisRound)) {
                if(lastCall || qualifies(thisRound[i], neighbors[[i]], containsParalogues[i], pending)) {
                    genes <- as.integer(names(neighbors[[i]]))
                    group <- list(
                        genes = genes,
                        organism = seqToOrg(object)[genes],
                        down = lapply(neighbors[[i]], `[[`, i='down'),
                        up = lapply(neighbors[[i]], `[[`, i='up')
                    )
                    newGroup <- neighborSplitting(
                        group, 
                        pangenome=object, 
                        kmerSize=kmerSize, 
                        lowerLimit=lowerLimit,
                        minFlank=as.integer(containsParalogues[i]),
                        forceParalogues=forceParalogues,
                        maxLengthDif=maxLengthDif
                    )
                    if(length(newGroup) != 1) {
                        newGroupNames <- seq(max(currentGrouping)+1, length.out=length(newGroup))
                        currentGrouping[unlist(newGroup)] <- rep(newGroupNames, lengths(newGroup))
                    }
                    pending[[thisRound[i]]] <- FALSE
                    hadChanges <- TRUE
                }
            }
            if(!hadChanges) {
                if(lastCall) {
                    stop('No convergence')
                }
                lastCall <- TRUE
            }
        }
        manualGrouping(object, match(currentGrouping, unique(currentGrouping)))
    }
)
# setMethod(
#     'similaritySplit', 'pgVirtual',
#     function(object, forceParalogues, kmerSize, lowerLimit) {
#         if(missing(kmerSize)) kmerSize <- object@.settings$kmerSize
#         if(missing(lowerLimit)) lowerLimit <- object@.settings$lowerLimit
#         
#         geneInfo <- data.frame(
#             1:nGenes(object), 
#             organism=orgNames(object)[seqToOrg(object)], 
#             geneGroup=groupNames(object)[seqToGeneGroup(object)]
#         )
#         
#         geneInfo
#     }
# )

### SPLITTING HELPER FUNCTIONS


#' Convert a sorted vector of groups to a vector of neighbors
#' 
#' This function takes a vector of gene groups, sorted by location, and converts
#' it to a vector of neighbors in either direction.
#' 
#' @param groups A vector of gene group indexes
#' 
#' @param dir Either 'f' for forward or 'b' for backward
#' 
#' @param n The number of neighbors to collect for each gene
#' 
#' @return A character vector with neighbors separated by ';'. If the neighbors
#' 'falls of' the end of the vector they will be substituted by NA
#' 
#' @noRd
#' 
collectNeighbors <- function(groups, dir, n) {
    if(dir=='b') {
        groups <- rev(groups)
    }
    groups <- c(groups[-1], NA_character_)
    res <- groups
    n <- n-1
    while(n) {
        groups <- c(groups[-1], NA_character_)
        res <- paste(res, groups, sep=';')
        n <- n-1
    }
    if(dir=='b') {
        rev(res)
    } else {
        res
    }
}
#' Create a matrix of neighborhood similarity
#' 
#' This function compares the neighborhood and organism membership of a group of
#' genes, converting it to a similarity matrix. Similarities are calculated over
#' a series of step. First if forceParalogues=TRUE gene pairs from the same
#' organism gets their similarity set to 0. Secondly if minFlink>0 the size of 
#' the gene group intersection on both sides is calculated and if below minFlank
#' the similarity is set to 0. Lastly for the remaining gene pairs the inverse
#' jackards distance for the neighborhood sequence for gene groups present in
#' both neighborhoods of the two genes is calculated and seet as the similarity.
#' 
#' @param geneGroup A list with genes and their assiciated up- and downstream
#' neighbors as well as the organism its from
#' 
#' @param minFlank The minimum number of corresponding genes on each side of the
#' gene needed to accept similarity
#' 
#' @param forceParalogues Should genes from the same organism be forced to have
#' zero similarity
#' 
#' @return A triangular matrix with the lower half filled with hamming distances
#' of the sequence of the shared genes.
#' 
#' @noRd
#' 
neighborhoodSimilarity <- function(geneGroup, minFlank=1, forceParalogues=TRUE) {
    backward <- geneGroup$down
    forward <- geneGroup$up
    res <- matrix(0, nrow=length(backward), ncol=length(backward))
    for(i in 1:length(backward)) {
        for(j in i:length(backward)) {
            if(j==i) next
            if(forceParalogues && geneGroup$organism[i] == geneGroup$organism[j]) next
            
            bIntersect <- intersect(backward[[i]], backward[[j]])
            if(length(bIntersect) < minFlank && length(backward[[i]]) > 0 && length(backward[[j]]) > 0) next
            
            fIntersect <- intersect(forward[[i]], forward[[j]])
            if(length(fIntersect) < minFlank && length(forward[[i]]) > 0 && length(forward[[j]]) > 0) next
            
            iVec <- c(bIntersect[order(match(bIntersect,backward[[i]]))],
                      fIntersect[order(match(fIntersect,backward[[i]]))])
            jVec <- c(bIntersect[order(match(bIntersect,backward[[j]]))],
                      fIntersect[order(match(fIntersect,backward[[j]]))])
            
            res[j, i] <- sum(iVec == jVec)
        }
    }
    res
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
#' @param geneGroup A list with genes and their assiciated up- and downstream
#' neighbors as well as the organism its from
#' 
#' @param pangenome The pangenome object. A pgVirtual subclass
#' 
#' @param kmerSize The size of the kmers used for similarities
#' 
#' @param lowerLimit The similarity threshold below which it is set to 0
#' 
#' @param maxLegthDif The maximum deviation in sequence length to allow. 
#' Between 0 and 1 it describes a percentage. Above 1 it describes a fixed 
#' length
#' 
#' @param ... Arguments passed on to \code{neighborhoodSimilarity}
#' 
#' @return A list containing the new groups as integer vectors holding the gene
#' index.
#' 
#' @importFrom kebabs linearKernel getExRep spectrumKernel
#' @importFrom igraph graph.adjacency maximal.cliques induced.subgraph E V
#' 
#' @noRd
#' 
neighborSplitting <- function(geneGroup, pangenome, kmerSize, lowerLimit, maxLengthDif, ...) {
    if(length(geneGroup$genes) == 1) return(list(geneGroup$genes))
    
    nMat <- neighborhoodSimilarity(geneGroup, ...)
    seqs <- genes(pangenome, subset=geneGroup$genes)
    sMat <- as.matrix(linearKernel(
        getExRep(seqs, spectrumKernel(kmerSize)), 
        sparse=F,  # To avoid strange crash on AWS
        diag=F, 
        lowerLimit = lowerLimit))
    if(!is.null(maxLengthDif)) {
        if(maxLengthDif < 1) {
            lMat <- outer(width(seqs), width(seqs), function(a,b) {abs(a-b)/pmax(a,b)}) < maxLengthDif
        } else {
            lMat <- outer(width(seqs), width(seqs), function(a,b) {abs(a-b)}) < maxLengthDif
        }
        sMat[!lMat] <- 0
    }
    dimnames(sMat) <- list(geneGroup$genes, geneGroup$genes)
    dimnames(nMat) <- list(geneGroup$genes, geneGroup$genes)
    nMat <- melt(nMat, varnames = c('from', 'to'), value.name = 'nWeight')
    sMat <- melt(sMat, varnames = c('from', 'to'), value.name = 'sWeight')
    aMat <- merge(nMat, sMat)
    aMat <- aMat[aMat$nWeight != 0 & aMat$sWeight != 0,]
    
    gr <- graph.data.frame(aMat, directed=FALSE, vertices = geneGroup$genes)
    
    chosenCliques <- list()
    
    unconnected <- clusters(gr)
    
    for(i in seq_len(unconnected$no)) {
        members <- which(unconnected$membership == i)
        subgr <- induced.subgraph(gr, members)
        
        while(ecount(subgr) != 0) {
            clique <- extractClique(subgr, length(unique(geneGroup$organism[members])))
            subgr <- subgr - clique
            chosenCliques <- append(chosenCliques, list(as.integer(clique)))
        }
        if(vcount(subgr) != 0) {
            chosenCliques <- append(chosenCliques, as.list(as.integer(V(subgr)$name)))
        }
    }
    
    chosenCliques
}
#' Extract the \emph{best} clique from a graph
#' 
#' This function first reduce the density of the graph by removing the worst 
#' edges until the k-core of the graph equals the maximum possible size of 
#' cliques in the graph. After this reduction all largest cliques are detected
#' and sorted by their minimum sWeight and nWeight edges. The clique with the
#' highest minimum weights are returned.
#' 
#' @param gr An undirected igraph object with edge attributes nWeight and 
#' sWeight
#' 
#' @param nDistinct The number of distinct members in the graph, i.e. the 
#' maximum size of a clique in the graph
#' 
#' @return A character vector with the name of the vertices included in the 
#' clique
#' 
#' @importFrom igraph degree vcount V graph.coreness E largest.cliques get.edgelist edge.attributes
#' 
#' @noRd
#' 
extractClique <- function(gr, nDistinct) {
    if(all(degree(gr) == vcount(gr)-1)) return(V(gr)$name)
    
    if(vcount(gr) > nDistinct && sort(graph.coreness(gr), decreasing=TRUE)[nDistinct] > nDistinct-1) {
        edges <- E(gr)
        edgeOrder <- order(edges$nWeight, edges$sWeight)
        lastUpper <- length(edgeOrder)
        lastLower <- 0
        breakPoint <- lastLower + floor((lastUpper-lastLower)/2)
        subgr <- gr - edges[edgeOrder[seq_len(breakPoint)]]
        coreness <- sort(graph.coreness(subgr), decreasing = TRUE)
        while(coreness[nDistinct] != nDistinct-1 && coreness[nDistinct+1 != nDistinct-1]) { # At least nDistinct must have nDistinct-1 neighbors
            if(coreness[nDistinct] > nDistinct-1) {
                lastLower <- breakPoint
            } else {
                lastUpper <- breakPoint
            }
            breakPoint <- lastLower + floor((lastUpper-lastLower)/2)
            subgr <- gr - edges[edgeOrder[seq_len(breakPoint)]]
            coreness <- sort(graph.coreness(subgr), decreasing = TRUE)
            if(lastLower == lastUpper) break
        }
        gr <- subgr
    }
    cliques <- largest.cliques(gr)
    
    edgelist <- cbind(get.edgelist(gr, names=FALSE), data.frame(edge.attributes(gr)))
    
    cliques[lengths(cliques) == 1] <- NULL
    cliqueStat <- do.call(rbind, lapply(cliques, function(clique) {
        edges <- edgelist[,1] %in% clique & edgelist[,2] %in% clique
        c(min(edgelist$sWeight[edges]), min(edgelist$nWeight[edges]))
    }))
    bestClique <- which(cliqueStat[,1] == max(cliqueStat[,1]))
    if(length(bestClique) != 1) {
        bestClique <- bestClique[which.max(cliqueStat[bestClique, 2])]
    }
    V(gr)$name[cliques[[bestClique]]]
}
#' @importFrom igraph get.edgelist edge.attributes neighborhood induced.subgraph largest.cliques degree vcount V E ecount
#' 
extractClique2 <- function(gr, nDistinct) {
    if(all(degree(gr) == vcount(gr)-1)) return(V(gr)$name)
    
    maxEdges <- ceiling(vcount(gr)/nDistinct)*nDistinct*(nDistinct-1)/2
    currentECount <- ecount(gr)
    if(currentECount > maxEdges) {
        edges <- E(gr)
        badEdges <- order(edges$nWeight, edges$sWeight)[seq_len(currentECount - maxEdges)]
        gr <- gr - E(gr)[badEdges]
    }
    
    edgelist <- cbind(get.edgelist(gr, names=FALSE), data.frame(edge.attributes(gr)))
    
    cliques <- maximal.cliques(gr)
    cliques[lengths(cliques) == 1] <- NULL
    cliqueStat <- do.call(rbind, lapply(cliques, function(clique) {
        edges <- edgelist[,1] %in% clique & edgelist[,2] %in% clique
        c(min(edgelist$sWeight[edges]), min(edgelist$nWeight[edges]))
    }))
    bestClique <- which(cliqueStat[,1] == max(cliqueStat[,1]))
    if(length(bestClique) != 1) {
        bestClique <- bestClique[which.max(cliqueStat[bestClique, 2])]
    }
    V(gr)$name[cliques[[bestClique]]]
}
#' Extract the neighborhood of each gene in multiple gene groups
#' 
#' This function finds the gene group membership of the flanking genes for each
#' gene in one or several gene groups.
#' 
#' @param groups The index of the groups for which the neighborhood of their
#' members must be found
#' 
#' @param currentGrouping The current group membership of the genes (often the 
#' result of seqToGeneGroup(pg))
#' 
#' @param pg A pgVirtual subclass object
#' 
#' @param vicinity The distance from the gene to extract group membership from
#' 
#' @return A list of list of lists. The elements of the outermost list are the 
#' gene groups queried the elements of the middle list are the genes in the 
#' group and the elements of the innermost list are 'up' and 'down' containing
#' the up- and downstream gene groups of each gene respectively.
#' 
#' @importFrom dplyr %>% group_by arrange do
#' 
#' @noRd
#' 
trailGroups2 <- function(groups, currentGrouping, pg, vicinity) {
    genes <- which(currentGrouping %in% groups)
    info <- geneLocation(pg)
    info$gene <- 1:nrow(info)
    info$group <- currentGrouping
    info$organism <- seqToOrg(pg)
    info$location <- paste(info$organism, info$contig, sep='>')
    info <- info[info$location %in% unique(info$location[genes]), , drop=FALSE]
    info <- info %>%
        group_by(location) %>%
        arrange(start, end) %>%
        do(trail={
            geneInd <- which(.$gene %in% genes)
            res <- lapply(geneInd, function(x) {
                trailSeq <- x + c(-1, 1)*vicinity
                if(trailSeq[1] < 1) trailSeq[1] <- 1
                if(trailSeq[2] > nrow(.)) trailSeq[2] <- nrow(.)
                if(.$strand[x] == -1) {
                    trailSeq <- rev(trailSeq)
                }
                downLength <- abs(trailSeq[1] - x)+1
                list(
                    down = .$group[seq(trailSeq[1], x)[-downLength]],
                    up = .$group[seq(x, trailSeq[2])[-1]]
                )
            })
            names(res) <- .$gene[geneInd]
            res
        })
    info <- unlist(info$trail, recursive=FALSE)
    info <- info[match(as.character(genes), names(info))]
    info <- split(info, currentGrouping[genes])
    info[match(as.character(groups), names(info))]
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
    groups <- data.frame(organism=seqToOrg(pangenome), geneGroup=seqToGeneGroup(pangenome))
    groups <- groups %>%
        group_by(geneGroup) %>%
        summarise(paralogues=anyDuplicated(organism)!=0)
    groups$paralogues
}
#' Test if a gene group qualifies for splitting
#' 
#' This function investigates several aspects of the gene group along with the
#' status of its neighboring gene groups to decide whether or not it is fit for
#' splitting.
#' 
#' @param group The gene group - currently unused
#' 
#' @param neighbors The neighborhood of each gene in the group as a list of 
#' lists. Each inner list have an 'up' and 'down' element with the index of the
#' upstream and downstream neighbors respectively.
#' 
#' @param hasParalogues logical. Does the group contain several genes from the 
#' same organism.
#' 
#' @param pending A logical vector with the current status of splitting for all
#' gene groups
#' 
#' @return TRUE if the gene group is fit for splitting and FALSE if not
#' 
#' @noRd
#' 
qualifies <- function(group, neighbors, hasParalogues, pending) {
    if(!hasParalogues) return(TRUE)
    down <- lapply(neighbors, `[[`, i='down')
    down[lengths(down)==0] <- NULL
    downResolved <- all(!pending[unique(unlist(down))], na.rm = TRUE)
    up <- lapply(neighbors, `[[`, i='up')
    up[lengths(up)==0] <- NULL
    upResolved <- all(!pending[unique(unlist(up))], na.rm = TRUE)
    if(all(downResolved, upResolved)) return(TRUE)
    if(length(Reduce(intersect, down)) == 0) return(TRUE)
    if(length(Reduce(intersect, up)) == 0) return(TRUE)
    FALSE
}