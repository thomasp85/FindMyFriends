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
        
        neighbors <- getNeighbors(object, flankSize)
        
        newGroups <- lapply(
            neighbors, 
            neighborSplitting, 
            pangenome=object, 
            kmerSize=kmerSize, 
            lowerLimit=lowerLimit,
            minFlank=minFlank,
            forceParalogues=forceParalogues,
            maxLengthDif=maxLengthDif
        )
        manualGrouping(object, unlist(newGroups, recursive=FALSE))
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

#' Extract the neighbors for each gene in a pangenome
#' 
#' This function extract neighborhood information for each gene in terms of the
#' gene groups the flanking genes are members of. The direction of the gene is
#' normalized so that forward and backward is flipped if the gene is located on
#' the complementary strand (strand=-1)
#' 
#' @param pangenome A pgVirtualLoc subclass object
#' 
#' @param flankSize The number of genes in each direction to extract
#' 
#' @return A list with an element for each gene group. Each element contains the
#' gene members as indexes, the backwards and forwards gene groups for each gene
#' member and the organism each gene member comes from.
#' 
#' @importFrom dplyr %>% group_by arrange mutate ungroup do
#' 
#' @noRd
#' 
getNeighbors <- function(pangenome, flankSize) {
    geneLocation <- geneLocation(pangenome)
    geneLocation$gene <- 1:nGenes(pangenome)
    geneLocation$organism <- orgNames(pangenome)[seqToOrg(pangenome)]
    geneLocation$geneGroup <- groupNames(pangenome)[seqToGeneGroup(pangenome)]
    
    neighbors <- geneLocation %>% 
        group_by(organism, contig) %>%
        arrange(start, end) %>%
        mutate(backward=collectNeighbors(geneGroup, 'b', flankSize), forward=collectNeighbors(geneGroup, 'f', flankSize)) %>%
        ungroup() %>%
        group_by(geneGroup) %>%
        do(groupInfo={
            list(
                genes=.$gene, 
                backward=ifelse(.$strand==1,.$backward, .$forward), 
                forward=ifelse(.$strand==1,.$forward, .$backward), 
                organism=.$organism
            )
        }) %>%
        arrange(match(geneGroup, groupNames(pangenome)))
    neighbors$groupInfo
}

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
    backward <- lapply(strsplit(geneGroup$backward, ';'), rev)
    forward <- strsplit(geneGroup$forward, ';')
    res <- matrix(0, nrow=length(backward), ncol=length(backward))
    for(i in 1:length(backward)) {
        for(j in i:length(backward)) {
            if(j==i) next
            if(forceParalogues && geneGroup$organism[i] == geneGroup$organism[j]) next
            
            bIntersect <- intersect(backward[[i]], backward[[j]])
            if(length(bIntersect) < minFlank) next
            
            fIntersect <- intersect(forward[[i]], forward[[j]])
            if(length(fIntersect) < minFlank) next
            
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
        sparse=T, 
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
    cliques <- maximal.cliques(gr)
    minWeights <- do.call(rbind, lapply(cliques, function(i) {
        edges <- E(induced.subgraph(gr, i))
        if(length(edges) != 0) {
            c(min(edges$sWeight), min(edges$nWeight))
        } else {
            c(0, 0)
        }
    }))
    verticeNames <- V(gr)$name
    visitedVertices <- c()
    chosenCliques <- list()
    
    while(TRUE) {
        for(i in order(minWeights[, 1], minWeights[, 2], decreasing = TRUE)) {
            vertices <- verticeNames[cliques[[i]]]
            if(length(intersect(visitedVertices, vertices)) == 0) {
                visitedVertices <- append(visitedVertices, vertices)
                chosenCliques <- append(chosenCliques, list(as.integer(vertices)))
                if(length(visitedVertices) == length(geneGroup$genes)) break
            }
        }
        
        missingVertices <- setdiff(as.character(geneGroup$genes), visitedVertices)
        
        if(length(missingVertices) != 0) {
            gr <- induced.subgraph(gr, missingVertices)
            cliques <- maximal.cliques(gr)
            minWeights <- do.call(rbind, lapply(cliques, function(i) {
                edges <- E(induced.subgraph(gr, i))
                if(length(edges) != 0) {
                    c(min(edges$sWeight), min(edges$nWeight))
                } else {
                    c(0, 0)
                }
            }))
            verticeNames <- V(gr)$name
        } else {
            break
        }
    }
    chosenCliques
}