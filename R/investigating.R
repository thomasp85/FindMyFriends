#' @include aaa.R
#' @include pgVirtual.R
#' @include pgVirtualLoc.R
NULL

#' @describeIn groupStat Group statistics for all pgVirtual subclasses
#' 
#' @param vicinity An integer given the number of flanking gene groups to 
#' traverse
#' 
#' @importFrom dplyr %>% group_by_ arrange_ mutate_ ungroup do
#' 
setMethod(
    'groupStat', 'pgVirtual',
    function(object, vicinity = 1) {
        oldOptions <- options(dplyr.show_progress = FALSE)
        on.exit({
            options(oldOptions)
        })
        if (hasGeneInfo(object)) {
            gInfo <- geneLocation(object)
            gInfo$width <- geneWidth(object)
        } else {
            gInfo <- data.frame(geneWidth(object))
        }
        gInfo$gene <- 1:nGenes(object)
        gInfo$organism <- orgNames(object)[seqToOrg(object)]
        gInfo$geneGroup <- groupNames(object)[seqToGeneGroup(object)]
        
        if (hasGeneInfo(object)) {
            info <- gInfo %>% 
                group_by_(~organism, ~contig) %>%
                arrange_(~start, ~end) %>%
                mutate_(backward = ~collectNeighbors(geneGroup, 'b', vicinity), 
                       forward = ~collectNeighbors(geneGroup, 'f', vicinity)) %>%
                ungroup() %>%
                group_by_(~geneGroup)
        } else {
            info <- gInfo %>%
                group_by_(~geneGroup)
        }
        info <- info %>%
            do(groupInfo = {
                list(
                    maxOrg = max(table(.$organism)),
                    minLength = min(.$width),
                    maxLength = max(.$width),
                    sdLength = sd(.$width),
                    genes = .$gene,
                    backward = if (is.null(.$backward)) NA else .$backward, 
                    forward = if (is.null(.$forward)) NA else .$forward
                )
            }) %>%
            arrange_(~match(geneGroup, groupNames(object)))
        info$groupInfo
    }
)
#' @describeIn orgStat Organism statistics for all pgVirtual subclasses
#' 
#' @param subset Name or indexes of organisms to include
#' 
#' @param getFrequency logical. Should amino/nucleic acid frequency be 
#' calculated
#' 
#' @importFrom dplyr bind_rows
#' @importFrom Biostrings alphabetFrequency
#' @importFrom Matrix colSums
#' 
setMethod(
    'orgStat', 'pgVirtual',
    function(object, subset, getFrequency = FALSE) {
        if (missing(subset)) {
            subset <- seq_along(object)
        } else if (inherits(subset, 'character')) {
            subset <- match(subset, orgNames(object))
        }
        orgInd <- findIn(as.integer(subset), seqToOrg(object))
        orgs <- split(orgInd, seqToOrg(object)[orgInd])
        info <- lapply(orgs, function(org) {
            stat <- data.frame(
                nGenes = length(org),
                minLength = min(geneWidth(object)[org]),
                maxLength = max(geneWidth(object)[org]),
                sdLength = sd(geneWidth(object)[org]),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            if (getFrequency) {
                cbind(stat, 
                      t(colSums(alphabetFrequency(genes(object, subset = org))))
                )
            } else {
                stat
            }
        })
        info <- as.data.frame(bind_rows(info))
        
        rownames(info) <- orgNames(object)[as.integer(names(orgs))]
        info$nGeneGroups <- colSums(
            pgMatrix(object)[, rownames(info), drop = FALSE] != 0)
        if (hasParalogueLinks(object)) {
            links <- split(1:nGeneGroups(object), groupInfo(object)$paralogue)
            parMat <- apply(
                pgMatrix(object)[, rownames(info), drop = FALSE], 
                2, 
                function(x) {
                    sapply(links, function(i) sum(x[i]))
                }
            )
        } else {
            parMat <- pgMatrix(object)[, rownames(info), drop = FALSE]
        }
        info$nParalogues <- colSums(parMat > 1)
        info
    }
)
#' @describeIn pcGraph Panchromosome creation for all pgVirtualLoc subclasses
#' 
#' @param slim Should the returned graph be stripped of all metadata and only
#' capture gene group connectivity. Defaults to FALSE
#' 
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table data.table
#' 
setMethod(
    'pcGraph', 'pgVirtualLoc',
    function(object, slim = FALSE) {
        neighbors <- getNeighbors(object, zeroInd = FALSE)
        geneGroups <- seqToGeneGroup(object)
        gGroupNames <- groupNames(object)
        
        neighbors[
            ,
            c('from', 'to') := .(id, up)
        ][
            to == 0,
            to := NA_integer_
        ][
            ,
            c('from', 'to') := .(geneGroups[from], geneGroups[to])
        ]
        
        edges <- neighbors[
            !is.na(to), 
            .(from, to, id)
        ]
        edges[
            from > to,
            c('from', 'to') := .(to, from)
        ][
            ,
            c('from', 'to') := .(gGroupNames[from], gGroupNames[to])
        ]
        if (slim) {
            return(graph_from_data_frame(unique(edges[, .(from, to)]), FALSE))
        }
        edges <- edges[
            ,
            'organisms' := seqToOrg(object)[id]
        ][
            ,
            .(weight = .N, organisms = list(unique(organisms))),
            by = .(from, to)
        ]
        
        neighbors <- neighbors[
            which(reverse),
            c('down', 'up') := .(up, down)
        ][
            down == 0,
            'down' := NA_integer_
        ][
            up == 0,
            'up' := NA_integer_
        ][
            ,
            .(down = gGroupNames[geneGroups[down]],
              up = gGroupNames[geneGroups[up]],
              name = gGroupNames[geneGroups[id]],
              organisms = orgNames(object)[seqToOrg(object)[id]])
        ][
            ,
            .(organisms = list(organisms),
              upstream = list(up), 
              downstream = list(down)),
            by = name
        ]
        vertices <- data.table(
            name = groupNames(object), 
            groupInfo(object)
        )
        vertices <- neighbors[vertices, on = 'name']
        
        graph_from_data_frame(edges, FALSE, vertices)
    }
)
globalVariables(
    c('id', 'up', 'to', 'from', '.N', 'organisms', 'reverse', 'name')
)
#' @describeIn variableRegions Variable region detection for all pgVirtualLoc
#' subclasses
#' 
#' @param flankSize The size of the neighborhood around vertices with outdegree
#' above 2 in where to search for cycles
#' 
setMethod(
    'variableRegions', 'pgVirtualLoc',
    function(object, flankSize) {
        .fillDefaults(defaults(object))
        
        graph <- pcGraph(object, slim = TRUE)
        regions <- locateCycles(graph, flankSize)
        regions <- mergeCycles(regions)
        summarizeCycles(regions, graph)
    }
)
# setMethod(
#     'stableThreads', 'pgVirtualLoc',
#     function(object, flankSize) {
#         graph <- getGraph(object)
#         regions <- locateCycles(graph, flankSize)
#         regions <- mergeCycles(regions)
#         regions <- summarizeCycles(regions, graph)
#         tempGraph <- graph
#         for(i in regions)
#     }
# )

#' @describeIn getNeighborhood Gene group neighborhoods for all pgVirtualLoc
#' subclasses
#' 
#' @param group Either the name or the index of the group whose neighborhood is
#' of interest
#' 
#' @param vicinity An integer giving the number of gene groups in both 
#' directions to collect
#' 
#' @importFrom igraph V V<-
#' 
setMethod(
    'getNeighborhood', 'pgVirtualLoc',
    function(object, group, vicinity = 4) {
        if (inherits(group, 'character')) {
            group <- which(groupNames(object) == group)
            if (length(group) != 1) {
                stop('Bad group name')
            }
        }
        groupGenes <- findIn(as.integer(group), seqToGeneGroup(object))
        neighborhoods <- trailGroups(groupGenes, pg = object, 
                                     vicinity = vicinity)
        neighborhoods <- lapply(neighborhoods, 
                                function(x) groupNames(object)[x])
        graph <- trailsToGraph(neighborhoods)
        groupName <- groupNames(object)[group]
        V(graph)$centerGroup <- FALSE
        V(graph)$centerGroup[V(graph)$name == groupName] <- TRUE
        graph
    }
)
#' @describeIn plotNeighborhood Gene group neighborhood plotting for all
#' pgVirtualLoc subclasses
#' 
#' @param group The name or index of a group.
#' 
#' @param vicinity An integer giving the number of gene groups in both 
#' directions to collect.
#' 
#' @importFrom igraph V E V<- E<- plot.igraph
#' @importFrom grDevices rgb
#' 
setMethod(
    'plotNeighborhood', 'pgVirtualLoc',
    function(object, group, vicinity = 4, ...) {
        gr <- getNeighborhood(object, group, vicinity)
        V(gr)$color <- 'steelblue'
        V(gr)$color[V(gr)$centerGroup] <- 'forestgreen'
        V(gr)$frame.color <- NA
        V(gr)$label.family <- 'sans'
        V(gr)$label.color <- 'black'
        V(gr)$label.cex <- 0.6
        V(gr)$size <- 15
        E(gr)$width <- 5
        weightScale <- scaleRange(E(gr)$weight, 0.75, 0)
        E(gr)$color <-  rgb(weightScale, weightScale, weightScale)
        
        plot.igraph(gr, ...)
        invisible(gr)
    }
)
#' @describeIn plotGroup Gene group similiarity plotting for all pgVirtual
#' subclasses
#' 
#' @param group Name or index of the gene group to plot
#' 
#' @param kmerSize The kmer size to use for similarity calculations
#' 
#' @param lowerLimit The lower threshold for similarity below which it will be
#' set to 0
#' 
#' @param rescale logical. Should the similarity be rescaled between lowerLimit
#' and 1
#' 
#' @param transform A transformation function to apply to the similarities
#' 
#' @import BiocGenerics
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' @importFrom igraph graph_from_adjacency_matrix V E plot.igraph
#' 
setMethod(
    'plotGroup', 'pgVirtual',
    function(object, group, kmerSize, lowerLimit, rescale, transform, ...) {
        .fillDefaults(defaults(object))
        
        groupGenes <- unlist(genes(object, 'group', group))
        sim <- linearKernel(getExRep(groupGenes, spectrumKernel(kmerSize), 
                                     sparse = TRUE), 
                            sparse = TRUE, diag = FALSE)
        sim <- transformSim(sim, lowerLimit, rescale, transform)
        
        gr <- graph_from_adjacency_matrix(sim, 'lower', weighted = TRUE)
        V(gr)$color <- 'steelblue'
        V(gr)$frame.color <- NA
        V(gr)$label.family <- 'sans'
        V(gr)$label.color <- 'black'
        V(gr)$label.cex <- 0.6
        V(gr)$size <- 15
        E(gr)$width <- scaleRange(E(gr)$weight, 1, 10)
        plot.igraph(gr, ...)
        invisible(gr)
    }
)

### HELPER FUNCTIONS

#' Extract neighbors for a set of genes
#' 
#' This function extract the trail of gene groups that leads to and from a set 
#' of genes (usually genes from the same gene group).
#' 
#' @param genes The indexes of the genes to trail
#' 
#' @param pg The pgVirtualLoc subclass to extract from
#' 
#' @param vicinity The distance from the gene in both directions to trail to.
#' 
#' @return A list of integer vectors with gene group indexes for each gene 
#' queried
#' 
#' @importFrom dplyr %>% filter_ group_by_ arrange_ do
#' 
#' @noRd
#' 
trailGroups <- function(genes, pg, vicinity) {
    oldOptions <- options(dplyr.show_progress = FALSE)
    on.exit({
        options(oldOptions)
    })
    info <- geneLocation(pg)
    locations <- unique(paste(seqToOrg(pg)[genes], info$contig[genes], 
                              sep = '>'))
    info$gene <- 1:nrow(info)
    info$group <- seqToGeneGroup(pg)
    info$organism <- seqToOrg(pg)
    info <- info %>%
        filter_(~paste(organism, contig, sep = '>') %in% locations) %>%
        group_by_(~organism, ~contig) %>%
        arrange_(~start, ~end) %>%
        do(trail = {
            geneInd <- which(.$gene %in% genes)
            res <- lapply(geneInd, function(x) {
                trailSeq <- x + c(-1, 1)*vicinity
                if (trailSeq[1] < 1) trailSeq[1] <- 1
                if (trailSeq[2] > nrow(.)) trailSeq[2] <- nrow(.)
                if (.$strand[x] == -1) {
                    rev(.$group[trailSeq[1]:trailSeq[2]])
                } else {
                    .$group[trailSeq[1]:trailSeq[2]]
                }
            })
            names(res) <- .$gene[geneInd]
            res
        })
    info <- unlist(info$trail, recursive = FALSE)
    info[match(names(info), as.character(genes))]
}

#' Convert the output of trailGroups to an igraph object
#' 
#' This function take trails/paths and convert it into a graph representation.
#' 
#' @param trails A list of paths through gene groups.
#' 
#' @return An igraph object
#' 
#' @importFrom dplyr %>% group_by_ summarise_
#' @importFrom igraph graph_from_data_frame
#' 
#' @noRd
#' 
trailsToGraph <- function(trails) {
    trails <- unlist(lapply(trails, function(x) c(x, NA)))
    edges <- data.frame(from = trails[-length(trails)], to = trails[-1])
    edges <- edges[!is.na(edges$from) & !is.na(edges$to),]
    edges <- edges %>%
        group_by_(~from, ~to) %>%
        summarise_(weight = ~length(from))
    graph_from_data_frame(edges)
}

#' Scale a vector of values between a low and a high
#' 
#' This a simple rescaling function that normalizes a vector of values to a 
#' certain range.
#' 
#' @param x A vector of values to normalize
#' 
#' @param low The lower bound of the new range
#' 
#' @param high The higher bound of the new range
#' 
#' @return A vector of the same length as x
#' 
#' @noRd
#' 
scaleRange <- function(x, low, high) {
    if (length(unique(x)) == 1) return(rep(mean(c(low, high)), length(x)))
    x <- (x - min(x))/diff(range(x))
    x*(high - low) + low
}
#' Detect small cycles in a graph
#' 
#' This function detects small cycles in a graph by investigating the 
#' neighborhood of all vertices with a degree above 2. The algorithm uses a 
#' breadth first search to find all putative cycles and then runs through these
#' to check if the center vertice is part of the cycle. For each neighborhood
#' investigated cycles are merged if they share more than one vertex.
#' 
#' @param graph An igraph object
#' 
#' @param maxLength The radius of the neighborhoods to search in
#' 
#' @return A list of vectors with each vector holding the members of a cycle
#' 
#' @importFrom igraph V degree make_ego_graph bfs neighbors as.undirected
#' 
#' @noRd
#' 
locateCycles <- function(graph, maxLength=4) {
    potentialSplits <- V(graph)$name[degree(graph) > 2]
    smallGr <- make_ego_graph(graph, order = maxLength, nodes = potentialSplits)
    cycles <- lapply(seq_along(potentialSplits), function(i) {
        tree <- bfs(smallGr[[i]], potentialSplits[i], father = TRUE)
        endPoints <- tree$order[!tree$order %in% tree$father]
        loops <- endPoints[degree(smallGr[[i]], endPoints) > 1]
        cycles <- list()
        if (length(loops) > 0) {
            for (j in loops) {
                route <- getRoute(j, as.numeric(tree$father))
                wayback <- V(smallGr[[i]])$name[route]
                links <- neighbors(smallGr[[i]], j)
                links <- links[links != tree$father[j]]
                for (k in links) {
                    altRoute <- getRoute(k, as.numeric(tree$father))
                    altWayback <- V(smallGr[[i]])$name[altRoute]
                    if (wayback[2] != altWayback[2]) {  # Check if they return to root by different nodes
                        cycle <- unique(c(wayback, altWayback))
                        if (length(cycles) != 0) {
                            overlap <- sapply(cycles, function(x,y) {
                                length(intersect(x,y)) > 1
                            }, y = cycle)
                            if (any(overlap)) {
                                cycle <- unique(c(unlist(cycles[overlap]), 
                                                  cycle))
                                cycles <- cycles[!overlap]
                            }
                        }
                        cycles <- append(cycles, list(cycle))
                    }
                }
            }
        }
        cycles
    })
    unlist(cycles, recursive = FALSE)
}
#' Merge small cycles into clusters
#' 
#' This function takes the output of locateCycles() and merges cycles if they 
#' share more than one vertex.
#' 
#' @param cycles A list as produced by locateCycles
#' 
#' @return A list of the same format as cycles, but possibly with fewer elements
#' 
#' @importFrom igraph components graph_from_adjacency_matrix
#' 
#' @noRd
#' 
mergeCycles <- function(cycles) {
    cycles <- unique(lapply(cycles, sort))
    adjMat <- matrix(0, ncol = length(cycles), nrow = length(cycles))
    for (i in 1:(length(cycles) - 1)) {
        for (j in (i + 1):length(cycles)) {
            if (length(intersect(cycles[[i]], cycles[[j]])) > 1) {
                adjMat[j, i] <- 1
            }
        }
    }
    gr <- graph_from_adjacency_matrix(adjMat, 'lower')
    cycleGroups <- components(gr)$membership
    lapply(split(cycles, cycleGroups), function(cycle) {unique(unlist(cycle))})
}
#' Gather statistics for small cycles
#' 
#' This function takes a list of cycles as produced by locateCycles() or 
#' mergeCycles() and calculates a range of statistics for it.
#' 
#' @param cycles A list as produced by mergeCycles
#' 
#' @param graph An igraph object from where the cycles have been detected
#' 
#' @return A list of the same length as cycles. Each element contains the 
#' following:
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
#' @importFrom igraph V induced_subgraph degree are_adjacent
#' 
#' @noRd
#' 
summarizeCycles <- function(cycles, graph) {
    graphNames <- V(graph)$name
    summary <- lapply(cycles, function(cycle, graph) {
        outsideNeighbors <- lapply(cycle, function(v, gr, cycle) {
            n <- graphNames[neighbors(gr, v)]
            n[!n %in% cycle]
        }, gr = graph, cycle = cycle)
        flank <- cycle[lengths(outsideNeighbors) != 0]
        if (length(flank) == 0) return(NA)
        cycleGraph <- induced_subgraph(graph, cycle)
        cyclic <- all(degree(cycleGraph) <= 2)
        if (length(flank) == 1) {
            type <- 'end'
        } else if (length(flank) > 2) {
            type <- 'hub'
        } else if (are_adjacent(cycleGraph, flank[1], flank[2])) {
            type <- 'ins/del'
        } else if (cyclic) {
            type <- 'frameshift'
        } else {
            type <- 'plastic'
        }
        isFlank <- V(cycleGraph)$name %in% flank
        V(cycleGraph)$flank <- isFlank
        list(
            type = type,
            members = cycle,
            flank = flank,
            connectsTo = outsideNeighbors[lengths(outsideNeighbors) != 0],
            graph = cycleGraph
        )
    }, graph = graph)
    summary[!is.na(summary)]
}
#' Convert breath first search into paths
#' 
#' This function is used to get the path used to reach a vertex in a bfs or dfs 
#' search. The search must have been run with father=TRUE to get the required
#' information.
#' 
#' @param start The vertex to get the path for
#' 
#' @param fathers The father element from the bfs or dfs search
#' 
#' @return A vector with the route used to reach start
#' 
#' @noRd
#' 
getRoute <- function(start, fathers) {
    nextV <- fathers[start]
    while (!is.na(nextV)) {
        start <- c(nextV, start)
        nextV <- fathers[nextV]
    }
    start
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
    if (dir == 'b') {
        groups <- rev(groups)
    }
    groups <- c(groups[-1], NA_character_)
    res <- groups
    n <- n - 1
    while (n) {
        groups <- c(groups[-1], NA_character_)
        res <- paste(res, groups, sep = ';')
        n <- n - 1
    }
    if (dir == 'b') {
        rev(res)
    } else {
        res
    }
}
#' @importFrom dplyr %>% group_by_ arrange_ mutate_ ungroup do transmute_ filter_
orgGraphs <- function(object) {
    oldOptions <- options(dplyr.show_progress = FALSE)
    on.exit({
        options(oldOptions)
    })
    gInfo <- geneLocation(object)
    gInfo$gene <- 1:nGenes(object)
    gInfo$organism <- orgNames(object)[seqToOrg(object)]
    gInfo$geneGroup <- groupNames(object)[seqToGeneGroup(object)]
    
    prep <- gInfo %>% 
        group_by_(~organism, ~contig) %>%
        arrange_(~start, ~end) %>%
        mutate_(up = ~c(gene[-1], NA), 
               reverse = ~gene < up, 
               down = ~c(NA, gene[-n()])) %>%
        ungroup()
    
    vertices <- prep %>%
        group_by_(~organism) %>%
        do(nodes = {data.frame(
            gene = .$gene, 
            geneGroup = .$geneGroup, 
            strand = .$strand, 
            up = ifelse(.$strand == 1, .$up, .$down), 
            down = ifelse(.$strand == 1, .$down, .$up),
            stringsAsFactors = FALSE
        )})
    nodes <- vertices$nodes
    names(nodes) <- vertices$organism
    
    graphs <- prep %>%
        group_by_(~organism) %>%
        transmute_(from = ~ifelse(reverse, gene, up), 
                  to = ~ifelse(reverse, up, gene)) %>%
        filter_(~!is.na(from) & !is.na(to)) %>%
        do(graph = {
            graph_from_data_frame(
                data.frame(from = .$from, to = .$to), 
                FALSE, 
                nodes[[.$organism[1]]]
            )
        })
    
    res <- graphs$graph
    names(res) <- graphs$organism
    res
}
