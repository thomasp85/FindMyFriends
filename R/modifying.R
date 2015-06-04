#' @include generics.R
#' @include aaa.R
#' @include pgVirtual.R

#' @describeIn addGenomes Genome addition for all pgVirtual subclasses
#' 
#' @param kmerSize The size of the kmers to use for comparing new genes to 
#' existing
#' 
#' @param lowerLimit The lower threshold for sequence similarity, below which it
#' is set to 0
#' 
#' @param algorithm The algorithm to use to cluster new genes with existing
#' 
#' @param flankSize The number of genes on each size of a gene to compare when
#' assessing neighborhood similarity
#' 
#' @param gpcParam A list of parameters to pass to \code{\link{gpcGrouping}} in
#' case newSet has not had gene grouping performed.
#' 
#' @param nsParam A list of parameters to pass to 
#' \code{\link{neighborhoodSplit}} in case newSet has not had gene grouping 
#' performed and inherits from pgVirtualLoc
#' 
#' @param klParam A list of parameters to pass to \code{\link{kmerLink}} In case
#' object had paralogues defined
#' 
#' @param pParam A BiocParallelParam object
#' 
#' @param nSplits The number of parallel jobs to split the parallelization into
#' 
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' @importFrom dplyr bind_rows
#' 
#' @noRd
#' 
setMethod(
    'addGenomes', c('pgVirtual', 'pgVirtual'),
    function(object, newSet, kmerSize, lowerLimit, algorithm, flankSize, gpcParam, nsParam, klParam, pParam, nSplits, ...) {
        .fillDefaults(defaults(object))
        
        if(class(object) != class(newSet)) stop('pangenome and newSet must be of the same class')
        if(translated(object) != translated(newSet)) stop('Translational status of pangenome and newSet must match')
        
        if(!hasGeneGroups(newSet)) {
            defaults(newSet) <- defaults(object)
            gpcParam$object <- newSet
            newSet <- do.call(gpcGrouping, gpcParam)
            if(hasGeneInfo(newSet)) {
                nsParam$object <- newSet
                newSet <- do.call(neighborhoodSplit, nsParam)
            }
        }
        
        fullSet <- c(getRep(object, 'random'), getRep(newSet, 'random'))
        er <- getExRep(fullSet, spectrumKernel(kmerSize))
        if(missing(pParam)) {
            sim <- linearKernel(er, sparse=TRUE, diag=FALSE, lowerLimit=lowerLimit)
        } else {
            sim <- lkParallel(er, pParam, nSplits, lowerLimit=lowerLimit)
        }
        similarGroups <- igGroup(sim, algorithm, ...)
        
        members <- switch(
            hasGeneInfo(object),
            'TRUE' = mergeNeighborhood(object, newSet, flankSize, similarGroups, sim),
            'FALSE' = mergeLargest(object, similarGroups)
        )
        
        newGroups <- max(members) - nGeneGroups(object)
        newInfo <- bind_rows(groupInfo(object), data.frame(description=rep(NA, newGroups)))
        newPG <- mergePangenomes(object, newSet, c(seqToGeneGroup(object), members), newInfo)
        newMat <- pgMatrix(newPG)
        newInfo$group <- apply(newMat, 1, function(x) {
            if(all(x!=0)) return('Core')
            if(sum(x!=0) == 1) return('Singleton')
            return('Accessory')
        })
        newInfo$nOrg <- apply(newMat!=0, 1, sum)
        newInfo$nGenes <- apply(newMat, 1, sum)
        groupInfo(newPG) <- newInfo
        if(hasParalogueLinks(object)) {
            klParam$object <- newPG
            newPG <- do.call(kmerLink, klParam)
        }
        newPG
    }
)
#' @describeIn removeGene Remove gene based on gene name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'missing', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- which(geneNames(object) %in% name)
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene and organism naem
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'character', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, name=name, organims=index)
    }
)
#' @describeIn removeGene Remove gene based on gene name and organism index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'numeric', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(n, o) {
            which(n==geneNames(object) & o==seqToOrg(object))
        }, n=name, o=organism))
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)
#' @describeIn removeGene Remove gene based on organism name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'character', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, organism=index)
    }
)
#' @describeIn removeGene Remove gene based on organism index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'numeric', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- which(seqToOrg(object) %in% organism)
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)
#' @describeIn removeGene Remove gene based on organism name and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'character', 'missing', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, organism=index, ind=ind)
    }
)
#' @describeIn removeGene Remove gene based on organism and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'numeric', 'missing', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(o, i) {
            which(seqToOrg(object)==o)[i]
        }, o=organism, i=ind))
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene group name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'character', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(group, groupNames(object))
        removeGene(object, group=index)
    }
)
#' @describeIn removeGene Remove gene based on gene group index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'numeric', 'missing'),
    function(object, name, organism, group, ind) {
        index <- which(seqToGeneGroup(object) %in% group)
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene group name and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'character', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- match(group, groupNames(object))
        removeGene(object, group=index, ind=ind)
    }
)
#' @describeIn removeGene Remove gene based on gene group and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'numeric', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(g, i) {
            which(seqToGeneGroup(object)==g)[i]
        }, g=group, i=ind))
        if(length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind=index)
        }
    }
)

#' @describeIn collapseParalogues Merge paralogue gene groups for all pgVirtual
#' subclasses
#' 
#' @param combineInfo The approach used to combine metadata from the collapsed
#' groups. Either 'merge' for merging, 'largest' for picking information from
#' the largest group, or a function that takes a data.frame of multiple rows and
#' converts it to a data.frame with one row and the same columns.
#' 
#' @importFrom dplyr bind_rows
#' 
setMethod(
    'collapseParalogues', 'pgVirtual',
    function(object, combineInfo='merge', ...) {
        if(!hasParalogueLinks(object)) stop('Paralogues must be detected before collapsing')
        paralogues <- split(1:nGeneGroups(object), groupInfo(object)$paralogue)
        groups <- split(1:nGenes(object), seqToGeneGroup(object))
        newGroups <- lapply(paralogues, function(x) {unlist(groups[x])})
        newInfo <- lapply(paralogues, function(x) {
            if(inherits(combineInfo, 'function')) {
                combineInfo(groupInfo(object)[x,], ...)
            } else {
                switch(
                    combineInfo,
                    merge = mergeInfo(groupInfo(object)[x,], ...),
                    largest = largestInfo(groupInfo(object)[x,], ...),
                    stop('Unknown combine function')
                )
            }
        })
        newInfo <- as.data.frame(bind_rows(newInfo))
        newInfo$paralogue <- NA
        newPG <- manualGrouping(object, newGroups)
        newPG <- addGroupInfo(newPG, newInfo)
        newPG
    }
)
#' @describeIn addGroupInfo Add gene group info safely for all pgVirtual 
#' subclasses
#' 
#' @param info A data.frame with information to add
#' 
#' @param key Either an integer vector with the index of each gene group the 
#' rows in info corresponds to, or the name of the column in info that holds the
#' indexes.
#' 
setMethod(
    'addGroupInfo', 'pgVirtual',
    function(object, info, key) {
        if(missing(key)) {
            key <- 1:nGeneGroups(object)
        } else if(inherits(key, 'character')) {
            keyInfo <- info[[key[1]]]
            info <- info[, names(info) != key[1], drop=FALSE]
            if(is.numeric(keyInfo)) {
                key <- keyInfo
            } else {
                key <- match(keyInfo, groupNames(object))
            }
        }
        if(!is.numeric(key)) stop('key must be convertible to an index')
        for(i in names(info)) {
            if(i %in% c('nGenes', 'nOrg', 'group')) next
            object <- setGroupInfo(object, i, info[[i]], key)
        }
        object
    }
)
#' @describeIn addOrgInfo Add organism info safely for all pgVirtual subclasses
#' 
#' @param info A data.frame with information to add
#' 
#' @param key Either an integer vector with the index of each organism the 
#' rows in info corresponds to, or the name of the column in info that holds the
#' indexes.
#' 
setMethod(
    'addOrgInfo', 'pgVirtual',
    function(object, info, key) {
        if(missing(key)) {
            key <- 1:nOrganisms(object)
        } else if(inherits(key, 'character')) {
            keyInfo <- info[[key[1]]]
            info <- info[, names(info) != key[1], drop=FALSE]
            if(is.numeric(keyInfo)) {
                key <- keyInfo
            } else {
                key <- match(keyInfo, orgNames(object))
            }
        }
        if(!is.numeric(key)) stop('key must be convertible to an index')
        for(i in names(info)) {
            if(i %in% c('nGenes')) next
            object <- setOrgInfo(object, i, info[[i]], key)
        }
        object
    }
)

### HELPER FUNCTIONS
#' Merge content from several rows
#' 
#' This function merges information by either concatenating the unique elements
#' or creating a list element.
#' 
#' @param data a data.frame whose rows should be merged
#' 
#' @param sep Separator to use for concatanation. If omitted rows will be merged 
#' into list elements
#' 
#' @return A data.frame with one row
#' 
#' @noRd
#' 
mergeInfo <- function(data, sep) {
    as.data.frame(lapply(data, function(x) {
        info <- unique(unlist(x))
        if(missing(sep)) {
            list(info)
        } else {
            paste(info, collapse=sep)
        }
    }), stringsAsFactors=FALSE)
}
#' Select the row with the largest group
#' 
#' This function merges the rows of a data.frame by selecting the row 
#' representing the largest gene group. On ties the first one is selected.
#' 
#' @param data A data.frame whose rows should be merged
#' 
#' @return A data.frame with one row
#' 
#' @noRd
#' 
largestInfo <- function(data) {
    largest <- which.max(data$nGenes)
    data[largest, , drop=FALSE]
}
#' Convert row-column pairs into vector indexes
#' 
#' This function takes a range of row-column pairs and returns a vector of 
#' indexes that can be used to directly access all the pairs.
#' 
#' @param row A vector of row indices
#' 
#' @param col A vector of column indices
#' 
#' @param nrow The number of rows in the matrix that needs to be indexed
#' 
#' @return A vector with the index for each row column pair
#' 
#' @noRd
#' 
pairToIndex <- function(row, col, nrow) {
    nrow*(col-1)+row
}
#' Remove a range of indexes
#' 
#' This function is used to reformat a list of indexes where certain elements
#' have been removed. E.g. c(1,1,2,5,4,5,2) is missing 3 so it could be 
#' rewritten to c(1,1,2,4,3,4,2) providing the object they are indexing into 
#' have been updated.
#' 
#' @param x A vector of indexes
#' 
#' @param index the indexes that have been removed
#' 
#' @return A reformatted x vector
#' 
#' @noRd
#' 
removeIndex <- function(x, index) {
    index <- sort(index, decreasing = TRUE)
    for(i in 1:length(index)) {
        x[x > index[i]] <- x[x > index[i]]-1
    }
    as.integer(x)
}
#' Merge gene groups based on neighborhood similarity
#' 
#' This function merges gene groups from two pangenomes based on similarity of
#' neighborhood and sequence.
#' 
#' @param pangenome The pangenome with the original gene groups
#' 
#' @param newPG The new pangenome with genes to add
#' 
#' @param flankSize The size of the neighborhood to take into account
#' 
#' @param members Grouping of gene groups from pangenome and newPG
#' 
#' @param similarities The similarity matrix used to group the gene groups
#' 
#' @return A vector with gene group membership for all genes in the two 
#' pangenomes
#' 
#' @importFrom dplyr mutate
#' 
#' @noRd
#' 
mergeNeighborhood <- function(pangenome, newPG, flankSize, members, similarities) {
#     newGroups <- seqToGeneGroup(newPG) + nGeneGroups(pangenome)
#     tempGrouping <- c(members[seqToGeneGroup(pangenome)], members[newGroups])
#     tempInfo <- groupInfo(pangenome)
#     tempInfo <- rbind(tempInfo, tempInfo[rep(1, nGeneGroups(newPG)),])
#     tempPG <- mergePangenomes(pangenome, newPG, tempGrouping, tempInfo)
#     tempPG <- neighborhoodSplitting(tempPG)
#     
#     
#     newMembers <- members[(nGeneGroups(pangenome)+1):length(members)]
#     oldMembers <- members[1:nGeneGroups(pangenome)]
#     
#     mergedGroups <- split(1:length(members), members)
#     oldGroupNo <- nGeneGroups(pangenome)
#     
#     oldNeighborhood <- getNeighborhood(pangenome)
#     newNeighborhood <- getNeighborhood(newPG)
#     
#     for(i in mergedGroups) {
#         newGroups <- i > oldGroupNo
#         if(sum(newGroups) == 0) next
#         
#         if(all(newGroups)) {
#             members[i] <- NA
#         } else {
#             oldGroups <- i[!newGroups]
#             newGroups <- i[newGroups]
#             groupMatch <- list()
#             for(j in newGroups) {
#                 neighborhood <- newNeighborhood[[j-oldGroupNo]]
#                 neighborhoodMatch <- sapply(oldNeighborhood[oldGroups], checkNeighborhood, newN=neighborhood, members=members)
#                 oldGroups[neighborhoodMatch]
#             }
#         }
#     }
#     
#     
#     
#     
#     sequenceInfo$gene <- 1:nrow(sequenceInfo)
#     neighbor <- sequenceInfo %>%
#         group_by(contig, organism) %>%
#         arrange(start, end) %>%
#         mutate(backward=collectNeighbors(gene, 'b', flankSize), forward=collectNeighbors(gene, 'f', flankSize)) %>%
#         mutate(backwards=rev(ifelse(strand==1, backward, forward)), forwards=ifelse(strand==1, forward, backward))
#     
#     mergedMembers <- lapply(newMembers, function(x) {
#         which(oldMembers == x)
#     })
#     ggIndex <- seqToGeneGroup(pangenome)
#     genes <- which(ggIndex %in% unlist(mergedMembers))
#     geneGroups <- ggIndex[genes]
#     trails <- trailGroups(genes, pangenome, flankSize)
#     for(i in 1:length(mergedMembers)) {
#         if(length(mergedMembers[[i]]) == 0) {
#             mergedMembers[[i]] <- NA
#         } else {
#             geneTrail <- unlist(strsplit(c(neighbor$backward[i], neighbor$forward[i]), ';'))
#             geneTrail <- as.integer(geneTrail[geneTrail != 'NA'])
#             geneTrail <- unlist(mergedMembers[geneTrail])
#             overlaps <- sapply(mergedMembers[[i]], function(x) {
#                 groupTrail <- trails[geneGroups == x]
#                 sum(unlist(groupTrail) %in% geneTrail)/(length(groupTrail)-0.01) # To favor big groups in tie
#             })
#             if(all(overlaps == 0)) {
#                 mergedMembers[[i]] <- NA
#             } else {
#                 mergedMembers[[i]] <- mergedMembers[[i]][which.max(overlaps)]
#             }
#         }
#     }
#     mergedMembers <- unlist(mergedMembers)
#     newGroups <- split(which(is.na(mergedMembers)), newMembers[is.na(mergedMembers)])
#     mergedMembers[unlist(newGroups)] <- rep(1:length(newGroups), sapply(newGroups, length))+nGeneGroups(pangenome)
#     mergedMembers
}
#' Merge gene groups into the largest
#' 
#' This function merges gene groups from two pangenomes choosing the largest
#' group in the first pangenome when ties occur.
#' 
#' @param pangenome The pangenome to merge new gene groups into
#' 
#' @param members The grouping of gene groups across two pangenomes
#' 
#' @return As mergeNeighborhood
#' 
#' @noRd
#' 
mergeLargest <- function(pangenome, members) {
    newMembers <- members[(nGeneGroups(pangenome)+1):length(members)]
    oldMembers <- members[1:nGeneGroups(pangenome)]
    groupSizes <- groupInfo(pangenome)$nGenes
    mergedMembers <- sapply(newMembers, function(x) {
        links <- which(oldMembers == x)
        if(length(links) == 0) {
            NA
        } else if(length(links) == 1) {
            links
        } else {
            links[which.max(groupSizes[links])]
        }
    })
    newGroups <- split(which(is.na(mergedMembers)), newMembers[is.na(mergedMembers)])
    mergedMembers[unlist(newGroups)] <- rep(1:length(newGroups), sapply(newGroups, length))+nGeneGroups(pangenome)
    mergedMembers
}