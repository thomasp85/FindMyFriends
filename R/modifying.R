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
#' @param pParam A BiocParallelParam object
#' 
#' @param nsParam A list of parameters to pass to 
#' \code{\link{neighborhoodSplit}} or FALSE to skip neighborhood splitting 
#' altogether. If object has had neighborhood splitting performed and nsParam is
#' set to FALSE it is bound to cause problems, so don't do that.
#' 
#' @param klParam A list of parameters to pass to \code{\link{kmerLink}} or 
#' FALSE to skip paralogue linking altogether. Independent of the value of 
#' klParam kmerLink will only be run if paralogue links have been defined on 
#' object beforehand.
#' 
#' @importFrom kebabs getExRep spectrumKernel linearKernel
#' @importFrom dplyr bind_rows
#' @importFrom BiocParallel bpworkers
#' 
setMethod(
    'addGenomes', c('pgVirtual', 'pgVirtual'),
    function(object, newSet, kmerSize, lowerLimit, pParam, nsParam = list(), 
             klParam = list()) {
        .fillDefaults(defaults(object))
        
        if (class(object) != class(newSet)) {
            stop('pangenome and newSet must be of the same class')
        }
        if (translated(object) != translated(newSet)) {
            stop('Translational status of pangenome and newSet must match')
        }
        
        if (!hasGeneGroups(newSet)) {
            gpcParam <- list(kmerSize = kmerSize, lowerLimit = lowerLimit)
            if (!missing(pParam)) {
                gpcParam$pParam <- pParam
            }
            defaults(newSet) <- defaults(object)
            gpcParam$object <- newSet
            newSet <- do.call(gpcGrouping, gpcParam)
        }
        
        fullSet <- c(getRep(object, 'random'), getRep(newSet, 'random'))
        er <- getExRep(fullSet, spectrumKernel(kmerSize))
        if (missing(pParam)) {
            sim <- linearKernel(er, sparse = TRUE, diag = FALSE, 
                                lowerLimit = lowerLimit)
        } else {
            sim <- lkParallel(er, pParam, bpworkers(pParam), 
                              lowerLimit = lowerLimit)
        }
        gr <- graph_from_adjacency_matrix(sim, mode = 'lower', 
                                          weighted = TRUE, 
                                          diag = FALSE)
        members <- components(gr)$membership
        groups <- c(seqToGeneGroup(object), 
                    seqToGeneGroup(newSet) + nGeneGroups(object))
        groups <- split(seq_along(groups), groups)
        groups <- lapply(split(groups, members), unlist)
        if (hasGeneInfo(object) && !(is.logical(nsParam) && !nsParam)) {
            info <- resizeDataFrame(groupInfo(object), length(groups))
            newPG <- mergePangenomes(object, newSet, convertGrouping(groups), 
                                     info)
            nsParam$object <- newPG
            nsParam$guideGroups <- seqToGeneGroup(object)
            newPG <- do.call(neighborhoodSplit, nsParam)
            groups <- convertGrouping(seqToGeneGroup(newPG))
        }
        newGroups <- matchGroups(groups, seqToGeneGroup(object))
        removedGroups <- which(!1:max(newGroups) %in% unique(newGroups))
        addedGroups <- max(newGroups) - nGeneGroups(object)
        addedGroupNames <- paste0(defaults(object)$groupPrefix,
                                  seq(object@.settings$nextGroup, 
                                      length.out = addedGroups))
        newInfo <- resizeDataFrame(groupInfo(object), max(newGroups), 
                                   addedGroupNames)[-removedGroups, ]
        newGroups <- removeIndex(newGroups)
        newPG <- mergePangenomes(object, newSet, newGroups, newInfo)
        
        if (hasParalogueLinks(object) && !(is.logical(klParam) && !klParam)) {
            oldLinks <- groupInfo(newPG)$paralogue
            klParam$object <- newPG
            newPG <- do.call(kmerLink, klParam)
            newLinks <- split(1:nGeneGroups(newPG), groupInfo(newPG)$paralogue)
            newLinks <- matchGroups(newLinks, oldLinks)
            groupInfo(newPG)$paralogue <- newLinks
        }
        
        reportGroupChanges(groupNames(newPG)[seqToGeneGroup(newPG)], 
                           groupNames(object)[seqToGeneGroup(object)])
        
        newPG
    }
)
#' @describeIn removeGene Remove gene based on gene name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'missing', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- which(geneNames(object) %in% name)
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene and organism name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'character', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, name = name, organism = index)
    }
)
#' @describeIn removeGene Remove gene based on gene name and organism index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'character', 'numeric', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(n, o) {
            which(n == geneNames(object) & o == seqToOrg(object))
        }, n = name, o = organism))
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
        }
    }
)
#' @describeIn removeGene Remove gene based on organism name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'character', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, organism = index)
    }
)
#' @describeIn removeGene Remove gene based on organism index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'numeric', 'missing', 'missing'),
    function(object, name, organism, group, ind) {
        index <- findIn(as.integer(organism), seqToOrg(object))
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
        }
    }
)
#' @describeIn removeGene Remove gene based on organism name and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'character', 'missing', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- match(organism, orgNames(object))
        removeGene(object, organism = index, ind = ind)
    }
)
#' @describeIn removeGene Remove gene based on organism and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'numeric', 'missing', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(o, i) {
            findIn(o, seqToOrg(object))[i]
        }, o = as.integer(organism), i = ind))
        index <- index[!is.na(index)]
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene group name
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'character', 'missing'),
    function(object, name, organism, group, ind) {
        index <- match(group, groupNames(object))
        removeGene(object, group = index)
    }
)
#' @describeIn removeGene Remove gene based on gene group index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'numeric', 'missing'),
    function(object, name, organism, group, ind) {
        index <- findIn(as.integer(group), seqToGeneGroup(object))
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
        }
    }
)
#' @describeIn removeGene Remove gene based on gene group name and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'character', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- match(group, groupNames(object))
        removeGene(object, group = index, ind = ind)
    }
)
#' @describeIn removeGene Remove gene based on gene group and gene index
#' 
setMethod(
    'removeGene', c('pgVirtual', 'missing', 'missing', 'numeric', 'numeric'),
    function(object, name, organism, group, ind) {
        index <- unlist(mapply(function(g, i) {
            findIn(g, seqToGeneGroup(object))[i]
        }, g = as.integer(group), i = ind))
        index <- index[!is.na(index)]
        if (length(index) == 0) {
            warning('No genes match criteria')
            object
        } else {
            removeGene(object, ind = index)
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
        if (!hasParalogueLinks(object)) {
            stop('Paralogues must be detected before collapsing')
        }
        paralogues <- split(1:nGeneGroups(object), groupInfo(object)$paralogue)
        groups <- split(1:nGenes(object), seqToGeneGroup(object))
        newGroups <- lapply(paralogues, function(x) {unlist(groups[x])})
        newInfo <- lapply(paralogues, function(x) {
            if (inherits(combineInfo, 'function')) {
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
        if (missing(key)) {
            key <- 1:nGeneGroups(object)
        } else if (inherits(key, 'character')) {
            keyInfo <- info[[key[1]]]
            info <- info[, names(info) != key[1], drop = FALSE]
            if (is.numeric(keyInfo)) {
                key <- keyInfo
            } else {
                key <- match(keyInfo, groupNames(object))
            }
        }
        if (!is.numeric(key)) stop('key must be convertible to an index')
        for (i in names(info)) {
            if (i %in% c('nGenes', 'nOrg', 'group')) next
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
        if (missing(key)) {
            key <- 1:nOrganisms(object)
        } else if (inherits(key, 'character')) {
            keyInfo <- info[[key[1]]]
            info <- info[, names(info) != key[1], drop = FALSE]
            if (is.numeric(keyInfo)) {
                key <- keyInfo
            } else {
                key <- match(keyInfo, orgNames(object))
            }
        }
        if (!is.numeric(key)) stop('key must be convertible to an index')
        for (i in names(info)) {
            if (i %in% c('nGenes')) next
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
    as.data.frame(lapply(data, function(x, sep) {
        info <- unique(unlist(x))
        if (missing(sep)) {
            I(list(info))
        } else {
            paste(info, collapse = sep)
        }
    }, sep = sep), stringsAsFactors = FALSE)
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
    data[largest, , drop = FALSE]
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
    nrow*(col - 1) + row
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
#' @return A reformatted x vector
#' 
#' @noRd
#' 
removeIndex <- function(x) {
    match(x, sort(unique(x)))
}
#' Match new groups to old ones based on members
#' 
#' This function takes a grouping of genes and compares it to an older, possibly
#' smaller, grouping (i.e. based on fewer genes), and returns a reindexed
#' grouping so that group indexes matches between new and old.
#' 
#' @param newGroups A list of integer vectors with the new groups
#' 
#' @param oldGroups An integer vector with the membership of genes in the old
#' grouping
#' 
#' @return An integer vector with the membership of genes in the new set
#' 
#' @noRd
#' 
matchGroups <- function(newGroups, oldGroups) {
    finalGrouping <- rep(NA, length(unlist(newGroups)))
    pendingGroups <- seq_along(newGroups)
    bestGroups <- lapply(newGroups, function(group) {
        table(oldGroups[group], useNA = 'no')
    })
    memberGroup <- rep(seq_along(bestGroups), lengths(bestGroups))
    oldGroupInd <- as.integer(unlist(lapply(bestGroups, names)))
    bestGroups <- unlist(bestGroups)
    bestGroupOrder <- order(bestGroups, decreasing = TRUE)
    for (i in seq_along(bestGroupOrder)) {
        oldGroup <- oldGroupInd[bestGroupOrder[i]]
        if (is.na(oldGroup)) next
        newGroup <- memberGroup[bestGroupOrder[i]]
        finalGrouping[newGroups[[newGroup]]] <- oldGroup
        oldGroupInd[oldGroupInd == oldGroup] <- NA
        oldGroupInd[memberGroup == newGroup] <- NA
        pendingGroups[newGroup] <- NA
    }
    missingGroups <- !is.na(pendingGroups)
    if (any(missingGroups)) {
        nextGroup <- max(finalGrouping, na.rm = TRUE) + 1
        addedGroups <- rep(seq(from = nextGroup, 
                               length.out = sum(missingGroups)),
                           lengths(newGroups[missingGroups]))
        finalGrouping[unlist(newGroups[missingGroups])] <- addedGroups
    }
    finalGrouping
}
#' Reports the change in grouping
#' 
#' This function inspects gene grouping before and after a change and reports on
#' the changes. If newGrouping is missing it reports on the last performed 
#' comparison; optionally writing it to a file if 'file' is specified. 
#' 
#' @param newGrouping An integer vector as produced by 
#' \code{\link{seqToGeneGroup}} with the grouping after the change
#' 
#' @param oldGrouping An integer vector as produced by 
#' \code{\link{seqToGeneGroup}} with the grouping before the change
#' 
#' @param file A file to write 
#' 
#' @return This function is called for its side effects
#' 
#' @examples 
#' # Show latest changes in grouping
#' reportGroupChanges()
#' 
#' # Alternatively write it to a file
#' reportGroupChanges(file = tempfile())
#' 
#' @export
#' 
reportGroupChanges <- function(newGrouping, oldGrouping, file) {
    if (missing(newGrouping)) {
        if (is.null(.pkg_variables$report)) {
            .pkg_variables$report <- c()
        }
        if (missing(file)) {
            for (i in .pkg_variables$report) {
                message(i)
            }
        } else {
            write(.pkg_variables$report, file = file)
        }
        return(invisible())
    }
    changes <- newGrouping[seq_along(oldGrouping)] != oldGrouping
    if (sum(changes) == 0) return()
    report <- c()
    groups <- split(which(changes), oldGrouping[changes])
    for (i in seq_along(groups)) {
        movedTo <- split(groups[[i]], newGrouping[groups[[i]]])
        for (j in seq_along(movedTo)) {
            report <- c(report, paste0('Gene ', 
                                       paste(movedTo[[j]], collapse = ', '), 
                                       ' moved from group ', names(groups)[i], 
                                       ' to ', names(movedTo)[j]))
        }
    }
    if (length(report) > 50) {
        message('More than 50 gene group changes. Use reportGroupChanges() to ',
                'see them all, or reportGroupChanges(file="your/file.txt" to ',
                'write them to a file.')
    } else {
        if (missing(file)) {
            for (i in report) {
                message(i)
            }
        } else {
            write(report, file = file)
        }
    }
    extraGroups <- unique(newGrouping[!newGrouping %in% unique(oldGrouping)])
    removedGroups <- unique(oldGrouping[!oldGrouping %in% unique(newGrouping)])
    if (length(extraGroups) != 0) {
        extraMessage <- paste0(length(extraGroups), ' new groups added')
        if (missing(file)) {
            message(extraMessage)
        } else {
            write(extraMessage, file = file, append = TRUE)
        }
        report <- c(report, extraMessage)
    }
    if (length(removedGroups) != 0) {
        removedMessage <- paste0(length(removedGroups), ' empty groups removed')
        if (missing(file)) {
            message(removedMessage)
        } else {
            write(removedMessage, file = file, append = TRUE)
        }
        report <- c(report, removedMessage)
    }
    assign('report', report, envir = .pkg_variables)
    invisible()
}
#' Resize the number of rows in a data frame
#' 
#' This function allows for quick resizing of data.frames either shrinking by
#' removing from the bottom, or growing by adding NA rows
#' 
#' @noRd
#' 
#' @importFrom dplyr bind_rows
#' 
resizeDataFrame <- function(df, nrows, newRowNames = NULL) {
    if (nrows == nrow(df)) {
        return(df)
    } else if (nrows < nrow(df)) {
        return(df[1:nrows, , drop = FALSE])
    }
    addedRows <- nrows - nrow(df)
    appender <- data.frame(rep(NA, addedRows))
    names(appender) <- names(df)[1]
    newDF <- as.data.frame(bind_rows(df, appender))
    if (is.null(newRowNames)) {
        newRowNames <- as.character((nrow(df) + 1):(nrow(newDF)))
    }
    rownames(newDF) <- make.unique(c(rownames(df), newRowNames))
    newDF
}