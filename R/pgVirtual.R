#' @include generics.R
#' @include aaa.R
NULL

#' Base class for pangenomic data
#' 
#' This virtual class is the superclass of all other pangenome classes in 
#' FindMyFriends. It is an empty shell that is mainly used for dispatch and 
#' checking that the promises of subclasses are held.
#' 
#' Subclasses of pgVirtual must implement the following methods in order for 
#' them to plug into FindMyFriends algorithms:
#' 
#' \describe{
#'  \item{seqToOrg(object)}{Returns the mapping from genes to organisms as an 
#'  integer vector with position mapped to gene and integer mapped to organism.}
#'  \item{seqToGeneGroup(object)}{As seqToOrg but mapped to gene group instead 
#'  of organism. If gene groups are yet to be defined return an empty vector.}
#'  \item{genes(object, split, subset)}{Return the underlying sequences. If 
#'  split is missing return an XStringSet, otherwise return an XStringSetList. 
#'  split can be either 'group', 'organism' or 'paralogue' and should group the 
#'  sequences accordingly. Subset should behave as if it was added as '[]' to 
#'  the results but allow you to avoid reading everything into memory if not 
#'  needed.}
#'  \item{geneNames(object)}{Return a character vector with the name of each 
#'  gene.}
#'  \item{geneNames<-(object, value)}{Set the name of each gene.}
#'  \item{geneWidth(object)}{Return an integer vector with the length (in 
#'  residues) of each gene.}
#'  \item{removeGene(object, name, organism, group, ind)}{\strong{Should only be
#'  implemented for signature: c(\emph{yourClass}, 'missing', 'missing', 
#'  'missing', 'integer')} Remove the genes at the given indexes and return the 
#'  object.}
#'  \item{orgNames(object)}{Return a character vector of organism names.}
#'  \item{orgNames<-(object, value)}{Set the name of the organisms.}
#'  \item{groupNames(object)}{Return a character vector of gene group names.}
#'  \item{groupNames<-(object, value)}{Set the name of gene groups.}
#'  \item{orgInfo(object)}{Return a data.frame with metadata about each 
#'  organim.}
#'  \item{orgInfo<-(object, value)}{Set a data.frame to be metadata about each
#'  organism.}
#'  \item{setOrgInfo(object, name, info, key)}{Set the metadata 'name', for the
#'  organisms corresponding to 'key' to 'info'}
#'  \item{groupInfo(object)}{Return a data.frame with metadata about each 
#'  gene group.}
#'  \item{groupInfo<-(object, value)}{Set a data.frame to be metadata about each
#'  gene group.}
#'  \item{setGroupInfo(object, name, info, key)}{Set the metadata 'name', for 
#'  the gene groups corresponding to 'key' to 'info'}
#'  \item{groupGenes(object, seqToGeneGroup)}{Sets the gene grouping of the 
#'  pangenome. 'seqToGeneGroup' should correspond to the output of the 
#'  seqToGeneGroup method (i.e. an integer vector with each element giving the
#'  group of the corresponding gene). This method \strong{must} include a 
#'  \code{callNextMethod(object)} as the last line.}
#'  \item{mergePangenomes(pg1, pg2, geneGrouping, groupInfo)}{Merge pg2 into pg1 
#'  preserving the indexing in pg1 and appending and modifying the indexing of 
#'  pg2. The geneGrouping argument is the new grouping of genes and groupInfo
#'  the new group info for the groups.}
#' }
#' 
#' Additionally subclasses can override the following methods for performance
#' gains. Otherwise they will be derived from the above methods.
#' 
#' \describe{
#'  \item{length(object)}{Return the number of organisms in the object.}
#'  \item{nOrganisms(object)}{As length.}
#'  \item{nGenes(object)}{Return the number of genes in the object.}
#'  \item{nGeneGroups(object)}{Return the number of gene groups}
#'  \item{hasGeneGroups}{Returns TRUE if gene groups have been defined}
#'  \item{pgMatrix}{Returns an integer matrix with organisms as columns and gene
#'  groups as rows and the corresponding number of genes in each element.}
#' }
#' 
#' Developers are encourages to consult the implementation of FindMyFriends own
#' classes when trying to implement new ones
#' 
#' @slot .settings A list containing settings pertaining to the object
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgVirtual',
    contains = 'VIRTUAL',
    slots = list(
        .settings = 'list'
    ),
    validity = function(object) {
        # Test for methods
        if (!hasMethod('seqToOrg', class(object))) {
            return('The method "seqToOrg" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('seqToGeneGroup', class(object))) {
            return('The method "seqToGeneGroup" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('genes', c(class(object), 'missing'))) {
            return('The method "genes" must be implemented for ', 
                   class(object), 
                   ', with signature "', 
                   class(object),
                   ', missing"')
        }
        if (!hasMethod('genes', c(class(object), 'character'))) {
            return('The method "genes" must be implemented for ', 
                   class(object), 
                   ', with signature "', 
                   class(object),
                   ', character"')
        }
        if (!hasMethod('geneNames', class(object))) {
            return('The method "geneNames" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('geneNames<-', class(object))) {
            return('The method "geneNames<-" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('geneWidth', class(object))) {
            return('The method "geneWidth" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('orgNames', class(object))) {
            return('The method "orgNames" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('orgNames<-', class(object))) {
            return('The method "orgNames<-" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('groupNames', class(object))) {
            return('The method "groupNames" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('groupNames<-', class(object))) {
            return('The method "groupNames<-" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('orgInfo', class(object))) {
            return('The method "orgInfo" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('orgInfo<-', class(object))) {
            return('The method "orgInfo<-" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('groupInfo', class(object))) {
            return('The method "groupInfo" must be implemented for ', 
                   class(object))
        }
        if (!hasMethod('groupInfo<-', class(object))) {
            return('The method "groupInfo<-" must be implemented for ', 
                   class(object))
        }
        
        # Test for object sanity
        if (length(object) != 0) {
            minOrg <- max(seqToOrg(object))
            oNames <- orgNames(object)
            oInfo <- orgInfo(object)
            if (length(oNames) != nrow(oInfo) || length(oNames) < minOrg) {
                return('Organism indexing mismatch')
            }
            if (hasGeneGroups(object)){
                minGroup <- max(seqToGeneGroup(object))
                gNames <- groupNames(object)
                gInfo <- groupInfo(object)
                if (length(gNames) != nrow(gInfo) || 
                    length(gNames) < minGroup) {
                    return('Gene group indexing mismatch')
                }
            }
        }
        return(TRUE)
    },
    prototype = list(
        .settings = list(translated = FALSE)
    )
)

### UTILITY FUNCTIONS

#' @describeIn defaults Default values for pgVirtual subclass objects
#' 
setMethod(
    "defaults", 'pgVirtual', 
    function(object) {
        object@.settings[!names(object@.settings) %in% c('translated', 
                                                         'nextGroup')]
    }
)
#' @describeIn defaults Set defaults for pgVirtual subclass objects
#' 
setMethod(
    'defaults<-', 'pgVirtual',
    function(object, value) {
        if ('translated' %in% names(value)) {
            stop('Translational status cannot be redefined')
        }
        if ('nextGroup' %in% names(value)) {
            stop('Group names are handled automatically')
        }
        value$translated <- translated(object)
        value$nextGroup <- object@.settings$nextGroup
        oldThreshold <- defaults(object)$coreThreshold
        object@.settings <- value
        if (value$coreThreshold != oldThreshold && nGeneGroups(object) != 0) {
            newGroups <- calcGroupInfo(split(seqToOrg(object), 
                                             seqToGeneGroup(object)), 
                                       nOrganisms(object), value$coreThreshold)
            groupInfo(object)$group <- newGroups$group
        }
        object
    }
)
#' @describeIn pgVirtual Length of a Pangenome, defined as the number of 
#' organisms it contain
#' 
#' @param x A pgVirtual subclass object
#' 
#' @return Length returns an integer giving the number of organisms
#' 
setMethod(
    'length', 'pgVirtual',
    function(x) {
        length(orgNames(x))
    }
)

#' @describeIn pgVirtual Basic information about the pangenome
#' 
#' @param object A pgVirtual subclass object
#' 
setMethod(
    'show', 'pgVirtual',
    function(object) {
        cat('An object of class ', class(object), '\n\n', sep = '')
        cat('The pangenome consists of', 
            nGenes(object), 
            'genes from', 
            nOrganisms(object), 
            'organisms\n')
        if (hasGeneGroups(object)) {
            cat(nGeneGroups(object), 'gene groups defined\n\n')
            groupSizes <- groupDistribution(object)
            cat('     Core|', createProgress(groupSizes['Core']), '\n', 
                sep = '')
            cat('Accessory|', createProgress(groupSizes['Accessory']), '\n', 
                sep = '')
            cat('Singleton|', createProgress(groupSizes['Singleton']), '\n\n', 
                sep = '')
        } else {
            cat('Gene groups not yet defined\n\n')
        }
        if (translated(object)) {
            cat('Genes are translated\n')
        } else {
            cat('Genes are untranslated\n')
        }
    }
)
#' @describeIn pgVirtual Create subsets of pangenomes based on index
#' 
#' @param i indices specifying genomes, either integer, numeric, character or 
#' logical, following the normal rules for indexing objects in R
#' 
setMethod(
    '[', c('pgVirtual', 'integer'),
    function(x, i) {
        if (anyDuplicated(i)) {
            warning('ignoring duplicated indices')
        }
        if (any(i > nOrganisms(x))) {
            warning('ignoring indices out of bound')
        }
        remove <- which(!seq_along(x) %in% i)
        removeGene(x, organism = remove)
    }
)
#' @describeIn pgVirtual Create subsets of pangenomes based on index
#' 
setMethod(
    '[', c('pgVirtual', 'numeric'),
    function(x, i) {
        i <- as.integer(i)
        if (any(is.na(i))) {
            stop('index not convertible to integer')
        }
        x[i]
    }
)
#' @describeIn pgVirtual Create subsets of pangenomes based on organism name
#' 
setMethod(
    '[', c('pgVirtual', 'character'),
    function(x, i) {
        i <- match(i, orgNames(x))
        if (any(is.na(i))) {
            stop('Organism names not present')
        }
        x[i]
    }
)
#' @describeIn pgVirtual Create subsets of pangenomes based on logical vector
#' 
setMethod(
    '[', c('pgVirtual', 'logical'),
    function(x, i) {
        i <- rep(i, length.out = nOrganisms(x))
        x[which(i)]
    }
)
#' @describeIn pgVirtual Extract sequences from a single organism
#' 
setMethod(
    '[[', 'pgVirtual',
    function(x, i) {
        if (length(i) != 1) {
            stop('i must be of length 1')
        }
        genes(x, split = 'organism', subset = i)[[1]]
    }
)
#' @describeIn nGenes The number of genes in the pangenome for pgVirtual 
#' subclasses.
#' 
setMethod(
    'nGenes', 'pgVirtual',
    function(object) {
        length(seqToOrg(object))
    }
)
#' @describeIn nOrganisms The number of organisms in the pangenome for pgVirtual
#' subclasses.
#' 
setMethod(
    'nOrganisms', 'pgVirtual',
    function(object) {
        length(object)
    }
)

#' @describeIn nGeneGroups The number of gene groups in the pangenome for
#' pgVirtual subclasses
#' 
setMethod(
    'nGeneGroups', 'pgVirtual',
    function(object) {
        nrow(groupInfo(object))
    }
)

#' @describeIn hasGeneGroups Gene group check for pgVirtual subclasses
#' 
setMethod(
    'hasGeneGroups', 'pgVirtual',
    function(object) {
        nGeneGroups(object) > 0
    }
)
#' @rdname internalGroupGenes
#' 
#' @param seqToGeneGroup Integer vector equivalent to the output of 
#' seqToGeneGroup on the new object.
#' 
setMethod(
    'groupGenes', 'pgVirtual',
    function(object, seqToGeneGroup) {
        seqToOrg <- seqToOrg(object)
        groups <- split(seqToOrg, seqToGeneGroup)
        nOrgs <- nOrganisms(object)
        groupInfoCalc <- calcGroupInfo(groups, nOrgs, defaults(object)$coreThreshold)
        groupInfo <- data.frame(
            description = NA, 
            group = NA, 
            paralogue = NA, 
            GO = NA, 
            EC = NA, 
            nOrg = rep(0, max(seqToGeneGroup)), 
            nGenes = rep(0, max(seqToGeneGroup)),
            row.names = paste0(defaults(object)$groupPrefix,
                              1:max(seqToGeneGroup))
        )
        groupInfo[as.integer(names(groups)), 
                  c('group', 'nOrg', 'nGenes')] <- groupInfoCalc
        groupInfo(object) <- groupInfo
        object@.settings$nextGroup <- max(seqToGeneGroup) + 1
        object
    }
)
#' @describeIn hasGeneInfo Checks whether gene location information is available
#' for pgVirtual subclasses
#' 
setMethod(
    'hasGeneInfo', 'pgVirtual',
    function(object) {
        inherits(object, 'pgVirtualLoc')
    }
)

#' @describeIn hasParalogueLinks Check for secondary gene grouping in pgVirtual
#' subclasses
#' 
setMethod(
    'hasParalogueLinks', 'pgVirtual',
    function(object) {
        if (hasGeneGroups(object)) {
            !any(is.na(groupInfo(object)$paralogue))
        } else {
            FALSE
        }
    }
)

#' @describeIn translated Sequence type check for pgVirtual subclasses
#' 
setMethod(
    'translated', 'pgVirtual',
    function(object) {
        object@.settings$translated
    }
)

#' @describeIn pgMatrix Get pangenome matrix for pgVirtual subclasses
#' 
setMethod(
    'pgMatrix', 'pgVirtual',
    function(object) {
        getPgMatrix(object)
    }
)
#' @rdname pgVirtual-class
#' 
#' @usage as(object, Class='ExpressionSet')
#' 
#' @param Class The class to coerce pgVirtual subclasses to. Outside of the
#' FindMyFriends class tree only 'ExpressionSet' and 'matrix' is implemented.
#' 
#' @name as
#' 
#' @importFrom Biobase ExpressionSet
#' 
setAs(
    'pgVirtual', 'ExpressionSet',
    function(from) {
        ExpressionSet(
            as.matrix(pgMatrix(from)),
            as(orgInfo(from), 'AnnotatedDataFrame'),
            as(groupInfo(from), 'AnnotatedDataFrame')
        )
    }
)
#' @rdname pgVirtual-class
#' 
#' @usage as(object, Class='matrix')
#' 
#' @name as
#' 
setAs(
    'pgVirtual', 'matrix',
    function(from) {
        as.matrix(pgMatrix(from))
    }
)
#' @describeIn getRep Get a representative sequence for each gene group for 
#' pgVirtual subclasses
#' 
setMethod(
    'getRep', c('pgVirtual', 'character'),
    function(object, method) {
        ind <- split(1:nGenes(object), seqToGeneGroup(object))
        rep <- switch(
            method,
            random = {
                ind <- sapply(ind, function(x) {x[sample(length(x), size = 1)]})
                genes(object, subset=ind)
            },
            longest = {
                ind <- sapply(ind, function(x, width) {
                    x[which.max(width[x])]
                }, width = geneWidth(object))
                genes(object, subset=ind)
            },
            shortest = {
                ind <- sapply(ind, function(x, width) {
                    x[which.min(width[x])]
                }, width = geneWidth(object))
                genes(object, subset=ind)
            },
            stop('Unknown method: ', method)
        )
        names(rep) <- groupNames(object)
        rep
    }
)

### Plotting functions

#' @describeIn plotStat Plot basic statistics for pgVirtual subclasses
#' 
#' @param sort logical. Should Genomes be sorted based on their number of genes
#' 
#' @param color A metadata name to color the organisms by
#' 
#' @importFrom ggplot2 ggplot theme_bw scale_y_continuous theme element_text geom_bar aes_string scale_fill_brewer element_blank scale_fill_manual coord_polar ggplotGrob
#' @importFrom grid grid.newpage grid.draw
#' 
setMethod(
    'plotStat', 'pgVirtual',
    function(object, sort = TRUE, color, ...) {
        data <- orgInfo(object)
        lev <- if (sort) {
            orgNames(object)[order(data$nGenes)]
        } else {
            orgNames(object)
        }
        data$organism <- factor(orgNames(object), levels = lev)
        p <- ggplot() + theme_bw() + scale_y_continuous('# genes')
        p <- p + theme(axis.text.x = element_text(angle = 45, 
                                                  vjust = 1, 
                                                  hjust = 1))
        if (missing(color)) {
            p <- p + geom_bar(aes_string(x = 'organism', y = 'nGenes'), stat = 'identity', 
                              data = data)
        } else {
            p <- p + geom_bar(aes_string(x = 'organism', y = 'nGenes', 
                                         fill = color), 
                              stat = 'identity', data = data)
            if (length(list(...)) != 0) {
                p <- p + scale_fill_brewer(...)
            }
        }
        if (hasGeneGroups(object)) {
            groupNames <- c('Singleton', 'Accessory', 'Core')
            groups <- table(groupInfo(object)$group)
            groups <- data.frame(nGenes = as.integer(groups), 
                                 group = factor(names(groups), 
                                                levels = groupNames))
            groups <- groups[order(groups$group),]
            labPos <- groups$nGenes/2 + cumsum(c(0, groups$nGenes[-nrow(groups)]))
            p1 <- ggplot() + theme_bw()
            p1 <- p1 + theme(axis.title = element_blank(), 
                             axis.text.y = element_blank(), 
                             axis.ticks.y = element_blank(),
                             axis.line = element_blank(),
                             panel.grid = element_blank(),
                             panel.border = element_blank())
            p1 <- p1 + geom_bar(aes_string(x = factor(1), 
                                    fill = 'group', weight = 'nGenes'), 
                                data = groups, width = 1)
            p1 <- p1 + scale_fill_manual('Group', breaks = groupNames, 
                                         values = c('goldenrod', 'forestgreen', 
                                                    'firebrick'), 
                                         drop = FALSE)
            p1 <- p1 + coord_polar(theta = 'y')
            p1 <- p1 + scale_y_continuous(breaks = labPos, 
                                          labels = groups$nGenes)
            
            p <- rbindGtable(ggplotGrob(p1), ggplotGrob(p))
            grid.newpage()
            grid.draw(p)
            invisible(p)
        } else {
            p
        }
    }
)
#' @describeIn plotEvolution Evolution plots for pgVirtual subclasses
#' 
#' @param ordering An ordering of the organisms by name or index, or alternative
#' one of 'bootstrap', 'random' or 'none'.
#' 
#' @param times The number of sampling for ordering='bootstrap'
#' 
#' @importFrom ggplot2 ggplot aes_string theme_bw scale_color_manual scale_y_continuous geom_smooth scale_x_continuous geom_line scale_x_discrete
#' 
setMethod(
    'plotEvolution', 'pgVirtual',
    function(object, ordering = 'bootstrap', times = 10) {
        if (length(ordering) != 1) {
            if (inherits(ordering, 'character')) {
                ordering <- match(ordering, orgNames(object))
            }
        }
        evol <- switch(
            as.character(ordering[1]),
            bootstrap = evolBoot(object, times = times),
            random = evolMan(object, sample(length(object))),
            none = evolMan(object, seq(length(object))),
            evolMan(object, ordering)
        )
        p <- ggplot(evol, aes_string(x = 'org', 
                              y = 'size', 
                              color = 'group', 
                              group = 'group')) + theme_bw()
        p <- p + scale_color_manual('',
                                    values = c(Singleton = 'goldenrod', 
                                               Accessory = 'forestgreen', 
                                               Core = 'firebrick', 
                                               Total = 'steelblue'), 
                                    breaks = c('Singleton', 
                                               'Accessory', 
                                               'Core', 
                                               'Total'))
        p <- p + scale_y_continuous('# Gene groups')
        if (ordering[1] == 'bootstrap') {
            p <- p + geom_smooth(size = 1.5)
            p <- p + scale_x_continuous('# Organisms')
        } else {
            p <- p + geom_line(size = 1.5)
            p <- p + scale_x_discrete('Organism added')
        }
        p
    }
)
#' @describeIn plotSimilarity Similarity heatmaps for pgVirtual subclasses
#' 
#' @param type The type of similarity calculation. Either 'pangenome' or 'kmer'
#' 
#' @param ordering The ordering of rows and column in the heatmap. Either 
#' integer og character vector with organism names or one of the following:
#' 'auto' or 'none'. For 'auto' The ordering will be based on a hierarchical 
#' clustering of the organisms based on Ward's distance.
#' 
#' @param kmerSize The size of the kmers to use for comparison
#' 
#' @param pParam An object of class BiocParallelParam
#' 
#' @param chunkSize Number of organisms to process at a time
#' 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string theme_bw theme element_text geom_raster scale_x_discrete scale_y_discrete coord_fixed scale_fill_distiller
#' @importFrom stats as.dist hclust
#' 
setMethod(
    'plotSimilarity', 'pgVirtual',
    function(object, type = 'pangenome', ordering = 'auto', kmerSize, pParam, 
             chunkSize = 100) {
        .fillDefaults(defaults(object))
        
        sim <- switch(
            type,
            pangenome = pgSim(object),
            kmer = kmerSim(object, kmerSize, chunkSize, pParam)
        )
        if (ordering[1] == 'auto') {
            dist <- as.dist(sqrt(1 - sim))
            ordering <- hclust(dist, method = 'ward.D2')$order
        } else if (ordering[1] == 'none') {
            ordering <- seq(length(object))
        } else if (inherits(ordering, 'character')) {
            ordering <- match(ordering, orgNames(object))
        }
        sim <- melt(sim, varnames = c('org1', 'org2'), 
                    value.name = 'Similarity')
        sim$org1 <- factor(sim$org1, levels = orgNames(object)[ordering])
        sim$org2 <- factor(sim$org2, levels = rev(orgNames(object)[ordering]))
        p <- ggplot(sim, aes_string(x = 'org1', y = 'org2', fill = 'Similarity')) + 
            theme_bw()
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                                  vjust = 1))
        p <- p + geom_raster()
        p <- p + scale_x_discrete('Organism', expand = c(0,0))
        p <- p + scale_y_discrete('Organism', expand = c(0,0))
        p <- p + coord_fixed()
        p <- p + scale_fill_distiller(palette = 7, guide = 'colorbar', 
                                      limits = c(0,1), direction = 1)
        p
    }
)
#' @describeIn plotTree Dendrogram plotting for pgVirtual subclasses
#' 
#' @param type The type of distance meassure. Either 'pangenome' or 'kmer'
#' 
#' @param circular logical. Should the dendrogram be plotted in a circular
#' fashion.
#' 
#' @param info Metadata to plot at the leafs of the dendrogram. Either the name
#' of an orgInfo column or a vector with information for each organism.
#' 
#' @param kmerSize The size of the kmers to use for comparison
#' 
#' @param dist The distance function to use. All possible values of method in 
#' dist() are allowed as well as 'cosine' for type='kmer'
#' 
#' @param clust The clustering function to use. Passed on to method in hclust()
#' 
#' @param pParam An object of class BiocParallelParam
#' 
#' @param chunkSize Number of organisms to process at a time
#' 
#' @importFrom ggdendro dendro_data label
#' @importFrom ggplot2 ggplot aes_string theme_bw theme element_blank geom_segment coord_polar scale_x_continuous scale_y_continuous element_text 
#' 
setMethod(
    'plotTree', 'pgVirtual',
    function(object, type = 'pangenome', circular = FALSE, info, kmerSize, dist, 
             clust, pParam, chunkSize = 100) {
        .fillDefaults(defaults(object))
        
        tree <- orgTree(object, type, kmerSize, dist, clust, chunkSize = 100, 
                        pParam)
        data <- dendro_data(tree)
        if (!missing(info)) {
            infoDat <- data.frame(x = label(data)$x - 0.5, 
                                  xend = label(data)$x + 0.5, 
                                  y = label(data)$y, 
                                  yend = label(data)$y)
            orgOrder <- match(label(data)$label, orgNames(object))
            if (length(info) == 1) {
                if (!any(names(orgInfo(object)) == info)) {
                    stop('info doesn\'t match any org info')
                }
                infoDat[[info]] <- orgInfo(object)[orgOrder, info]
            } else {
                infoDat[['.PLCHLDR']] <- info[orgOrder]
                info <- '.PLCHLDR'
            }
        }
        if (circular) {
            data$segments$y <- -data$segments$y
            data$segments$yend <- -data$segments$yend
        }
        p <- ggplot(data$segments, 
                    aes_string(x = 'x', xend = 'xend', y = 'y', yend = 'yend')) + theme_bw()
        p <- p + theme(axis.title = element_blank(), 
                       axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(),
                       axis.line = element_blank(),
                       panel.grid = element_blank(),
                       panel.border = element_blank())
        p <- p + geom_segment()
        if (!missing(info)) {
            p <- p + geom_segment(aes_string(color = info), data = infoDat, 
                                  size = I(4))
            if (info == '.PLCHLDR') {
                p <- p + theme(legend.title = element_blank())
            }
        }
        if (circular) {
            p <- p + coord_polar()
            p <- p + scale_x_continuous(breaks = data$labels$x, 
                                        labels = data$labels$label, 
                                        limits = c(0.5, 
                                                  max(label(data)$x) + 0.5))
        } else {
            p <- p + scale_x_continuous(breaks = data$labels$x, 
                                        labels = data$labels$label)
            p <- p + scale_y_continuous(expand = c(0,0))
            p <- p + theme(axis.text.x = element_text(angle = 45, 
                                                      hjust = 1, 
                                                      vjust = 1), 
                           axis.ticks.x = element_blank())
        }
        p
    }
)

### HELPER FUNCTIONS

#' Calculate bootstrap evolution data
#' 
#' This function bootstraps combinations of organisms at each size between 1 and
#' the length of the pangenome. For each combination the size of the different 
#' pangenome groups are recorded.
#' 
#' @param pangenome The pgVirtual subclass object
#' 
#' @param times The number of bootstrap samplings per size
#' 
#' @return A data.frame with org, Singleton, Accessory, Core and Total columns.
#' One row per sampling.
#' 
#' @noRd
#' 
evolBoot <- function(pangenome, times=10) {
    size <- length(pangenome)
    mat <- pgMatrix(pangenome)
    res <- lapply(seq(size), function(n) {
        res1 <- lapply(seq(times), function(i) {
            ind <- sample(size, size = n, replace = TRUE)
            subMat <- mat[, ind, drop = FALSE]
            data.frame(panGroups(subMat, defaults(pangenome)$coreThreshold), 
                       stringsAsFactors = FALSE)
        })
        data.frame(org = n, do.call(rbind, res1), stringsAsFactors = FALSE)
    })
    do.call(rbind, res)
}
#' Calculate evolution data based on manual ordering
#' 
#' This function takes an ordering of organism and calculates the sizes of the
#' different pangenome groups at each progression.
#' 
#' @param pangenome The pgVirtual subclass object
#' 
#' @param order An integer vector with the order of organism addition
#' 
#' @return A data.frame with org, Singleton, Accessory, Core and Total columns.
#' One row per element in order
#' 
#' @noRd
#' 
evolMan <- function(pangenome, order) {
    mat <- pgMatrix(pangenome)
    res <- lapply(seq(along = order), function(i) {
        ind <- order[1:i]
        subMat <- mat[, ind, drop = FALSE]
        data.frame(org = factor(orgNames(pangenome)[order[i]]), 
                   panGroups(subMat, defaults(pangenome)$coreThreshold), 
                   stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, res)
    levels(res$org) <- orgNames(pangenome)[order]
    res
}
#' Calculate pangenome group sizes
#' 
#' This function calculates the sizes of the different pangenome groups based on
#' a pangenome matrix.
#' 
#' @param mat A pangenome matrix with gene groups as rows and organisms as 
#' columns. Empty gene groups are allowed
#' 
#' @return A data.frame with one row and the columns Singleton, Accessory, Core
#' and Total.
#' 
#' @noRd
#' 
panGroups <- function(mat, coreThreshold = 1) {
    mat <- as(mat, 'nsparseMatrix')
    nGenes <- Matrix::rowSums(mat)
    data.frame(group = c('Singleton', 'Accessory', 'Core', 'Total'),
               size = c(sum(nGenes == 1),
                        sum(nGenes > 1 & nGenes / ncol(mat) < coreThreshold),
                        sum(nGenes / ncol(mat) >= coreThreshold),
                        sum(nGenes != 0)),
               stringsAsFactors = FALSE)
}
#' Calculate pangenomebased organism similarity
#' 
#' This function creates a similarity matrix for the organisms in a pangenome
#' based on the presence/absence pattern of their collective gene groups. The
#' similarity is defined as the number of gene groups they have in common,
#' divided by the total number of unique gene groups present in either.
#' 
#' @param pangenome The pgVirtual subclass object
#' 
#' @return A symmetric matrix with columns and rows for each organism.
#' 
#' @noRd
#' 
pgSim <- function(pangenome) {
    if (!hasGeneGroups(pangenome)) stop('Gene groups must be defined')
    mat <- pgMatrix(pangenome)
    panSim(mat@p, mat@i, colnames(mat))
}
#' Calculate pangenome-based organism distance
#' 
#' This function calculates the distance between organisms based on the 
#' pangenome matrix. It simply calls dist on the pangenome matrix.
#' 
#' @param pangenome The pgVirtual subclass object
#' 
#' @param method Passed on to dist()
#' 
#' @return A distance matrix
#' 
#' @importFrom stats dist
#' 
#' @noRd
#' 
pgDist <- function(pangenome, method) {
    if (!hasGeneGroups(pangenome)) stop('Gene groups must be defined')
    dist(t(pgMatrix(pangenome)), method = method)
}
#' Calculate percent distribution of gene groups
#' 
#' This function simply returns the percentage of Core, Accessory and Singleton
#' gene groups in a pangenome.
#' 
#' @param pangenome A pgVirtual subclass
#' 
#' @return A named numeric vector
#' 
#' @noRd
#' 
groupDistribution <- function(pangenome) {
    ans <- c(Core = 0, Accessory = 0, Singleton = 0)
    gDist <- table(groupInfo(pangenome)$group)
    gDist <- gDist/sum(gDist) * 100
    ans[names(gDist)] <- gDist
    ans
}
#' Create a horizontal progress bar scaled to 50
#' 
#' This function create a string of '=' and possibly ':' (in case of uneven
#' numbers)
#' 
#' @param percent The percentage to create the progress bar for
#' 
#' @return A character string
#' 
#' @noRd
#' 
createProgress <- function(percent) {
    equals <- round(percent)/2
    half <- round(percent) %% 2 != 0
    paste(c(rep('=', equals), if (half) ':' else ''), collapse = '')
}