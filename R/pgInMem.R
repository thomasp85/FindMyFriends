#' @include pgVirtual.R
NULL

#' FindMyFriends standard base class for pangenomic data
#' 
#' This virtual class is the superclass of the standard pangenome classes in
#' FindMyFriends. It defines storage for everything except gene information,
#' which is delegated to its subclasses.
#' 
#' As gene storage is not defined in this class the following methods must be
#' defined by subclasses:
#' 
#' \describe{
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
#' }
#' 
#' @slot seqToOrg An integer vector that reference all genes to a specific 
#' organism.
#' 
#' @slot seqToGeneGroup An integer vector that references all genes to a 
#' specific gene group.
#' 
#' @slot groupInfo A data.frame storing metadata information about gene groups.
#' 
#' @slot orgInfo A data.frame storing metadata information about organisms
#' 
#' @family Pangenome_classes
#' 
#' @export
#' 
setClass(
    'pgInMem',
    contains = c('VIRTUAL', 'pgVirtual'),
    slots = list(
        seqToOrg = 'integer',
        seqToGeneGroup = 'integer',
        groupInfo = 'data.frame',
        orgInfo = 'data.frame'
    ),
    validity = function(object) {
        if (length(object@seqToGeneGroup) != 0) {
            if (length(object@seqToOrg) != length(object@seqToGeneGroup)) {
                return('Gene indexes of different length')
            }
        }
        return(TRUE)
    },
    prototype = list(
        seqToOrg = integer(),
        seqToGeneGroup = integer(),
        groupInfo = data.frame(),
        orgInfo = data.frame()
    )
)

### REQUIRED

# Methods genes, geneNames, geneNames<- and geneWidth are passed on to subclasses

#' @describeIn seqToOrg Gene to organism indexing for pgInMem subclasses
#' 
setMethod(
    'seqToOrg', 'pgInMem',
    function(object) {
        object@seqToOrg
    }
)
#' @describeIn seqToGeneGroup Gene to genegroup indexing for pgInMem subclasses
#' 
setMethod(
    'seqToGeneGroup', 'pgInMem',
    function(object) {
        object@seqToGeneGroup
    }
)
#' @describeIn orgNames Get organism names for pgInMem subclasses
#' 
setMethod(
    'orgNames', 'pgInMem',
    function(object) {
        rownames(object@orgInfo)
    }
)
#' @describeIn orgNames Set organism names for pgInMem subclasses
#' 
setMethod(
    'orgNames<-', 'pgInMem',
    function(object, value) {
        rownames(object@orgInfo) <- value
        object
    }
)
#' @describeIn groupNames Get gene group names for pgInMem subclasses
#' 
setMethod(
    'groupNames', 'pgInMem',
    function(object) {
        rownames(object@groupInfo)
    }
)
#' @describeIn groupNames Set gene group names for pgInMem subclasses
#' 
setMethod(
    'groupNames<-', 'pgInMem',
    function(object, value) {
        rownames(object@groupInfo) <- value
        object
    }
)
#' @describeIn orgInfo Get organism metadata for pgInMem subclasses
#' 
setMethod(
    'orgInfo', 'pgInMem',
    function(object) {
        object@orgInfo
    }
)
#' @describeIn orgInfo Set organism metadata for pgInMem subclasses
#' 
setMethod(
    'orgInfo<-', 'pgInMem',
    function(object, value) {
        object@orgInfo <- value
        object
    }
)
#' @describeIn groupInfo Get gene group metadata for pgInMem subclasses
#' 
setMethod(
    'groupInfo', 'pgInMem',
    function(object) {
        object@groupInfo
    }
)
#' @describeIn groupInfo Set gene group metadata for pgInMem subclasses
#' 
setMethod(
    'groupInfo<-', 'pgInMem',
    function(object, value) {
        object@groupInfo <- value
        object
    }
)
#' @rdname internalGroupGenes
#' 
setMethod(
    'groupGenes', 'pgInMem',
    function(object, seqToGeneGroup) {
        object <- callNextMethod()
        object@seqToGeneGroup <- seqToGeneGroup
        object
    }
)
#' @describeIn removeGene Gene removal base function for pgInMem subclasses
#' 
setMethod(
    'removeGene', c('pgInMem', 'missing', 'missing', 'missing', 'numeric'),
    function(object, name, organism, group, ind) {
        currentOrgs <- unique(seqToOrg(object))
        object@seqToOrg <- object@seqToOrg[-ind]
        removedOrgs <- currentOrgs[!currentOrgs %in% unique(object@seqToOrg)]
        if (hasGeneGroups(object)) {
            currentGroups <- unique(seqToGeneGroup(object))
            object@seqToGeneGroup <- object@seqToGeneGroup[-ind]
            keptGroups <- currentGroups %in% unique(object@seqToGeneGroup)
            removedGroups <- currentGroups[!keptGroups]
            
            if (length(removedGroups) != 0) {
                object@seqToGeneGroup <- removeIndex(object@seqToGeneGroup)
                groupInfo(object) <- groupInfo(object)[-removedGroups, , 
                                                       drop = FALSE]
            }
        }
        if (length(removedOrgs) != 0) {
            object@seqToOrg <- removeIndex(object@seqToOrg)
            orgInfo(object) <- orgInfo(object)[-removedOrgs, , drop = FALSE]
        }
        if (inherits(object, 'pgInMemLoc')) {
            object@geneLocation <- object@geneLocation[-ind, , drop = FALSE]
        }
        if (inherits(object, 'pgFull')) {
            object@sequences <- object@sequences[-ind]
        } else if (inherits(object, 'pgLM')) {
            object@seqIndex <- object@seqIndex[-ind, , drop = FALSE]
        }
        orgInfo(object)$nGenes <- 0
        nGenesOrg <- lengths(split(seq_len(nGenes(object)), seqToOrg(object)))
        orgInfo(object)$nGenes[as.integer(names(nGenesOrg))] <- nGenesOrg
        object
    }
)
#' @rdname internalMetadata
#' 
setMethod(
    'setGroupInfo', 'pgInMem',
    function(object, name, info, key) {
        if (is.factor(info)) info <- as.character(info)
        if (is.null(object@groupInfo[[name]])) {
            object@groupInfo[[name]] <- NA
        } else if (is.factor(object@groupInfo[[name]])) {
            object@groupInfo[[name]] <- as.character(object@groupInfo[[name]])
        }
        object@groupInfo[[name]][key] <- info
        object
    }
)
#' @rdname internalMetadata
#' 
setMethod(
    'setOrgInfo', 'pgInMem',
    function(object, name, info, key) {
        if (is.factor(info)) info <- as.character(info)
        if (is.null(object@orgInfo[[name]])) {
            object@orgInfo[[name]] <- NA
        } else if (is.factor(object@orgInfo[[name]])) {
            object@orgInfo[[name]] <- as.character(object@orgInfo[[name]])
        }
        object@orgInfo[[name]][key] <- info
        object
    }
)
#' @rdname internalMergePangenomes
#' 
#' @importFrom dplyr bind_rows
#' 
setMethod(
    'mergePangenomes', c('pgInMem', 'pgInMem'),
    function(pg1, pg2, geneGrouping, groupInfo, ...) {
        if (class(pg1) != class(pg2)) {
            stop('pangenomes must be instances of the same class')
        }
        seqToOrg <- c(pg1@seqToOrg, nOrganisms(pg1) + pg2@seqToOrg)
        orgInfo <- as.data.frame(bind_rows(pg1@orgInfo, pg2@orgInfo))
        new(
            class(pg1),
            .settings = pg1@.settings,
            seqToOrg = as.integer(seqToOrg),
            seqToGeneGroup = as.integer(geneGrouping),
            orgInfo = orgInfo,
            groupInfo = groupInfo,
            ...
        )
    }
)
