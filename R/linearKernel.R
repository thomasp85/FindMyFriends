#' @importFrom Matrix sparseMatrix
#' 
lkFMFmat <- function(x, selx, order, lowerLimit = 0, upperLimit = 1) {
    if (missing(selx)) {
        selx <- seq_len(dim(x)[1])
    }
    if (!missing(order)) {
        selx <- selx[order]
    }
    lk <- lkMatrix(x@p, x@j, x@x, selx - 1, lowerLimit, upperLimit)
    lk$member <- lk$member + 1
    if (!missing(order)) {
        lk$member[order] <- lk$member
    }
    lkSize <- max(lk$member)
    sm <- sparseMatrix(j = lk$i, p = lk$p, x = lk$x, 
                       dims = c(lkSize, lkSize), symmetric = TRUE, 
                       index1 = FALSE)[lk$member, lk$member]
    dimnames(sm) <- dimnames(x)[c(1, 1)]
    sm
}
lkFMF <- function(x, selx, order, lowerLimit = 0, upperLimit = 1) {
    
    if (missing(selx)) {
        selx <- seq_len(dim(x)[1])
    }
    if (!missing(order)) {
        selx <- selx[order]
    }
    members <- lkMembers(x@p, x@j, x@x, selx - 1, lowerLimit, upperLimit)
    if (!missing(order)) {
        members[order] <- members
    }
    members
}
#' @importFrom Matrix t
#' 
clustersFromAdjMatrix <- function(x) {
    if (x@uplo == 'U') x <- t(x)
    getClusters(x@i, x@p, x@x)
}