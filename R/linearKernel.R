#' @importFrom Matrix sparseMatrix
#' @export
#' 
lkFMF <- function(x, selx, order, lowerLimit = 0, upperLimit = 1) {
    
    if (missing(selx)) {
        selx <- seq_len(dim(x)[1])
    }
    if (!missing(order)) {
        selx <- selx[order]
    }
    lk <- linearKernel(x@p, x@j, x@x, selx - 1, lowerLimit, upperLimit)
    lk$member <- lk$member + 1
    if (!missing(order)) {
        lk$member[order] <- lk$member
    }
    lkSize <- max(lk$member)
    sm <- sparseMatrix(i = lk$i, p = lk$p, x = lk$x, dims = c(lkSize, lkSize), 
                       symmetric = TRUE, index1 = FALSE)[lk$member, lk$member]
    dimnames(sm) <- dimnames(x)[c(1, 1)]
    sm
}
#' @importFrom Matrix t
#' @export
clustersFromAdjMatrix <- function(x) {
    if (x@uplo == 'U') x <- t(x)
    FindMyFriends:::getClusters(dim(x)[1], x@i, x@p, x@x)
}