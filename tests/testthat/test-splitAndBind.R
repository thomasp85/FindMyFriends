context("Custom split and bind functions")

matrixA <- matrix(1:4, ncol=2, dimnames=list(NULL, c('a', 'b')))
matrixB <- matrix(5:8, ncol=2, dimnames=list(NULL, c('c', 'a')))
matrixC <- matrix(9:14, ncol=3, dimnames=list(NULL, c('c', 'd', 'e')))
pg <- .loadPgExample()
geneSeqs <- genes(pg, subset=1:10)

test_that("rbindMat works", {
    expect_equal(rbindMat(matrixA, matrixB, matrixC), rbindMat(list(matrixA, matrixB, matrixC)))
    expect_equal(colnames(rbindMat(matrixA, matrixB)), c('a', 'b', 'c'))
    expect_true(all(is.na(rbindMat(matrixA, matrixB)[1:2, 'c'])))
    expect_true(all(rbindMat(matrixA, matrixB, fill=0)[1:2, 'c']==0))
})

test_that("splitStringSet works", {
    splitted <- splitStringSet(geneSeqs, rep(1:5, 2))
    expect_is(splitted, 'AAStringSetList')
    expect_equal(length(splitted), 5)
    expect_equal(geneSeqs[6], splitted[[1]][2])
})

test_that("paralogueInd works", {
    expect_equal(paralogueInd(rep(1:5, 2), rep(1:2, length.out=5)), c(1, 2, 1, 2, 1, 1, 2, 1, 2, 1))
})