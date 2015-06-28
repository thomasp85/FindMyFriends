context("Parallel versions of linearKernel")

pg <- .loadPgExample()

test_that("Split-combine works", {
    chunks <- getChunks(9, 3)
    expect_named(chunks, c('combs', 'chunks'))
    expect_named(chunks$combs, c('col', 'row', 'origInd'))
    expect_equal(chunks$combs$col, c(1,1,2,1,2,3))
    expect_equal(chunks$combs$row, c(2,3,3,1,2,3))
    expect_equal(chunks$combs$origInd, 1:6)
    expect_equal(chunks$chunks[,'from'], c(1,4,7))
    expect_equal(chunks$chunks[,'to'], c(3,6,9))
    result <- Matrix::Matrix(1:81, ncol=9, nrow=9, sparse=TRUE)
    result[upper.tri(result)] <- 0
    splits <- lapply(1:nrow(chunks$combs), function(i) {
        result[chunks$chunks[chunks$combs$row[i], 1]:chunks$chunks[chunks$combs$row[i], 2],
               chunks$chunks[chunks$combs$col[i], 1]:chunks$chunks[chunks$combs$col[i], 2]]
    })
    expect_equal(result, weaveChunks(splits, chunks))
})

test_that("linearKernel can be run in parallel", {
    pParam <- BiocParallel::SerialParam()
    er <- kebabs::getExRep(genes(pg, subset=1:9))
    expect_equal(lkParallel(er, pParam, 3), kebabs::linearKernel(er, sparse=T, diag=F))
})

test_that("Low memory mode works", {
    pParam <- BiocParallel::SerialParam()
    pg <- removeGene(pg, organism = 2:length(pg))
    res <- lkParallelLM(pg, 3, pParam, 3)
    res2 <- kebabs::linearKernel(kebabs::getExRep(genes(pg)), sparse=T, diag=F)
    expect_equal(res, res2)
})