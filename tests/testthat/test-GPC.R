## TBD
# Include gpcGrouping, recurseCompare, recurseCompPar, fillTree, cutK

context("GPC and helpers")

pg <- .loadPgExample()

test_that("gpcGrouping works", {
    skip_on_os('win')
    set.seed(1)
    pg1 <- gpcGrouping(pg, kmerSize = 4)
    set.seed(1)
    pg2 <- gpcGrouping(pg, lowMem=T, kmerSize = 4)
    cacheDir <- tempdir()
    set.seed(1)
    pg3 <- gpcGrouping(pg, cacheDB=cacheDir, kmerSize = 4)
    pg4 <- gpcGrouping(pg, cacheDB=cacheDir, kmerSize = 4)
    expect_equal(pg1, pg2)
    expect_equal(pg2, pg3)
    expect_equal(pg3, pg4)
})
