## TBD
# Include gpcGrouping, recurseCompare, recurseCompPar, fillTree, cutK

context("GPC and helpers")

pg <- .loadPgExample()

test_that("gpcGrouping works", {
    set.seed(1)
    pg1 <- gpcGrouping(pg)
    set.seed(1)
    pg2 <- gpcGrouping(pg, lowMem=T)
    cacheDir <- tempdir()
    set.seed(1)
    pg3 <- gpcGrouping(pg, cacheDB=cacheDir)
    pg4 <- gpcGrouping(pg, cacheDB=cacheDir)
    expect_equal(pg1, pg2)
    expect_equal(pg2, pg3)
    expect_equal(pg3, pg4)
})
