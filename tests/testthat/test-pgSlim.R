context("pgSlim class methods")

pg1 <- .loadPgExample(withGroups = T)
pg2 <- .loadPgExample(withGroups = T, geneLoc = T)

test_that("Coercion works", {
    expect_is(as(pg1, 'pgSlim'), 'pgSlim')
    expect_is(as(pg2, 'pgSlim'), 'pgSlim')
    expect_is(as(pg2, 'pgSlimLoc'), 'pgSlimLoc')
    expect_error(as(pg1, 'pgSlimLoc'))
})

test_that("genes setters and getters are blocked", {
    pg1 <- as(pg1, 'pgSlim')
    expect_error(genes(pg1))
    expect_error(genes(pg1, split = 'group'))
    expect_error(geneNames(pg1))
    expect_error(geneNames(pg1) <- 1:nGenes(pg1))
    expect_error(geneWidth(pg1))
})