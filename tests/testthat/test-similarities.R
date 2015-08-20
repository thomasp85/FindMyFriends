context("Similarity calculations work")

pg <- .loadPgExample()

test_that("kmerSimilarity works", {
    sim1 <- kmerSimilarity(pg, lowMem=T, kmerSize=4, lowerLimit=0.8)
    sim2 <- kmerSimilarity(pg, lowMem=F, kmerSize=4, lowerLimit=0.8)
    sim3 <- kmerSimilarity(pg, lowMem=F, kmerSize=4, lowerLimit=0.8, pParam=BiocParallel::SerialParam(), nSplits=4)
    expect_equal(sim1, sim2)
    expect_equal(sim1, sim3)
    expect_is(sim1, 'dgCMatrix')
    expect_equal(dim(sim1), rep(nGenes(pg), 2))
    expect_equal(sum(sim1), 4448.01858751133)
})