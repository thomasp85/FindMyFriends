context("Similarity calculations work")

pg <- .loadPgExample()

test_that("kmerSimilarity works", {
    skip_on_os('win')
    sim1 <- kmerSimilarity(pg, lowMem=T, kmerSize=4, lowerLimit=0.8)
    sim2 <- kmerSimilarity(pg, lowMem=F, kmerSize=4, lowerLimit=0.8)
    sim3 <- kmerSimilarity(pg, lowMem=F, kmerSize=4, lowerLimit=0.8, pParam=BiocParallel::SerialParam())
    expect_equal(sim1, sim2)
    expect_equal(sim2, sim3)
    expect_is(sim2, 'dgCMatrix')
    expect_equal(dim(sim2), rep(nGenes(pg), 2))
    expect_equal(sum(sim2), 4448.01858751133)
})