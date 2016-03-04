context("Organism comparison")

pg <- .loadPgExample(withGroups = TRUE)
pgNoGroup <- .loadPgExample()

test_that("Organism exrep works", {
    orgER <- orgExRep(pg, 3)
    expect_is(orgER, 'ExplicitRepresentationDense')
    expect_equal(as.numeric(orgER[1,1]), 0.0119486627166315)
    expect_equal(as.numeric(orgExRep(pg, 3, pParam=BiocParallel::SerialParam())[1,1]), 0.0119486627166315)
})
test_that("Organism distance works", {
    expect_equal(kmerDist(pg, 3, method='cosine')[1], 0.0393149129234893)
    expect_equal(kmerDist(pg, 3, method='euclidean')[1], 0.0555996830599097)
    expect_equal(pgDist(pg, 'euclidean')[1], 16.9410743460974)
    expect_error(pgDist(pgNoGroup, 'euclidean'))
})
test_that("Organism similarity works", {
    expect_equal(kmerSim(pg, 3)[2], 0.998454337621818)
    expect_equal(pgSim(pg)[2], 0.849137931034483)
    expect_error(pgSim(pgNoGroup))
})