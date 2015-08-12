context("Gene grouping")

pg <- .loadPgExample()

test_that("Manual grouping works", {
    fauxGroups <- rep(1:10, length.out=nGenes(pg))
    fauxGroups2 <- split(1:length(fauxGroups), fauxGroups)
    expect_equal(manualGrouping(pg, fauxGroups), manualGrouping(pg, fauxGroups2))
    expect_equal(nGeneGroups(manualGrouping(pg, fauxGroups)), 10)
    expect_equal(seqToGeneGroup(manualGrouping(pg, fauxGroups)), fauxGroups)
})

test_that("Graph grouping works", {
    set.seed(1)
    gr <- erdos.renyi.game(100, 0.015)
    members <- igMembers(gr, 'infomap')
    expect_equal(100, length(members))
    expect_is(members, 'numeric')
    expect_equal(members[c(1, 25, 50, 100)], c(7, 3, 2, 5))
    expect_equal(members, igMembers(gr, 'cluster_infomap'))
})