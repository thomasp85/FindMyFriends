context("Gene grouping")

pg <- .loadPgExample()

test_that("Manual grouping works", {
    fauxGroups <- rep(1:10, length.out=nGenes(pg))
    fauxGroups2 <- split(seq_along(fauxGroups), fauxGroups)
    expect_equal(manualGrouping(pg, fauxGroups), manualGrouping(pg, fauxGroups2))
    expect_equal(nGeneGroups(manualGrouping(pg, fauxGroups)), 10)
    expect_equal(seqToGeneGroup(manualGrouping(pg, fauxGroups)), fauxGroups)
})

test_that("Graph grouping works", {
    set.seed(1)
    gr <- igraph::erdos.renyi.game(100, 0.015)
    set.seed(1)
    members <- igMembers(gr, 'infomap')
    expect_equal(100, length(members))
    expect_is(members, 'membership')
    expect_equal(members[c(1, 25, 50, 100)], c(7, 3, 2, 5))
    set.seed(1)
    expect_equal(members, igMembers(gr, 'cluster_infomap'))
})