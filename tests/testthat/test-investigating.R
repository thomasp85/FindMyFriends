context("Result investigation")

pg <- .loadPgExample(geneLoc = TRUE, withGroups = TRUE, withNeighborhoodSplit = TRUE)

test_that("Groupstatistics works", {
    stats1 <- groupStat(pg, 1)
    stats2 <- groupStat(pg, 2)
    expect_is(stats1, 'list')
    expect_is(stats2, 'list')
    expect_equal(length(stats1), length(stats2))
    expect_equal(length(stats1), nGeneGroups(pg))
    expect_named(stats1[[1]], c('maxOrg', 'minLength', 'maxLength', 'sdLength', 'genes', 'backward', 'forward'))
    expect_equal(length(unlist(strsplit(stats1[[1]]$backward, ';'))), length(stats1[[1]]$genes))
    expect_equal(length(unlist(strsplit(stats2[[1]]$backward, ';'))), 2 * length(stats2[[1]]$genes))
})

test_that("Organism statistics work", {
    stats1 <- orgStat(pg)
    stats2 <- orgStat(pg, 1)
    stats3 <- orgStat(pg, orgNames(pg)[1])
    stats4 <- orgStat(pg, 1:3)
    stats5 <- orgStat(pg, getFrequency=TRUE)
    expect_is(stats1, 'data.frame')
    expect_equal(stats2, stats3)
    expect_equal(stats1[1:3,], stats4)
    expect_named(stats1, c('nGenes', 'minLength', 'maxLength', 'sdLength', 'nGeneGroups', 'nParalogues'))
    expect_equal(dim(stats5), c(5, 37))
    expect_equal(as.numeric(stats2), c(1336, 30, 1824, 121.171824254577, 1336, 0))
})

test_that("Panchromosome works", {
    pc <- pcGraph(pg)
    expect_is(pc, 'igraph')
    expect_equal(igraph::gorder(pc), nGeneGroups(pg))
    expect_equal(igraph::gsize(pc), 3529)
    expect_named(igraph::vertex.attributes(pc), c('name', 'nMembers', 'organisms', 'strands'))
    expect_named(igraph::edge.attributes(pc), c('weight', 'organisms'))
})

# test_that("Variable regions detection works", {
#     vr <- variableRegions(pg)
#     expect_is(vr, 'list')
#     expect_named(vr[[1]], c('type', 'members', 'flank', 'connectsTo', 'graph'))
#     expect_equal(length(vr), 232)
# })

test_that("Neighborhood extraction works", {
    neigh <- getNeighborhood(pg, 4, 4)
    expect_is(neigh, 'igraph')
    expect_equal(igraph::gorder(neigh), 11)
    expect_equal(igraph::gsize(neigh), 12)
    expect_named(igraph::vertex.attributes(neigh), c('name', 'centerGroup'))
    expect_named(igraph::edge.attributes(neigh), c('weight'))
})

test_that("Group trailing works", {
    tr1 <- trailGroups(10, pg, 4)
    tr2 <- trailGroups(10:12, pg, 4)
    tr3 <- trailGroups(10, pg, 7)
    expect_is(tr1, 'list')
    expect_is(tr2, 'list')
    expect_is(tr3, 'list')
    expect_equal(tr1[[1]], tr2[[1]])
    expect_equal(tr1[[1]], c(814,  565,  813,  564,  563,  562, 1037, 1367, 1237))
    expect_true(all(tr1[[1]] %in% tr3[[1]]))
    expect_equal(length(tr3[[1]]), 7*2+1)
    expect_equal(length(tr2), 3)
})

test_that("Conversion of trails to graphs work", {
    fakeTrail <- list(
        c(1,2,3,4,5),
        c(1,2,3,4,5),
        c(1,2,6,4,5),
        c(1,2,6,4,7)
    )
    trGraph <- trailsToGraph(fakeTrail)
    trueGraph <- data.frame(from = c("1", "2", "2", "3", "4", "4", "6"),
                            to = c("2", "3", "6", "4", "5", "7", "4"),
                            weight = c(4L, 2L, 2L, 2L, 3L, 1L, 2L),
                            stringsAsFactors = FALSE)
    expect_equal(igraph::as_data_frame(trGraph), trueGraph)
})

test_that("Simple scaling works", {
    r1 <- scaleRange(rep(1, 5), 1, 5)
    r2 <- scaleRange(1:5, 10, 50)
    expect_equal(r1, rep(3, 5))
    expect_equal(r2, seq(10, 50, by=10))
})

# test_that("Cycle locating works", {
#     pc <- pcGraph(pg)
#     cycles <- locateCycles(pc, maxLength=4)
#     expect_is(cycles, 'list')
#     expect_equal(length(cycles), 568)
#     expect_is(unlist(cycles), 'character')
#     expect_equal(length(unlist(cycles)), 3010)
#     merged <- mergeCycles(cycles)
#     expect_is(merged, 'list')
#     expect_equal(length(merged), 232)
#     expect_is(unlist(merged), 'character')
#     expect_equal(length(unlist(merged)), 1113)
#     summaryCycles <- summarizeCycles(merged, pc)
#     expect_is(summaryCycles, 'list')
#     expect_equal(length(summaryCycles), length(merged))
# })

test_that("Route from BFS can be extracted", {
    route <- c(2, 3, 4, 5, NA, 5, 6, 7, 8, 1)
    expect_equal(getRoute(1, route), 5:1)
    expect_equal(getRoute(5, route), 5)
})