context("Gene group splitting")

pg1 <- .loadPgExample(geneLoc = TRUE, withGroups = TRUE)

neighbors <- getNeighbors(pg1, zeroInd = TRUE)

test_that("neighborhoodSplit works", {
    skip_on_os('win')
    set.seed(1)
    pg2 <- neighborhoodSplit(pg1, lowerLimit=0.8)
    expect_equal(nGeneGroups(pg2), 3136)
    expect_equal(sum(seqToGeneGroup(pg2)), 8794625)
})

test_that("collectNeighbors works", {
    expect_equal(collectNeighbors(1:5, 'f', 4), c("2;3;4;5", "3;4;5;NA", "4;5;NA;NA", "5;NA;NA;NA", "NA;NA;NA;NA"))
    expect_equal(collectNeighbors(1:5, 'b', 4), c("NA;NA;NA;NA", "1;NA;NA;NA", "2;1;NA;NA", "3;2;1;NA", "4;3;2;1"))
})

test_that("neighborhoodSim works", {
    grouping <- seqToGeneGroup(pg1)
    
    NS1 <- neighborhoodSim(which(grouping == 1), grouping, seqToOrg(pg1), 4, 
                           neighbors$down, neighbors$up, neighbors$reverse, geneWidth(pg1), 0.1, forceParalogues=T)
    expect_is(NS1, 'list')
    expect_named(NS1, c('i', 'p', 'x'))
    expect_equal(length(NS1$i), 31)
    expect_equal(length(NS1$p), 21)
    expect_equal(length(NS1$x), 31)
    expect_equal(sum(unlist(NS1)), 990)
})

test_that("neighborSplitting works", {
    newGroups1 <- neighborSplitting(1, pg1, seqToOrg(pg1), neighbors, seqToGeneGroup(pg1), geneWidth(pg1), 0.1, TRUE, 4, 5, 0.75,  rep(1L, nGenes(pg1)))
    newGroups2 <- neighborSplitting(1, pg1, seqToOrg(pg1), neighbors, seqToGeneGroup(pg1), geneWidth(pg1), 0.1, FALSE, 4, 5, 0.75,  rep(1L, nGenes(pg1)))
    newGroups3 <- neighborSplitting(1, pg1, seqToOrg(pg1), neighbors, seqToGeneGroup(pg1), geneWidth(pg1), 0.1, TRUE, 7, 5, 0.75,  rep(1L, nGenes(pg1)))
    expect_is(newGroups1, 'list')
    expect_equal(length(newGroups1), 10)
    expect_equal(unname(sapply(newGroups1, sum)), c(3803L, 5808L, 5911L, 2492L, 2652L, 2720L, 2796L, 3685L, 3336L, 
                                            2677L))
    expect_is(newGroups2, 'list')
    expect_equal(length(newGroups2), 2)
    expect_equal(unname(sapply(newGroups2, sum)), c(3803L, 32077L))
    expect_is(newGroups3, 'list')
    expect_equal(length(newGroups3), 10)
    expect_equal(unname(sapply(newGroups3, sum)), c(5808L, 6402L, 2492L, 3803L, 2317L, 2662L, 2748L, 3635L, 3336L, 
                                            2677L))
})

test_that("Clique extraction works", {
    edges <- data.frame(from = c(1L, 2L, 3L, 5L, 6L, 6L, 7L, 7L, 8L), 
                        to = c(9L, 10L, 11L, 12L, 13L, 14L, 13L, 14L, 15L), 
                        nSim = c(8L, 8L, 8L, 8L, 8L, 2L, 2L, 8L, 8L), 
                        sSim = c(0.955752212389382, 1, 0.977375565610857, 1, 1, 0.756839003805933, 0.756839003805933, 0.999999999999998, 0.999999999999998), 
                        gSim = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L))
    members <- 15
    cliques <- extractCliques(edges, members)
    expect_is(cliques, 'integer')
    expect_equal(cliques, c(7L, 1L, 6L, 8L, 2L, 3L, 4L, 5L, 7L, 1L, 6L, 2L, 3L, 4L, 5L))
})

test_that("anyParalogues works", {
    expect_equal(sum(anyParalogues(pg1)), 93)
})
