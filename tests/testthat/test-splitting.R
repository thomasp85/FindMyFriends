context("Gene group splitting")

pg1 <- .loadPgExample(geneLoc = TRUE, withGroups = TRUE)

geneGroup <- list(genes = c(567L, 594L, 809L, 929L, 987L, 1035L, 1055L, 1085L, 1229L, 1552L, 1723L, 1733L, 1761L, 1925L, 2165L, 2281L, 2574L, 2600L, 2677L, 2834L, 3765L), 
                  organism = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L), 
                  down = list(`567` = c(1081L, 395L, 396L, 2L), 
                              `594` = c(1284L, 788L, 789L, 2L), 
                              `809` = c(306L, 305L, 747L, 2L), 
                              `929` = c(256L, 730L, 1019L, 2L), 
                              `987` = c(719L, 718L, 1005L, 2L), 
                              `1035` = c(212L, 704L, 213L, 2L), 
                              `1055` = c(201L, 1244L, 987L, 2L), 
                              `1085` = c(23L, 20L, 29L, 2L), 
                              `1229` = c(955L, 139L, 138L, 954L), 
                              `1552` = c(1305L, 875L, 442L, 2L), 
                              `1723` = c(522L, 1106L, 523L, 2L), 
                              `1733` = c(837L, 838L, 1105L, 2L), 
                              `1761` = c(539L, 540L, 829L, 2L), 
                              `1925` = c(1081L, 2714L, 1309L, 2L), 
                              `2165` = c(305L, 747L, 2671L, 2L), 
                              `2281` = c(729L, 17L, 17L, 2L), 
                              `2574` = c(955L, 139L, 138L, 954L), 
                              `2600` = c(121L, 1216L, 122L, 2L), 
                              `2677` = c(660L, 92L, 928L, 2L), 
                              `2834` = c(924L, 1185L, 598L, 2L), 
                              `3765` = c(19L, 24L, 2L, 2L)), 
                  up = list(`567` = c(81L, 794L, 397L, 398L), 
                            `594` = c(790L, 1079L, 1080L, 389L), 
                            `809` = c(1049L, 1253L, 304L, 2628L), 
                            `929` = c(17L, 1249L, 729L, 1018L), 
                            `987` = c(717L, 229L, 2620L, 716L), 
                            `1035` = c(705L, 214L, 215L, 1221L), 
                            `1055` = c(988L, 202L, 989L, 703L), 
                            `1085` = c(24L, 19L, 21L, 23L), 
                            `1229` = c(681L, 137L, 136L, 680L), 
                            `1552` = c(876L, 441L, 440L, 439L), 
                            `1723` = c(2552L, 2736L, 524L, 525L), 
                            `1733` = c(1187L, 528L, 527L, 526L), 
                            `1761` = c(541L, 828L, 542L, 543L), 
                            `1925` = c(11L, 42L, 1235L, 51L), 
                            `2165` = c(1049L, 304L, 2670L, 1222L), 
                            `2281` = c(1019L, 730L, 256L, 257L), 
                            `2574` = c(681L, 137L, 136L, 680L), 
                            `2600` = c(948L, 123L, 124L, 125L), 
                            `2677` = c(93L, 661L, 929L, 662L), 
                            `2834` = c(886L, 1152L, 887L, 1153L), 
                            `3765` = c(29L, 20L, 41L, 23L))
)

test_that("neighborhoodSplit works", {
    skip_on_os('win')
    set.seed(1)
    pg2 <- neighborhoodSplit(pg1, lowerLimit=0.8)
    expect_equal(nGeneGroups(pg2), 3158)
    expect_equal(sum(seqToGeneGroup(pg2)), 9740008)
})

test_that("collectNeighbors works", {
    expect_equal(collectNeighbors(1:5, 'f', 4), c("2;3;4;5", "3;4;5;NA", "4;5;NA;NA", "5;NA;NA;NA", "NA;NA;NA;NA"))
    expect_equal(collectNeighbors(1:5, 'b', 4), c("NA;NA;NA;NA", "1;NA;NA;NA", "2;1;NA;NA", "3;2;1;NA", "4;3;2;1"))
})

test_that("neighborhoodSimilarity works", {
    NS1 <- neighborhoodSimilarity(geneGroup, minFlank=1, forceParalogues=T)
    NS2 <- neighborhoodSimilarity(geneGroup, minFlank=0, forceParalogues=T)
    NS3 <- neighborhoodSimilarity(geneGroup, minFlank=0, forceParalogues=F)
    expect_equal(dim(NS1), c(21, 21))
    expect_is(NS1, 'matrix')
    expect_equal(sum(NS1), 15)
    expect_equal(sum(NS2), 120)
    expect_equal(sum(NS3), 185)
})

test_that("neighborSplitting works", {
    newGroups1 <- neighborSplitting(geneGroup, pangenome=pg1, kmerSize=4, lowerLimit=0.8, minFlank=1, forceParalogues=T, maxLengthDif=0.1)
    newGroups2 <- neighborSplitting(geneGroup, pangenome=pg1, kmerSize=4, lowerLimit=0.8, minFlank=0, forceParalogues=T, maxLengthDif=0.1)
    newGroups3 <- neighborSplitting(geneGroup, pangenome=pg1, kmerSize=4, lowerLimit=0.8, minFlank=0, forceParalogues=F, maxLengthDif=0.1)
    expect_is(newGroups1, 'list')
    expect_equal(length(newGroups1), 18)
    expect_equal(sapply(newGroups1, sum), c(567, 594, 2974, 929, 987, 1035, 1055, 4850, 3803, 1552, 1723, 1733, 1761, 1925, 2281, 2600, 2677, 2834))
    expect_is(newGroups2, 'list')
    expect_equal(length(newGroups2), 10)
    expect_equal(sapply(newGroups2, sum), c(6739, 6596, 3635, 2912, 2690, 2327, 3336, 2290, 1552, 3803))
    expect_is(newGroups3, 'list')
    expect_equal(length(newGroups3), 2)
    expect_equal(sapply(newGroups3, sum), c(32077, 3803))
})

test_that("Clique extraction works", {
    gr <- data.frame(
        from = c("1089", "1081", "1089", "1089", "1081", "1089", "2437", "2444", "1089", "2437", "2444", "2451"), 
        to = c("2437", "2444", "2444", "2451", "3761", "3761", "3761", "3761", "3771", "3771", "3771", "3771"), 
        nWeight = c(3, 3, 4, 5, 7, 3, 3, 4, 6, 3, 4, 4), 
        sWeight = c(0.866071428571428, 0.807183003750947, 0.905838704209395, 0.892857142857142, 0.928571428571428, 0.901785714285713, 0.928571428571428, 0.843057803917655, 0.892857142857142, 0.821428571428571, 0.834089103875978, 0.999999999999999)
    )
    gr <- igraph::graph_from_data_frame(gr, directed=FALSE)
    clique <- extractClique(gr, 5)
    expect_is(clique, 'character')
    expect_equal(clique, c("1089", "3771", "2451"))
})

test_that("trailGroups2 works", {
    trail1 <- trailGroups2(10, seqToGeneGroup(pg1), pg1, 4)
    trail2 <- trailGroups2(c(10, 21), seqToGeneGroup(pg1), pg1, 4)
    trail3 <- trailGroups2(10, seqToGeneGroup(pg1), pg1, 6)
    expect_is(trail1, 'list')
    expect_equal(trail1[[1]], trail2[[1]])
    expect_equal(length(trail2), 2)
    expect_equal(length(trail1), 1)
    expect_equal(length(trail1[[1]]), 9)
    expect_named(trail1[[1]][[1]], c('down', 'up'))
    expect_equal(length(trail1[[1]][[1]][[1]]), 4)
    expect_equal(length(trail3[[1]][[1]][[1]]), 6)
})

test_that("anyParalogues works", {
    expect_equal(sum(anyParalogues(pg1)), 93)
})
