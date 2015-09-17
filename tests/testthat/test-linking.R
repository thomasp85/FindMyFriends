context("Gene group linking")

pg <- .loadPgExample(withGroups = T, withNeighborhoodSplit = T)

test_that("kmerLink works", {
    set.seed(1)
    pg2 <- kmerLink(pg, lowerLimit=0.8)
    set.seed(1)
    pg3 <- kmerLink(pg, lowerLimit=0.8, pParam=BiocParallel::SerialParam())
    expect_equal(groupInfo(pg2)$paralogue, groupInfo(pg3)$paralogue)
    expect_equal(length(unique(groupInfo(pg2)$paralogue)), 3038)
    if (R.version$arch == 'i386') {
        expect_equal(groupInfo(pg2)$paralogue[1:6], c(2, 2, 1, 111, 55, 112))
    } else {
        expect_equal(groupInfo(pg2)$paralogue[1:6], c(2, 2, 1, 116, 55, 117))
    }
})