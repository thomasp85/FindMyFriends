context("Example pangenome loading")

pgFull <- .loadPgExample()
pgLM <- .loadPgExample(lowMem=TRUE)
pgFullLoc <- .loadPgExample(geneLoc=TRUE)
pgLMLoc <- .loadPgExample(lowMem=TRUE, geneLoc=TRUE)
pgFullGr <- .loadPgExample(withGroups=TRUE)
pgLMGr <- .loadPgExample(lowMem=TRUE, withGroups=TRUE)
pgFullLocGr <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE)
pgLMLocGr <- .loadPgExample(lowMem=TRUE, geneLoc=TRUE, withNeighborhoodSplit=TRUE)
pgFullPara <- .loadPgExample(withGroups=TRUE, withParalogues=TRUE)
pgFullLocPara <- .loadPgExample(geneLoc=TRUE, withNeighborhoodSplit=TRUE, withParalogues=TRUE)
pgLMPara <- .loadPgExample(lowMem=TRUE, withGroups=TRUE, withParalogues=TRUE)
pgLMLocPara <- .loadPgExample(lowMem=TRUE, geneLoc=TRUE, withNeighborhoodSplit=TRUE, withParalogues=TRUE)

test_that("test object are correct class", {
    expect_is(pgFull, 'pgFull')
})
