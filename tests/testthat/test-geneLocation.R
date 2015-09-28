context("Parsing of gene location")

pg <- .loadPgExample()
genenames <- head(geneNames(pg))
badFunction <- function(desc) {
    data.frame(contig=desc, stringsAsFactors=FALSE)
}

test_that("Prodigal parser works", {
    parsedNames <- prodigalParse(genenames)
    expect_named(parsedNames, c('contig', 'start', 'end', 'strand'))
    expect_equal(nrow(parsedNames), length(genenames))
    expect_is(parsedNames$contig, 'character')
    expect_is(parsedNames$start, 'integer')
    expect_is(parsedNames$end, 'integer')
    expect_is(parsedNames$strand, 'integer')
    expect_false(any(sapply(parsedNames, is.na)))
})

test_that("Gene location wrapper works", {
    parsedNames <- prodigalParse(genenames)
    expect_equal(getSeqInfo(parsedNames, genenames), parsedNames)
    expect_equal(getSeqInfo('prodigal', genenames), parsedNames)
    expect_equal(getSeqInfo(prodigalParse, genenames), parsedNames)
    expect_equal(getSeqInfo(split(parsedNames, 1:6), genenames), parsedNames)
    expect_error(getSeqInfo(badFunction, genenames))
    expect_error(getSeqInfo('nonsense', genenames))
    expect_error(getSeqInfo(parsedNames, genenames[1:4]))
})