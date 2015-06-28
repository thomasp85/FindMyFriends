context("Using object defaults")

pg <- .loadPgExample()
defaultTest <- function(x, defaults) {
    .fillDefaults(defaults)
    x
}

test_that("Setter and getter works", {
    expect_equal(defaults(pg)$kmerSize, defaultArgs$kmerSize)
    defaults(pg)$kmerSize <- 1
    expect_equal(defaults(pg)$kmerSize, 1)
    defaults(pg)$newDef <- 'test'
    expect_equal(defaults(pg)$newDef, 'test')
    expect_error(defaults(pg)$translated <- FALSE)
})

test_that("defaults can be used with .fillDefaults", {
    expect_error(defaultTest(defaults=list(y=1)))
    expect_equal(defaultTest(defaults=list(x=1)), 1)
    expect_equal(defaultTest(x=2,defaults=list(x=1)), 2)
})