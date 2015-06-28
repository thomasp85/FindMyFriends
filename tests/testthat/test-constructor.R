context("Pangenome constructor")

location <- tempdir()
unzip(system.file('extdata', 'Mycoplasma.zip', package='FindMyFriends'),
      exdir=location)
genomeFiles <- list.files(location, full.names=TRUE, pattern='*.fasta')

test_that("Constructor gives the right classes", {
    expect_is(pangenome(genomeFiles, TRUE), 'pgFull')
    expect_is(pangenome(genomeFiles, TRUE, 'prodigal'), 'pgFullLoc')
    expect_is(pangenome(genomeFiles, TRUE, lowMem=TRUE), 'pgLM')
    expect_is(pangenome(genomeFiles, TRUE, 'prodigal', lowMem=TRUE), 'pgLMLoc')
})

test_that("Defaults are getting assigned", {
    expect_true(translated(pangenome(genomeFiles, TRUE)))
    expect_equal(defaults(pangenome(genomeFiles, TRUE, testDef='test'))$testDef, 'test')
})

test_that("Basic structure is sane", {
    expect_equal(length(pangenome(genomeFiles, TRUE)), length(genomeFiles))
    expect_equal(nGenes(pangenome(genomeFiles, TRUE)), 13649)
    expect_named(orgInfo(pangenome(genomeFiles, TRUE)), 'nGenes')
})