context("pgLM class methods")

pg <- .loadPgExample(lowMem=T, withGroups = T, withParalogues = T)
location <- tempdir()
unzip(system.file('extdata', 'Mycoplasma.zip', package='FindMyFriends'),
      exdir=location)
genomeFiles <- list.files(location, full.names=TRUE, pattern='*.fasta')[1:5]
realGenes <- Biostrings::readAAStringSet(genomeFiles)

test_that("genes getter works", {
    expect_equal(genes(pg), realGenes)
    subset <- sample(nGenes(pg), 10)
    expect_equal(as.character(genes(pg, subset=subset)), as.character(realGenes[subset]))
    expect_is(genes(pg, split = 'organism'), 'AAStringSetList')
    expect_equal(length(genes(pg, split='organism')), nOrganisms(pg))
    expect_equal(as.character(genes(pg, split='organism')[[3]]), as.character(genes(pg, split='organism', subset=3)[[1]]))
    expect_equal(as.character(genes(pg, split='organism')[[3]]), as.character(realGenes[seqToOrg(pg)==3]))
    expect_is(genes(pg, split = 'group'), 'AAStringSetList')
    expect_equal(length(genes(pg, split='group')), nGeneGroups(pg))
    expect_equal(as.character(genes(pg, split='group')[[3]]), as.character(genes(pg, split='group', subset=3)[[1]]))
    expect_equal(as.character(genes(pg, split='group')[[3]]), as.character(realGenes[seqToGeneGroup(pg)==3]))
    expect_is(genes(pg, split = 'paralogue'), 'AAStringSetList')
    expect_equal(length(genes(pg, split='paralogue')), length(unique(groupInfo(pg)$paralogue)))
    expect_equal(as.character(genes(pg, split='paralogue')[[3]]), as.character(genes(pg, split='paralogue', subset=3)[[1]]))
    expect_equal(as.character(genes(pg, split='paralogue')[[3]]), as.character(realGenes[seqToGeneGroup(pg) %in% which(groupInfo(pg)$paralogue == 3)]))
})

test_that("Gene names setter and getter works", {
    expect_equal(geneNames(pg), names(realGenes))
    geneNames(pg) <- as.character(1:nGenes(pg))
    expect_equal(geneNames(pg), as.character(1:nGenes(pg)))
    geneNames(pg) <- names(realGenes)
})

test_that("Gene width works", {
    expect_equal(geneWidth(pg), Biostrings::width(realGenes))
})

