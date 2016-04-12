context("Pangenome modifications")

pg <- .loadPgExample(withGroups = TRUE)
pg1 <- .loadPgExample(withGroups = TRUE, withParalogues = TRUE)
groupInfo(pg1)$description <- as.character(1:nGeneGroups(pg))

test_that("removeGene works", {
    expect_equal(removeGene(pg, name=geneNames(pg)[100]), removeGene(pg, ind=100))
    expect_warning(removeGene(pg, name='not a gene name'))
    
    expect_equal(removeGene(pg, name=geneNames(pg)[20], organism=orgNames(pg)[1]), removeGene(pg, name=geneNames(pg)[20], organism=1))
    expect_equal(removeGene(pg, name=geneNames(pg)[20], organism=1), removeGene(pg, ind=20))
    expect_warning(removeGene(pg, name=geneNames(pg)[20], organism=2))
    
    expect_equal(removeGene(pg, organism=orgNames(pg)[3]), removeGene(pg, organism=3))
    expect_equal(removeGene(pg, organism=3), removeGene(pg, ind=which(seqToOrg(pg)==3)))
    expect_warning(removeGene(pg, organism=20))
    
    expect_equal(removeGene(pg, organism=orgNames(pg)[2:3], ind=20), removeGene(pg, organism=2:3, ind=20))
    expect_equal(removeGene(pg, organism=2:3, ind=20), removeGene(pg, ind=c(which(seqToOrg(pg)==2)[20], which(seqToOrg(pg)==3)[20])))
    expect_equal(removeGene(pg, organism=2:3, ind=c(20, 300)), removeGene(pg, ind=c(which(seqToOrg(pg)==2)[20], which(seqToOrg(pg)==3)[300])))
    expect_warning(removeGene(pg, organism=9, ind=30))
    expect_warning(removeGene(pg, organism=2, ind=50000))
    
    expect_equal(removeGene(pg, group='OG4'), removeGene(pg, group=4))
    expect_equal(removeGene(pg, group=4), removeGene(pg, ind=which(seqToGeneGroup(pg)==4)))
    expect_warning(removeGene(pg, group='not a group name'))
    expect_warning(removeGene(pg, group = -1))
    
    expect_equal(removeGene(pg, group='OG4', ind=8), removeGene(pg, group=4, ind=8))
    expect_equal(removeGene(pg, group=c('OG4', 'OG9'), ind=c(8, 2)), removeGene(pg, group=c(4, 9), ind=c(8, 2)))
    expect_equal(removeGene(pg, group=c(4, 9), ind=c(8, 2)), removeGene(pg, ind=c(which(seqToGeneGroup(pg)==4)[8], which(seqToGeneGroup(pg)==9)[2])))
    expect_warning(removeGene(pg, group=9, ind=sum(seqToGeneGroup(pg)==9)+1))
})

test_that("collapseParalogues works", {
    expect_error(collapseParalogues(pg), 'Paralogues must be detected before collapsing')
    expect_error(collapseParalogues(pg1, combineInfo='not a method'), 'Unknown combine function')
    pg2 <- collapseParalogues(pg1, sep='/')
    pg3 <- collapseParalogues(pg1)
    expect_is(groupInfo(pg2)$description, 'character')
    expect_is(groupInfo(pg3)$description, 'list')
    expect_equal(nGeneGroups(pg2), nGeneGroups(pg3))
    expect_equal(nGeneGroups(pg2), length(unique(groupInfo(pg1)$paralogue)))
    paraGroup <- which(seqToGeneGroup(pg1) %in% which(groupInfo(pg1)$paralogue==1))
    expect_true(all(paraGroup %in% which(seqToGeneGroup(pg2) == seqToGeneGroup(pg2)[paraGroup[1]])))
})

test_that("addXInfo works", {
    gInfo <- data.frame(test1=1:nGeneGroups(pg), test2=as.character(1:nGeneGroups(pg)), test3=I(as.list(1:nGeneGroups(pg))), stringsAsFactors=FALSE)
    oInfo <- data.frame(test1=1:nOrganisms(pg), test2=as.character(1:nOrganisms(pg)), test3=I(as.list(1:nOrganisms(pg))), stringsAsFactors=FALSE)
    pg2 <- addGroupInfo(pg, info=gInfo, key=1:nGeneGroups(pg))
    pg3 <- addGroupInfo(pg, info=gInfo, key='test1')
    pg4 <- addGroupInfo(pg, info=gInfo[1:5,], key=1:5)
    pg5 <- addOrgInfo(pg, info=oInfo, key=1:nOrganisms(pg))
    pg6 <- addOrgInfo(pg, info=oInfo, key='test1')
    pg7 <- addOrgInfo(pg, info=oInfo[2:3,], key=2:3)
    expect_equal(groupInfo(pg2)$test1, gInfo$test1)
    expect_equal(groupInfo(pg3)$test2, gInfo$test2)
    expect_equal(I(groupInfo(pg4)$test3)[1:5], gInfo$test3[1:5])
    expect_true(all(is.na(groupInfo(pg4)$test1[6:nGeneGroups(pg4)])))
    expect_equal(orgInfo(pg5)$test1, oInfo$test1)
    expect_equal(orgInfo(pg6)$test2, oInfo$test2)
    expect_equal(I(orgInfo(pg7)$test3)[2:3], oInfo$test3[2:3])
    expect_true(all(is.na(orgInfo(pg7)$test1[-(2:3)])))
})

test_that("Info merging works", {
    testData <- groupInfo(pg1)[1:4, ]
    mergedList <- mergeInfo(testData)
    mergedSep <- mergeInfo(testData, sep='/')
    mergedLargest <- largestInfo(testData)
    expect_is(mergedList, 'data.frame')
    expect_is(mergedSep, 'data.frame')
    expect_is(mergedLargest, 'data.frame')
    expect_named(mergedList, names(testData))
    expect_named(mergedSep, names(testData))
    expect_named(mergedLargest, names(testData))
    expect_is(mergedList$description, 'AsIs')
    expect_equal(length(mergedList$description[[1]]), nrow(testData))
    expect_is(mergedSep$description, 'character')
    expect_equal(length(strsplit(mergedSep$description, '/')[[1]]), 4)
    expect_equal(mergedLargest, testData[which.max(testData$nGenes),])
})

test_that("Pair-to-index works", {
    sampleData <- matrix(rnorm(20), ncol=4)
    locations <- data.frame(row=c(3, 2, 5, 2), col=c(1, 1, 4, 3))
    index <- pairToIndex(locations$row, locations$col, nrow(sampleData))
    expect_equal(sampleData[index[1]], sampleData[locations$row[1], locations$col[1]])
    expect_equal(sampleData[index[2]], sampleData[locations$row[2], locations$col[2]])
    expect_equal(sampleData[index[3]], sampleData[locations$row[3], locations$col[3]])
    expect_equal(sampleData[index[4]], sampleData[locations$row[4], locations$col[4]])
})

test_that("Removing index works", {
    expect_equal(removeIndex(c(1,1,2,5,4,5,2)), c(1,1,2,4,3,4,2))
})
