#' @include aaa.R
NULL

#' Read GFF file into GRanges object
#' 
#' 
#' 
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom GenomeInfoDb seqlevels isCircular
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' 
#' @export
#' 
readGFF <- function(file) {
    regions <- scan(file, character(), sep = '\n', quiet = TRUE)
    fastaStart <- which(regions == '##FASTA')
    if (length(fastaStart) == 0) {
        stop("Raw sequence must be included in GFF file")
    }
    if (!grepl('##gff-version\\s+3', regions[1])) {
        stop("Only GFF3 supported")
    }
    regions <- regions[1:(fastaStart - 1)]
    regionInfo <- regions[grepl('##sequence-region', regions)]
    regions <- regions[!grepl('^#', regions)]
    regions <- read.table(text = regions, sep = '\t', stringsAsFactors = FALSE)
    attributes <- parseAttributes(regions[,9])
    attributes <- attributes[, !names(attributes) %in% c("seqnames", "ranges", 
                                                        "strand", "seqlevels", 
                                                        "seqlengths", 
                                                        "isCircular", "start", 
                                                        "end", "width", 
                                                        "element", "type", 
                                                        "score", "phase")]
    gff <- GRanges(seqnames = Rle(regions[,1]), 
            ranges = IRanges(start = regions[,4], 
                             end = regions[,5]), 
            strand = regions[,7], type = regions[,3], 
            score = regions[, 6], phase = regions[,8], attributes)
    if (length(regionInfo) != 0) {
        regionInfo <- strsplit(sub('##sequence-region\\s+', '', regionInfo), 
                               '\\s+')
        regLengths <- sapply(regionInfo, function(inf) {
            regL <- as.integer(inf[3]) - as.integer(inf[2]) + 1
            names(regL) <- inf[1]
            regL
        })
        seqlengths(gff) <- regLengths[match(seqlevels(gff), names(regLengths))]
    }
    if (!is.null(gff$Is_circular)) {
        circ <- sapply(split(gff$Is_circular, seqnames(gff)), `[`, i = 1)
        isCircular(gff) <- circ[match(seqlevels(gff), names(circ))] == 'true'
    }
    gff
}
#' Read gff file into DNAStringSet
#' 
#' This function parses a gff file into a DNAStringSet provided that sequence
#' information is included in the bottom of the file (not a requirement for gff
#' files). The features can be filtered upfront by type (such as CDS and exon). 
#' GFF annotation information is stored alongside the DNAStringSet in the 
#' elementMetadata slot as a DataFrame (unfortunately it is not possible to 
#' assign GRanges object to elementMetadata).
#' 
#' @param file The path to a GFF file.
#' 
#' @param type Character vector or NULL. The feature types to include in the 
#' output. NULL means include all.
#' 
#' @return A DNAStringSet with the sequences of the features annotated in the 
#' GFF file.
#' 
#' @importFrom Biostrings readDNAStringSet extractAt reverseComplement
#' @importFrom GenomicRanges seqnames split strand
#' 
#' @export
#' 
gff2StringSet <- function(file, type = 'CDS') {
    gff <- readGFF(file)
    if (!is.null(type)) {
        gff <- gff[gff$type %in% type]
    }
    sequences <- try(readDNAStringSet(file, seek.first.rec = TRUE), TRUE)
    if (inherits(sequences, 'try-error')) {
        stop(file, 'does not include sequence information')
    }
    seqNames <- sub('^(\\S+)\\s.*$', '\\1', names(sequences))
    gff <- split(gff, seqnames(gff))
    gff <- gff[match(seqNames, names(gff))]
    names(gff) <- NULL
    sequences <- unlist(extractAt(sequences, as(gff, 'RangesList')))
    gff <- unlist(gff)
    oppStrand <- strand(gff) == '-'
    sequences[oppStrand] <- reverseComplement(sequences[oppStrand])
    if (is.null(gff$ID)) {
        gff$ID <- 1:length(gff)
    } else {
        gff$ID <- make.unique(gff$ID)
    }
    names(sequences) <- paste0(seqnames(gff), '_', gff$ID)
    sequences@elementMetadata <- as(gff, 'DataFrame')
    sequences
}
#' Import annotation from an .annot file
#' 
#' This function makes it easy to import annotation create in Blast2GO or other
#' programs supporting .annot exporting of results.
#' 
#' @param file The .annot file to import
#' 
#' @return A data.frame ready to merge with a pangenome object using 
#' \code{\link{addGroupInfo}} with the \code{key} argument set to 'name'.
#' 
#' @examples 
#' # Get path to file
#' annot <- system.file('extdata', 'examplePG', 'example.annot', package='FindMyFriends')
#' 
#' # Parse the file
#' readAnnot(annot)
#' 
#' @export
#' 
readAnnot <- function(file) {
    data <- read.table(file, header = FALSE, sep = '\t', fill = TRUE, 
                       stringsAsFactors = FALSE)
    names(data) <- c('name', 'annot', 'desc')
    data <- data %>% 
        mutate(ontology = grepl('GO:', annot)) %>%
        group_by(name) %>%
        summarise(description = desc[1], GO = I(list(annot[ontology])), 
                  EC = I(list(annot[!ontology])))
    data <- as.data.frame(data)
    data$GO <- unclass(data$GO)
    data$EC <- unclass(data$EC)
    data
}
#' Extracts location information from prodigal generated files
#' 
#' This function parses prodigal generated sequence names and extract the 
#' contig, start, stop and strand of each protein.
#' 
#' @param desc A character vector of sequence descriptions
#' 
#' @return A data.frame with the columns contig, start, end and strand
#' 
#' @noRd
#' 
prodigalParse <- function(desc) {
    info <- do.call(rbind, strsplit(desc, ' # '))[,-5]
    info[, 1] <- sub('_\\d+$', '', info[,1])
    info <- data.frame(
        contig = info[, 1], 
        start = as.integer(info[,2]), 
        end = as.integer(info[,3]), 
        strand = as.integer(info[,4]), 
        stringsAsFactors = FALSE)
    info
}

## HELPERS

#' Parse key-val into data.frame
#' 
#' This function takes a set of strings of the format: key=val;key=val, where 
#' '=' is set by valSep and ';' is set by attrSep, and parses it into a 
#' data.frame with a column for each unique key and a row for each element in 
#' the character vector. Missing key in some elements will be filled with NA.
#' 
#' @param strings A character vector of strings to parse
#' 
#' @param valSep The token separating the key and value
#' 
#' @param attrSep The token separating key-value pairs
#' 
#' @return A data.frame with a row per element in strings, and a column for each
#' unique key
#' 
#' @noRd
#' 
parseAttributes <- function(strings, valSep = '=', attrSep = ';') {
    attr <- strsplit(strings, split = attrSep)
    ind <- rep(1:length(strings), lengths(attr))
    attr <- unlist(attr)
    attrKey <- sub(paste0('^(.*?)', valSep, '.*$'), '\\1', attr)
    attrVal <- sub(paste0('^.*?', valSep, '(.*)$'), '\\1', attr)
    attrNames <- unique(attrKey)
    attr <- lapply(attrNames, function(name) {
        val <- rep(NA, length(strings))
        elements <- which(attrKey == name)
        val[ind[elements]] <- attrVal[elements]
        type.convert(val, as.is = TRUE)
    })
    attr <- as.data.frame(attr, stringsAsFactors = FALSE)
    names(attr) <- attrNames
    attr
}
#' @importFrom GenomicRanges seqnames start end strand
#' 
createDescriptionString <- function(gRange, sep=';') {
    paste0(seqnames(gRange), sep, 
           start(gRange), sep, 
           end(gRange), sep,
           strand(gRange))
}