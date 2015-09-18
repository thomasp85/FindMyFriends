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