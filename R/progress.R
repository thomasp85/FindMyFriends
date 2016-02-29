makeProgress <- function(end, name, freq, namelength = nchar(name), maxwait = 10) {
    nameMod <- unlist(strsplit(name, split = ''))
    nameMod <- nameMod[seq_len(namelength)]
    nameMod[is.na(nameMod)] <- ' '
    nameMod <- paste(nameMod, collapse = '')
    if (end/freq < 50) {
        freq <- floor(end/50)
    }
    prog <- list(
        end = end,
        name = nameMod,
        freq = freq,
        prog = 0,
        last = 0,
        maxwait = maxwait,
        lasttime = as.numeric(Sys.time())
    )
    createBar(prog)
    list2env(prog, parent = emptyenv())
}

progress <- function(progress) {
    progress$prog <- progress$prog + 1
    progress$last <- progress$last + 1
    currentTime <- as.numeric(Sys.time())
    if (progress$end == progress$prog) {
        createBar(progress)
    } else if (progress$last > progress$freq || currentTime - progress$lasttime > progress$maxwait) {
        progress$lasttime <- currentTime
        progress$last <- 0
        createBar(progress)
    }
    flush.console()
}
createBar <- function(progress) {
    p <- round(progress$prog / progress$end * 50)
    whites <- paste(rep(' ', nchar(progress$end) - nchar(progress$prog)), collapse = '')
    pString <- paste(c(rep('=', p), rep(' ', 50 - p)), collapse = '')
    pString <- paste0(progress$name, ' |', pString, '| ')
    pString <- paste0('\r', pString, whites, progress$prog, '/', progress$end)
    cat(pString)
    flush.console()
}