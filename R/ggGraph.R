# ggGraph <- function(graph, layout) {
#     if(!inherits(graph, 'igraph')) {
#         stop('graph must be an igraph object')
#     }
#     if(inherits(layout, 'function')) {
#         layout <- layout(graph)
#     } else if (!inherits(layout, 'matrix')) {
#         stop('layout must by a matrix of positions or a function returning one')
#     }
#     if(ncol(layout) < 2) {
#         stop('layout must have at least 2 columns')
#     }
#     if(nrow(layout) != length(V(graph))) {
#         stop('layout must have the same number of rows as the number of vertices')
#     }
#     vFrame <- data.frame(x=layout[,1], y=layout[,2], label=V(graph)$name, type=V(graph)$centerGroup)
#     edges <- get.edgelist(graph, FALSE)
#     eFrame <- data.frame(x=vFrame$x[edges[,1]], xend=vFrame$x[edges[,2]], y=vFrame$y[edges[,1]], yend=vFrame$y[edges[,2]], weight=E(graph)$weight)
#     p <- ggplot() + theme_minimal()
#     p <- p + geom_segment(aes(x=x, xend=xend, y=y, yend=yend, size=weight), arrow=arrow(length=unit(1, 'cm')), data=eFrame, colour='grey')
#     p <- p + geom_point(aes(x=x, y=y, colour=type), data=vFrame, size=12) + scale_colour_manual('', labels=c('Center', 'Neighborhood'), breaks=c(TRUE, FALSE), values=c('steelblue', 'forestgreen'))
#     p <- p + geom_text(aes(x=x, y=y, label=label), data=vFrame, size=3, family='sans')
#     p
# }
# 
# geom_segment2 <- function (mapping = NULL, data = NULL, stat = "identity",
#                           position = "identity", arrow = NULL, lineend = "butt", na.rm = FALSE, startAdjust = NULL, endAdjust = NULL, ...) {
#     
#     GeomSegment2$new(mapping = mapping, data = data, stat = stat,
#                     position = position, arrow = arrow, lineend = lineend, na.rm = na.rm, startAdjust = startAdjust, endAdjust = endAdjust, ...)
# }
# 
# GeomSegment2 <- proto(ggplot2:::Geom, {
#     objname <- "segment2"
#     
#     draw <- function(., data, scales, coordinates, arrow = NULL,
#                      lineend = "butt", na.rm = FALSE, startAdjust = startAdjust, endAdjust = endAdjust, ...) {
#         
#         data <- remove_missing(data, na.rm = na.rm,
#                                c("x", "y", "xend", "yend", "linetype", "size", "shape"),
#                                name = "geom_segment2")
#         if (empty(data)) return(zeroGrob())
#         
#         if (is.linear(coordinates)) {
#             return(with(coord_transform(coordinates, data, scales), {
#                 segmentsGrob2(x, y, xend, yend, default.units="native", startAdjust=startAdjust, endAdjust=endAdjust,
#                              gp = gpar(col=alpha(colour, alpha), fill = alpha(colour, alpha),
#                                        lwd=size * .pt, lty=linetype, lineend = lineend),
#                              arrow = arrow)
#             }
#             ))
#         }
#         
#         data$group <- 1:nrow(data)
#         starts <- subset(data, select = c(-xend, -yend))
#         ends <- rename(subset(data, select = c(-x, -y)), c("xend" = "x", "yend" = "y"),
#                        warn_missing = FALSE)
#         
#         pieces <- rbind(starts, ends)
#         pieces <- pieces[order(pieces$group),]
#         
#         GeomPath$draw_groups(pieces, scales, coordinates, arrow = arrow, ...)
#     }
#     
#     
#     default_stat <- function(.) StatIdentity
#     required_aes <- c("x", "y", "xend", "yend")
#     default_aes <- function(.) aes(colour="black", size=0.5, linetype=1, alpha = NA)
#     guide_geom <- function(.) "path"
# })
# segmentsGrob2 <- function(x0 = unit(0, "npc"), y0 = unit(0, "npc"), x1 = unit(1, "npc"), y1 = unit(1, "npc"), startAdjust = unit(0, 'npc'), endAdjust = unit(0, 'npc'), default.units = "npc", arrow = NULL, name = NULL, gp = gpar(), vp = NULL) {
#     if (!is.unit(x0)) 
#         x0 <- unit(x0, default.units)
#     if (!is.unit(x1)) 
#         x1 <- unit(x1, default.units)
#     if (!is.unit(y0)) 
#         y0 <- unit(y0, default.units)
#     if (!is.unit(y1)) 
#         y1 <- unit(y1, default.units)
#     grid.draw(grob(x0 = x0, y0 = y0, x1 = x1, y1 = y1, startAdjust=startAdjust, endAdjust=endAdjust, arrow = arrow, name=name, gp=gp, vp=vp, cl="segments2"))
# }
# drawDetails.segments2 <- function(x, recording=TRUE) {
#     devSize <- dev.size()
#     transformation <- matrix(c(devSize[1], 0, 0, devSize[2]), ncol=2)
#     newVec <- cbind(as.numeric(x$x1)-as.numeric(x$x0), as.numeric(x$y1)-as.numeric(x$y0)) %*% transformation
#     segAngle <- atan2(newVec[,2], newVec[, 1])
#     xAdjust <- cos(segAngle)
#     yAdjust <- sin(segAngle)
#     if(!is.null(x$startAdjust)) {
#         x$x0 <- x$x0 + unit(as.numeric(x$startAdjust)*xAdjust, attributes(x$startAdjust)$unit)
#         x$y0 <- x$y0 + unit(as.numeric(x$startAdjust)*yAdjust, attributes(x$startAdjust)$unit)
#     }
#     if(!is.null(x$endAdjust)) {
#         x$x1 <- x$x1 - unit(as.numeric(x$endAdjust)*xAdjust, attributes(x$endAdjust)$unit)
#         x$y1 <- x$y1 - unit(as.numeric(x$endAdjust)*yAdjust, attributes(x$endAdjust)$unit)
#     }
#     if(!is.null(x$arrow)) {
#         x$x0 <- x$x0 + unit(ifelse(x$arrow$ends %in% c(1,3), xAdjust*(x$size * ggplot2:::.pt)/96, 0), 'inch')
#         x$y0 <- x$y0 + unit(ifelse(x$arrow$ends %in% c(1,3), yAdjust*(x$size * ggplot2:::.pt)/96, 0), 'inch')
#         x$x1 <- x$x1 - unit(ifelse(x$arrow$ends %in% c(2,3), xAdjust*(x$size * ggplot2:::.pt)/96, 0), 'inch')
#         x$y1 <- x$y1 - unit(ifelse(x$arrow$ends %in% c(2,3), yAdjust*(x$size * ggplot2:::.pt)/96, 0), 'inch')
#     }
#     grid.segments(x0 = x$x0, y0 = x$y0, x1 = x$x1, y1 = x$y1, arrow = x$arrow, name=x$name, gp=x$gp, vp=x$vp)
# }
# geom_circle <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", ...) {
#     GeomCircle$new(mapping = mapping, data = data, stat = stat, position = position, ...)
# }
# 
# GeomCircle <- proto(ggplot2:::Geom, {
#     objname <- "circle"
#     
#     default_stat <- function(.) StatIdentity
#     default_pos <- function(.) PositionIdentity
#     default_aes <- function(.) aes(colour=NA, fill="grey50", size=0.5, radius=5, linetype=1, alpha = NA)
#     
#     required_aes <- c("x", "y")
#     
#     draw <- function(., data, scales, coordinates, ...) {
#         return(with(coord_transform(coordinates, data, scales), {
#             circleGrob(x, y, unit(radius, 'mm'), default.units="native",
#                          gp = gpar(col=alpha(colour, alpha), fill = alpha(fill, alpha),
#                                    lwd=size * .pt, lty=linetype))
#         }
#         ))
#     }
#     guide_geom <- function(.) "polygon"
#     
# })