## this is a modification of the biplot.defualt function
biplot.default.2 <- function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
                  xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
                  arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                  x.plot.col = NULL, arrow.col = NULL,
                  ...)
{
    n <- nrow(x)
    p <- nrow(y)
    if (missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if (is.null(xlabs))
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if (missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if (is.null(ylabs))
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
    if (length(cex) == 1L)
        cex <- c(cex, cex)
    if (missing(col)) {
        col <- par("col")
        if (!is.numeric(col))
            col <- match(col, palette(), nomatch = 1L)
        col <- c(col, col + 1L)
    }
    else if (length(col) == 1L)
        col <- c(col, col)
    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),
                                    abs(max(x, na.rm = TRUE)))

    x.ranges <- apply(x, 2, range)
    y.ranges <- apply(y, 2, range)
    ratio <- max(y.ranges[,1] / x.ranges[,1],
                 y.ranges[,2] / x.ranges[,2]) / expand

    ## Since the weight vectors start from (0,0), make sure zero is
    ## within both x and y of x.ranges
    if (x.ranges[1,1] > 0) { # both xlims are positive
        x.ranges[1,1] <- 0 # get zero in the range
    } else if (x.ranges[2,1] < 0) { # both xlims are negative
        x.ranges[2,1] <- 0 # get zero in the range
    }

    if (x.ranges[1,2] > 0) { # both y lims are positive
        x.ranges[1,2] <- 0 # get zero in the range
    } else if (x.ranges[2,2] < 0) { # both y lims are negative
        x.ranges[2,2] <- 0 # get zero in the range
    }

    if (missing(xlim) && missing(ylim)) {
        xlim <- x.ranges[,1]
        ylim <- x.ranges[,2]
    } else if (missing(xlim)) {
        xlim <- x.ranges[,1]
    } else if (missing(ylim)) {
        ylim <- x.ranges[,2]
    }
    print("xlim")
    print(xlim)
    print("ylim")
    print(ylim)
    ## ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    print("ratio")
    print(ratio)
    on.exit(par(op))
    op <- par(pty = "s")
    if (!is.null(main))
        op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
    plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim *
             ratio, xlab = "", ylab = "", col = col[1L], ...)
    axis(3, col = col[2L], ...)
    axis(4, col = col[2L], ...)
    box(col = col[1L])
    if (var.axes) {
        if (missing(arrow.col)) {
            ## text(y, labels = ylabs, cex = cex[2L], col = col[2L], ...)
            arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = col[2L],
                   length = arrow.len)
        } else {
            ## text(y, labels = ylabs, cex = cex[2L], col = arrow.col, ...)
            arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = arrow.col,
                   length = arrow.len)
        }
    }
    par(new = TRUE)
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    if (missing(x.plot.col)) {
        points(x, cex = cex[1L], col = col[1L], pch=16, ...)
    } else {
        points(x, cex = cex[1L], col=x.plot.col, pch=16, ...)
    }
    invisible()
}
