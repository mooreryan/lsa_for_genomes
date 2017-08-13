Sys.setlocale(locale = "en_US.UTF-8")
Sys.setenv(LANG = "en")

biplot2 <- function(doc.scores,
                    term.loadings,
                    doc.names,
                    num.terms.to.keep = 0,
                    topic.x = 1,
                    topic.y = 2,
                    doc.color=rgb(0, 0, 0, 0.7),
                    term.color=rgb(1, 0, 0, 0.35))
{
    ## Lims for square plot that includes the origin
    all.points <- unlist(c(0,
                           doc.scores[, c(topic.x, topic.y)],
                           term.loadings[, c(topic.x, topic.y)]))
    lims <- c(min(all.points),
              max(all.points)) * 1.1

    par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)

    plot(doc.scores[, c(topic.x, topic.y)],
         xlab="",
         ylab="",
         bty="n",
         axes=F,
         type="n",
         xlim=lims,
         ylim=lims)

    grid(lwd=1)
    box()
    axis(1)
    axis(2, las=1)
    title(xlab = paste("Topic", topic.x, sep=" "),
          ylab = paste("Topic", topic.y, sep=" "),
          line=2.5)

    if (num.terms.to.keep == 0) {
        n <- nrow(term.loadings)
    } else {
        n <- round(num.terms.to.keep / 2)
    }

    abs.term.loadings <- apply(term.loadings, MARGIN=2, FUN=abs)
    top.term.loadings.x <- order(abs.term.loadings[, topic.x], decreasing=T)[1:n]
    top.term.loadings.y <- order(abs.term.loadings[, topic.y], decreasing=T)[1:n]

    top.term.loadings <- rbind(term.loadings[top.term.loadings.x, ],
                               term.loadings[top.term.loadings.y, ])

    arrows(0,
           0,
           x1=top.term.loadings[, topic.x],
           y1=top.term.loadings[, topic.y],
           col=term.color,
           length=0.08)

    offset <- abs(lims[1] - lims[2]) / 35

    actual.scores <- cbind(doc.scores[, topic.x],
                           doc.scores[, topic.y])

    offset.scores <- cbind(doc.scores[, topic.x],
                           doc.scores[, topic.y] - offset)

    points(actual.scores,
           pch=16,
           cex=0.8,
           col=doc.color)
    text(offset.scores,
         labels=doc.names,
         cex=0.8,
         col=doc.color)
}

US <- read.table("/Users/moorer/projects/lsa_for_genomes/output/metadata_groups/original/redsvd/svd.US", header=F, sep=" ")

VS <- read.table("/Users/moorer/projects/lsa_for_genomes/output/metadata_groups/original/redsvd/svd.VS", header=F, sep=" ")

n <- nrow(VS)

## Projections scaled to unit variance
doc.projections <- VS / sqrt(n - 1)

term.loadings <- US / sqrt(n - 1)

pdf("test.pdf", width=8, height=8)
biplot2(doc.projections,
        term.loadings,
        c("E. coli", "S. flexneri", "M. Mazei"),
        num.terms.to.keep=50,
        topic.x=1,
        topic.y=2)
invisible(dev.off())
