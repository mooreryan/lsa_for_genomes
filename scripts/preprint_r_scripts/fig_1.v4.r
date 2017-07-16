library("lsa")
library("Matrix")
library("ape")

source("functions.r")

hw <- 8
pdf("fig_1.2.pdf", height=hw, width=hw)
par(mfrow=c(2,2))
choices <- 1L:2L
arrow.weight <- 0.8
col.1 <- "#E394F1" # purple
col.1.dark <- "#921E7D"
col.2 <- "#F3C25A"
col.2.dark <- "#F38400" # orange
col.3 <- "#A1CAF1" # light.blue
col.3.dark <- "#0075F1"

enterobacteria <- c("E_coli", "S_flexneri")
methanogens <- c("M_maripaludis", "M_mazei")
cyanobacteria <- c("P_marinus", "S_elongatus")
color.of <- function(label) {
    if(label == enterobacteria[1]) {
        col.1
    } else if(label == enterobacteria[2]) {
        col.1.dark
    } else if(label == methanogens[1]) {
        col.2
    } else if(label == methanogens[2]) {
        col.2.dark
    } else if(label == cyanobacteria[1]) {
        col.3
    } else if(label == cyanobacteria[2]) {
        col.3.dark
    } else {
        "black"
    }
}

label.col <- function(node) {
    if (is.leaf(node)) {
        ## fetch label
        label <- attr(node, "label")
        attrs <- attributes(node)
        print("AAPLE")
        print(attrs)

        ## set label color
        attr(node, "nodePar") <- c(attrs$nodePar, list(lab.col = color.of(label)))
        print("pie")
        print(attributes(node))
    }

    node
}

dir <- "/Users/moorer/Documents/Lab Notebook/LatentSemanticAnalysis/2017 july preprint/reference genomes/mga_orfs/"
## dat <- readMM("/Users/moorer/Documents/Lab Notebook/LatentSemanticAnalysis/2017 july preprint/reference genomes/mga_orfs/all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt")

sing.vals <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.singular_values.txt", sep=""))
sing.vals <- as.matrix(sing.vals)

transformed.doc.matrix <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt", sep=""))
transformed.doc.matrix <- as.matrix(transformed.doc.matrix)

idx.to.doc <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.idx_to_doc_map.txt", sep=""),
                         col.names="name")

transformed.doc.matrix.dist <- read.table(paste(dir,
                                                "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt",
                                                sep=""),
                                          col.names=idx.to.doc$name)
transformed.doc.matrix.dist <- as.dist(transformed.doc.matrix.dist)
transformed.doc.matrix.cluster <- hclust(transformed.doc.matrix.dist,
                                         method="average")

rows.are.terms <- as.matrix(read.table(paste(dir,
                                             "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.rows_are_terms.txt",
                                             sep="/")))

p <- as.phylo(transformed.doc.matrix.cluster)
colors <- unlist(lapply(p$tip.label, color.of))
plot(p,
     label.offset = 0.01,
     tip.color=colors,
     direction = "rightwards",
     align.tip.label = FALSE,
     main = "LSA hierarchical clustering",
     frame.plot=F)

var.explained <- sing.vals ^ 2 / sum(sing.vals ^ 2) * 100
sum.of.var <- c();for(i in 1:length(var.explained)) { sum.of.var <- c(sum.of.var, sum(var.explained[1:i]))}
plot(sum.of.var, type="o", xlab="Topic", ylab="% variance expained", main="Cumulative variance explained")

x <- transformed.doc.matrix[, 1:2]
y <- rows.are.terms[, 1:2]
y.mag <- apply(y, MARGIN=2, FUN=function(x){ abs(x) })
top.n <- 5
top.5.terms.topic.1 <- y.mag[order(y.mag[,1], decreasing=T),][1:top.n,]
top.5.terms.topic.2 <- y.mag[order(y.mag[,2], decreasing=T),][1:top.n,]
top.5.terms <- rbind(top.5.terms.topic.1, top.5.terms.topic.2)
biplot.default.2(x,
                 top.5.terms,
                 main="Gene weights, genomes in topic space",
                 xlab="Topic 1",
                 ylab="Topic 2",
                 ylim=c(-5,150), # this is for the topic vectors in term space
                 x.plot.col=colors,
                 arrow.col=rgb(0.5,0,0,0.5),
                 cex=1.25
                 )

x <- transformed.doc.matrix[, 2:3]
y <- rows.are.terms[, 2:3]
y.mag <- apply(y, MARGIN=2, FUN=function(x){ abs(x) })
top.n <- 5
top.5.terms.topic.2 <- y.mag[order(y.mag[,1], decreasing=T),][1:top.n,]
top.5.terms.topic.3 <- y.mag[order(y.mag[,2], decreasing=T),][1:top.n,]
top.5.terms <- rbind(top.5.terms.topic.2, top.5.terms.topic.3)
## swap the columns of the matrices so that topic 2 is on the y axis both times
top.5.terms.rev <- cbind(top.5.terms[,2], top.5.terms[,1])
x.rev <- cbind(x[,2], x[,1])
biplot.default.2(x.rev,
                 top.5.terms.rev,
                 main="Gene weights, genomes in topic space",
                 xlab="Topic 3",
                 ylab="Topic 2",
                 xlim=c(-30,110),
                 ylim=c(-5,150),
                 x.plot.col=colors,
                 arrow.col=rgb(0.5,0,0,0.5),
                 cex=1.25,
                 )
dev.off()
