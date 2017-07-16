## install.packages("lsa")
## install.packages("Matrix")
## install.packages("ape")


library("lsa")
library("Matrix")
library("ape")

source("functions.r")

dir <- "/work/moorer/orfs_minus_bad_iron_mat"
choices <- 1L:2L
arrow.weight <- 0.8
top.n <- 5
## col.1 <- "#E394F1" # purple
## col.1.dark <- "#921E7D"
## col.2 <- "#F3C25A"
## col.2.dark <- "#F38400" # orange
## col.3 <- "#A1CAF1" # light.blue
## col.3.dark <- "#0075F1"

estuary <- c("estuarine_columbia_river_estuary_faa")
lake <- c("freshwater_lake_trout_bog_lake_1_faa",
           "saline_lake_ace_lake_antarctica_faa")
soil <- c("soil_arctic_peat_barrow_env_obs_alaska_faa",
          "soil_grasslands_andelo_coastal_reserve_california_faa",
          "soil_permafrost_arctic_organic_soil_faa",
          "soil_permafrost_arctic_permafrost_soil_faa")
surface.iron.mat <- c("iron_mat_479_BS4_faa",
                      "iron_mat_481_BS4_faa",
                      "iron_mat_483_BS63_faa")
vent.orifice.iron.mat <- c("iron_mat_482_BS8_faa",
                           "iron_mat_476_BS1_faa")
iron.rich.groundwater <- c("groundwater_crystal_geyser_aquifers_utah_faa")
other.groundwater <- c("groundwater_sediment_aquifer_colorado_faa")
deadish <- c("pelagic_marine_north_sea_faa")


color.of <- function(label) {
    if(label %in% estuary) {
        "mediumaquamarine"
    } else if(label %in% lake) {
        "blue"
    } else if(label %in% soil) {
        "green"
    } else if(label %in% surface.iron.mat) {
        "orange"
    } else if(label %in% vent.orifice.iron.mat) {
        "darkorange3"
    } else if(label %in% iron.rich.groundwater) {
        "red"
    } else if(label %in% deadish) {
        "brown"
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

## dat <- readMM(paste(dir,

sing.vals <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.singular_values.txt", sep="/"))
sing.vals <- as.matrix(sing.vals)

transformed.doc.matrix <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt", sep="/"))
transformed.doc.matrix <- as.matrix(transformed.doc.matrix)

idx.to.doc <- read.table(paste(dir, "all_prepped.fa.clu.tsv.sorted.seanie_lsa.idx_to_doc_map.txt", sep="/"),
                         col.names="name")

transformed.doc.matrix.dist <- read.table(paste(dir,
                                                "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt",
                                                sep="/"),
                                          col.names=idx.to.doc$name)
transformed.doc.matrix.dist <- as.dist(transformed.doc.matrix.dist)
transformed.doc.matrix.cluster <- hclust(transformed.doc.matrix.dist,
                                         method="average")

rows.are.terms <- as.matrix(read.table(paste(dir,
                                             "all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.rows_are_terms.txt",
                                             sep="/")))



## Top left panel
p <- as.phylo(transformed.doc.matrix.cluster)
colors <- unlist(lapply(p$tip.label, color.of))

## Top right panel
var.explained <- sing.vals ^ 2 / sum(sing.vals ^ 2) * 100
sum.of.var <- c();for(i in 1:length(var.explained)) { sum.of.var <- c(sum.of.var, sum(var.explained[1:i]))}


## Bottom left panel
x <- transformed.doc.matrix[, 1:2]
y <- rows.are.terms[, 1:2]
y.mag <- apply(y, MARGIN=2, FUN=function(x){ abs(x) })
top.5.terms.topic.1 <- y.mag[order(y.mag[,1], decreasing=T),][1:top.n,]
top.5.terms.topic.2 <- y.mag[order(y.mag[,2], decreasing=T),][1:top.n,]
top.5.terms <- rbind(top.5.terms.topic.1, top.5.terms.topic.2)


## Bottom right panel
fig2.x <- transformed.doc.matrix[, 3:4]
fig2.y <- rows.are.terms[, 3:4]
fig2.y.mag <- apply(fig2.y, MARGIN=2, FUN=function(x){ abs(x) })
fig2.top.5.terms.topic.3 <- fig2.y.mag[order(fig2.y.mag[,1], decreasing=T),][1:top.n,]
fig2.top.5.terms.topic.4 <- fig2.y.mag[order(fig2.y.mag[,2], decreasing=T),][1:top.n,]
fig2.top.5.terms <- rbind(fig2.top.5.terms.topic.3, fig2.top.5.terms.topic.4)
## swap the columns of the matrices so that topic 2 is on the y axis both times
## fig2.top.5.terms.rev <- cbind(fig2.top.5.terms[,2], fig2.top.5.terms[,1])
## fig2.x.rev <- cbind(fig2.x[,2], fig2.x[,1])


transparent.red <- rgb(0.5, 0, 0, 0.5)
transparent.blue <- rgb(0, 0, 0.5, 0.5)
arrow.colors <- c(
    rep(transparent.blue,  times=top.n),
    rep(transparent.red, times=top.n)
    )

## Plots all down here for easy printing
hw <- 15
pdf(paste(dir, "fig_1ish.3.pdf", sep="/"), height=hw, width=hw)
par(mfrow=c(2,2))
## PLOT
plot(p,
     label.offset = 0.01,
     tip.color=colors,
     direction = "rightwards",
     align.tip.label = FALSE,
     main = "LSA hierarchical clustering",
     frame.plot=F)
## PLOT
plot(sum.of.var, type="o", xlab="Topic", ylab="% variance expained", main="Cumulative variance explained")
## PLOT
biplot.default.2(x,
                 top.5.terms,
                 main="Gene weights, samples in topic space",
                 xlab="Topic 1",
                 ylab="Topic 2",
                 x.plot.col=colors,
                 arrow.col=arrow.colors,
                 cex=1.25
                 )
## PLOT
biplot.default.2(fig2.x,
                 fig2.top.5.terms,
                 main="Gene weights, samples in topic space",
                 xlab="Topic 3",
                 ylab="Topic 4",
                 x.plot.col=colors,
                 arrow.col=arrow.colors,
                 cex=1.25,
                 )
dev.off()
