## compare with aai

par(mfrow=c(1,1))

dat <- read.table("/Users/moorer/Documents/Lab Notebook/LatentSemanticAnalysis/2017 july preprint/reference genomes/mga_orfs/aai/aai_scores/dist.mat", header=T, sep="\t", col.names=c("E. coli", "M. maripaludis", "M. mazei", "P. marinus", "S. elongatus", "S. flexneri plasmid", "S. flexneri"))

dist <- as.dist(1-as.matrix(dat))
cluster <- hclust(dist, method="complete")
cluster$labels <- c("E. coli", "M. maripaludis", "M. mazei", "P. marinus", "S. elongatus", "S. flexneri plasmid", "S. flexneri")

plot(cluster, main="AAI clustring")
