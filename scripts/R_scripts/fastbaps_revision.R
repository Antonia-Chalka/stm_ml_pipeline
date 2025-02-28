library(fastbaps)
library(phytools)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ggnewscale)

set.seed(100)

iqtree <- treeio::read.newick("in/gubbins.node_labelled.final_tree.tre")
iqtree <- drop.tip(iqtree, c("Reference"))
iqtree.rooted <- midpoint.root(iqtree)

sparse.data <- import_fasta_sparse_nt("in/core.aln") #TODO set as params
sparse.data <- optimise_prior(sparse.data, type = "baps")
baps.hc <- fast_baps(sparse.data)

plot.df <- data.frame(clusters = best_baps_partition(sparse.data, as.phylo(baps.hc)), stringsAsFactors = FALSE) %>% rownames_to_column(var="assembly")
assembly_names <- data.frame(assembly = iqtree.rooted$tip.label)
clusters_clean <- right_join(plot.df,assembly_names) 

write.table(clusters_clean, file = "./out/fastbaps_clusters.tsv", sep="\t", quote = FALSE, row.names = FALSE)
