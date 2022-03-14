library(tidyverse)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=1) {
  length(args)
  stop("Exactly 1 arguments must be supplied (input file).n", call.=FALSE)
} 

#### Load AMR & metadata ####
amr_all <- read.table(args[1], sep="\t", quote="",header=TRUE) %>% 
  filter(Element.type =="AMR")

print("reading files ok")
################## AMR ANTIBIOTIC GENE DATASET #############################
# All hosts, all assemblies 
amr_gene_metadata_all <- amr_all  %>% 
  count(Gene.symbol, Filename) %>% 
  pivot_wider(names_from = Gene.symbol, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
write.table(amr_gene_metadata_all, file="./amr_gene_all.tsv", sep="\t")

################## AMR ANTIBIOTIC CLASS DATASET #############################
# All hosts, all assemblies 
amr_class_metadata_all <- amr_all %>% 
  count(Class, Filename) %>% 
  pivot_wider(names_from = Class, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
write.table(amr_class_metadata_all, file="./amr_class_all.tsv", sep="\t")
