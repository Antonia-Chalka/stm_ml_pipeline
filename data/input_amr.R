library(tidyverse)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=4) {
  length(args)
  stop("Exactly 4 arguments must be supplied (input file).n", call.=FALSE)
} 

#### Load AMR & metadata ####
amr_all <- read.table(args[1], sep="\t", quote="",header=TRUE) %>% 
  filter(Element.type =="AMR")
metadata_all <- read.csv(args[2]) %>% 
  rename(Filename=args[2],Source.Host=args[3])

print("reading files ok")
################## AMR ANTIBIOTIC GENE DATASET #############################
# All hosts, all assemblies 
amr_gene_metadata_all <- inner_join(amr_all, metadata_all, by="Filename") %>% 
  count(Gene.symbol, Filename, Source.Host) %>% 
  pivot_wider(names_from = Gene.symbol, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
write.table(amr_gene_metadata_all, file="./amr_gene_all.tsv", sep="\t")

# Only livestock assemblies 
amr_gene_metadata_bps <- inner_join(amr_all, metadata_all, by="Filename") %>% 
  count(Gene.symbol, Filename, Source.Host) %>%
  filter(Source.Host != "Human") %>%
  pivot_wider(names_from = Gene.symbol, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")

write.table(amr_gene_metadata_bps, file="./amr_gene_bps.tsv", sep="\t")

################## AMR ANTIBIOTIC CLASS DATASET #############################
# All hosts, all assemblies 
amr_class_metadata_all <- inner_join(amr_all, metadata_all, by="Filename") %>% 
  count(Class, Filename, Source.Host) %>% 
  pivot_wider(names_from = Class, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
write.table(amr_class_metadata_all, file="./amr_class_all.tsv", sep="\t")

# Only livestock assemblies 
amr_class_metadata_bps <- inner_join(amr_all, metadata_all, by="Filename") %>% 
  count(Class, Filename, Source.Host) %>% 
  filter(Source.Host != "Human") %>%
  pivot_wider(names_from = Class, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
write.table(amr_class_metadata_bps, file="./amr_class_bps.tsv", sep="\t")
