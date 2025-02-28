library(tidyverse)

important_pv <- colnames( read.table("in/model_input/pv_all.tsv", header=TRUE,sep="\t")  %>% select(!Source.Host))
important_igr <- colnames(read.table("in/model_input/igr_all.tsv", header=TRUE,sep="\t") %>% select(!Source.Host))

igr_pv <- read.table("in/roary_piggy_combined.tab",quote="",sep="\t") %>% separate(V2, sep="(_\\+_\\+_)", c("gene_1","gene_2")) %>% 
  mutate(gene_1_promoted=if_else(substring(gene_1,1,1)=="*","Y","N" ),
         gene_1=str_replace(gene_1,"\\*",""),
         gene_2_promoted=if_else(substring(gene_2,1,1)=="*","Y","N" ),
         gene_2=str_replace(gene_2,"\\*","")) %>%
  mutate(gene_1=str_replace_all(gene_1,"~","."),
         gene_2=str_replace_all(gene_2,"~","."))

important_igr_with_pv_data <- igr_pv %>%  filter(V1 %in% important_igr)

pv_associated_with_sig_igr <- important_igr_with_pv_data$gene_1 
pv_associated_with_sig_igr <- append(pv_associated_with_sig_igr, important_igr_with_pv_data$gene_2)
pv_associated_with_sig_igr <- unique(pv_associated_with_sig_igr)

captured_pvs <- intersect(important_pv, pv_associated_with_sig_igr)
noncapturedpv <- setdiff(important_pv,pv_associated_with_sig_igr)

igr_with_1_sig_pv <- important_igr_with_pv_data %>% filter(gene_1 %in% important_pv | gene_2 %in% important_pv)
igr_with_2_sig_pv <- important_igr_with_pv_data %>% filter(gene_1 %in% important_pv & gene_2 %in% important_pv)


switched <- read.table("in/switched_region_divergences.csv",quote="",sep=",", header=TRUE) %>% 
  separate(Gene, sep="(_\\+_\\+_)", c("gene_1","gene_2","igr_1","igr_2")) %>%
  mutate(igr_2=str_replace(igr_2,"_aligned","")) %>%
  mutate(gene_1=str_replace(gene_1,"\\*",""),
         gene_2=str_replace(gene_2,"\\*","")) %>%
  mutate(gene_1=str_replace_all(gene_1,"~","."),
         gene_2=str_replace_all(gene_2,"~","."))

sig_switched_igr <- switched %>% filter(igr_1 %in% important_igr | igr_2 %in% important_igr)

pvs_associated_with_sig_switched_igr <- sig_switched_igr %>% filter(gene_1 %in% important_pv | gene_2 %in% important_pv)

