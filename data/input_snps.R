library(tidyverse)

set.seed(100)

#### Load snippy data ####
snippy_all <- read.table("./data/raw/core.tab", sep="\t", header=TRUE) %>%
  select(-SAL_LA7644AA_AS.scaffold)

# load metadata
metadata_all <- read.csv("./data/raw/metadata/metadata_quast_all_merged_filtered_discarded_modified.csv") %>%
  select("Filename","Source.Host")

metadata_nonclonal <- read.csv("./data/raw/metadata/metadata_quast_all_merged_filtered_discarded_modified_noclonal.csv") %>% select(!X) %>%
  select("Filename","Source.Host")


############################################### Abudance filtering ###############################################


# Count SNP abudance
snp_all_counts <-  snippy_all %>% 
  pivot_longer(-c("CHR","POS")) %>% 
  filter(name!="REF") %>% 
  count(CHR, POS, value) 
#66, 730 snps


# Exclude SNPs present in <1% or >99% of genomes

# calculate thresholsd (num of assemblies)
abudance_lower <- round(nrow(metadata_all)/100 *1)
abudance_upper <- round(nrow(metadata_all)/100 *99)


# filter for >1 & <99
snp_all_include <- snp_all_counts %>%
  filter(n > abudance_lower & n < abudance_upper)
# Above will filter single snps bases
# eg could be 3278 19 16

# get position and chromosome (unique id) ( + remove snps with < 1 variant)
snp_include_list <- snp_all_include %>%  select(CHR,POS) %>%
  count(CHR,POS) %>%
  filter(n > 1) %>%
  select(-n)
# 1546 snp final


# filter snps and get metadata & combine chr and pos
snp_filter <-  inner_join(snippy_all, snp_include_list, by = c("CHR" = "CHR", "POS" = "POS")) %>%
  unite("CHR_POS", CHR:POS, remove = TRUE) %>%
  select(-REF)

# Transpose
snp_filter_trans <- data.table::transpose(snp_filter)  
colnames(snp_filter_trans) <- snp_filter$CHR_POS # set snps as column names
snp_filter_trans$Assembly <- colnames(snp_filter) # get assemblies to column
snp_filter_trans = snp_filter_trans[-1,] # remove header row

########### Combine with metadata and output ###########

# All
snp_abudance_all <- inner_join(snp_filter_trans, metadata_all, by=c("Assembly" = "Filename")) %>%
  column_to_rownames(var="Assembly")
write.table(snp_abudance_all, file="./data/out/round_2/input_data/snp_abudance_all.tsv", sep="\t")

# Non clonal
snp_abudance_nonclonal <- inner_join(snp_filter_trans, metadata_nonclonal, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly")
write.table(snp_abudance_nonclonal, file="./data/out/round_2/input_data/snp_abudance_nonclonal.tsv", sep="\t")

# BPS
snp_abudance_bps <- inner_join(snp_filter_trans, metadata_all, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly") %>%
  filter(Source.Host != "Human")
write.table(snp_abudance_bps, file="./data/out/round_2/input_data/snp_abudance_bps.tsv", sep="\t")

# Non clonal bps
snp_abudance_nonclonal_bps <- inner_join(snp_filter_trans, metadata_nonclonal, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly") %>%
  filter(Source.Host != "Human")
write.table(snp_abudance_nonclonal_bps, file="./data/out/round_2/input_data/snp_abudance_nonclonal_bps.tsv", sep="\t")




