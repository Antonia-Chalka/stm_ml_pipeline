library(tidyverse)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  length(args)
  stop("Exactly 4 arguments must be supplied (input file).n", call.=FALSE)
}

#### Load snippy data ####
snippy_all <- read.table(args[1], sep="\t", header=TRUE)

# load metadata
metadata_nonclonal <- read.csv(args[2]) %>%
  rename(Filename=args[3],Source.Host=args[4]) %>%
  select("Filename","Source.Host")

############################################### Abudance filtering ###############################################
# Exclude SNPs present in <1% or >99% of genomes

# Count SNP abudance
snp_all_counts <-  snippy_all %>% 
  pivot_longer(-c("CHR","POS")) %>% 
  filter(name!="REF") %>% 
  count(CHR, POS, value) 
# calculate thresholsd (num of assemblies)
abudance_lower <- round(nrow(metadata_nonclonal)/100 *1)
abudance_upper <- round(nrow(metadata_nonclonal)/100 *99)
# filter for >1 & <99
snp_all_include <- snp_all_counts %>%
  filter(n > abudance_lower & n < abudance_upper)
# get position and chromosome (unique id) ( + remove snps with < 1 variant)
snp_include_list <- snp_all_include %>%  select(CHR,POS) %>%
  count(CHR,POS) %>%
  filter(n > 1) %>%
  select(-n)
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
# Non clonal
snp_abudance_nonclonal <- inner_join(snp_filter_trans, metadata_nonclonal, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly")
write.table(snp_abudance_nonclonal, file="snp_abudance_all.tsv", sep="\t")

# Non clonal bps
snp_abudance_nonclonal_bps <- snp_abudance_nonclonal %>%
  filter(Source.Host != "Human")
write.table(snp_abudance_nonclonal_bps, file="snp_abudance_bps.tsv", sep="\t")

# Non clonal bps
snp_abudance_nonclonal_human <- snp_abudance_nonclonal %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))
write.table(snp_abudance_nonclonal_human, file="snp_abudance_human.tsv", sep="\t")
