library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  length(args)
  stop("Exactly 2 arguments must be supplied (input file).n", call.=FALSE)
}

# arg[1] core.tab
#arg2 core.ref.tab

#### Load snippy data ####
snippy_in <- read.table(args[1], sep="\t", header=TRUE) %>%
  unite("CHR_POS", CHR:POS, remove = TRUE) 

snippy_core_ref<- read.table(args[2], sep="\t", header=TRUE) %>% #todo change to args[5]
  unite("CHR_POS", CHR:POS, remove = TRUE)

############################################### Adding missing SNPs ###############################################

snippy_all <- right_join(snippy_in,snippy_core_ref)

# Transpose
snp_all_trans <- data.table::transpose(snippy_all)  
colnames(snp_all_trans) <- snippy_all$CHR_POS # set snps as column names
snp_all_trans$Assembly <- colnames(snippy_all) # get assemblies to column
snp_all_trans = snp_all_trans[-1,] # remove header row

snp_all_trans <- snp_all_trans %>% 
  fill(everything(),  .direction = "downup") %>%
  filter(!Assembly=="REF") %>%
  column_to_rownames(var="Assembly")


########### Combine with metadata and output ###########
write.table(snp_all_trans, file="snp_abudance_all.tsv", sep="\t")
