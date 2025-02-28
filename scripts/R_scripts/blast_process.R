library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=2) {
  length(args)
  stop("Exactly 2 arguments must be supplied (input file).n", call.=FALSE)
} 

# arg1 igr_filter.fasta_results.tsv  pv_filter.fasta_results.tsv
# arg2  igr_blast_all.tsv

blast <- read.table(args[1],sep="\t",quote="")

blast_pv <- blast %>% 
  filter(V11<0.01) %>%
  distinct(V1,V2) %>%
  pivot_wider(names_from = V1, values_from = V1,
              values_fn = list(V1 = length), values_fill = list(V1 = 0)) %>%
  rename(Assembly = V2) %>%
  column_to_rownames(var="Assembly")

write.table(blast_pv, file=paste(args[2], "_blast.tsv" ,sep=""), sep="\t")
