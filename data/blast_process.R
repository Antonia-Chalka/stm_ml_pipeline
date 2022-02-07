library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=5) {
  length(args)
  stop("Exactly 5 arguments must be supplied (input file).n", call.=FALSE)
} 

blast <- read.table(args[1],sep="\t",quote="")
metadata <- read.csv(args[2], header=TRUE) %>%
  rename(Filename=args[3],Source.Host=args[4]) %>% 
  select("Filename","Source.Host")

blast_pv <- blast %>% 
  filter(V11<0.01) %>%
  distinct(V1,V2) %>%
  pivot_wider(names_from = V1, values_from = V1,
              values_fn = list(V1 = length), values_fill = list(V1 = 0)) %>%
  mutate(Assembly = paste(V2, "fasta",sep=".")) %>% # TODO CHange "fasta" to param filename
  select(!V2) %>%
  left_join(metadata, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly")

write.table(blast_pv, file=args[5], sep="\t")
