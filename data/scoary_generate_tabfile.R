library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least 4 argumentS: if not, return an error
if (length(args)!=4) {
  length(args)
  stop("Exactly 4 argument must be supplied (input file).n", call.=FALSE)
} 

scoary_data <- read.csv(args[1]) %>% 
  rename(Filename=args[2],Source.Host=args[3]) %>% 
  add_column(Present=1) %>%
  pivot_wider(names_from = Source.Host, values_from = Present, values_fill = 0) %>%  
  mutate(Filename_base=str_remove(Filename, paste(".",args[4],sep=""))) %>% 
  select(!Filename) %>%
  column_to_rownames('Filename_base')
write.table(scoary_data,file="./scoary_traitfile.csv", quote=FALSE,sep=",",row.names = TRUE,col.names=NA)
