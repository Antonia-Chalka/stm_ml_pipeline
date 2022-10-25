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
  select(Filename,Source.Host) %>%
  mutate(Threshold=ifelse(Source.Host>= as.double(args[4]),1,0)) %>% #todo change to input arg
  column_to_rownames('Filename') %>%
  select(Threshold)
write.table(scoary_data,file="./scoary_traitfile.csv", quote=FALSE,sep=",",row.names = TRUE,col.names=NA)
