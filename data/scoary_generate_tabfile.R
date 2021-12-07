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
  rename(Filename=args[2],Source.Host=args[3])

scoary_data <- dcast(setDT(scoary_data), Filename~Source.Host, fun=length, drop=c(TRUE,FALSE))

scoary_data <- scoary_data %>%  
  mutate(Filename_base=str_remove(Filename, args[4])) %>% 
  select(!Filename) %>%
  column_to_rownames('Filename_base')
write.table(scoary_data,file="./scoary_traitfile.csv", quote=FALSE,sep=",",row.names = TRUE,col.names=NA)
