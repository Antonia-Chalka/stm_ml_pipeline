library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=6) {
  length(args)
  stop("Exactly 6 arguments must be supplied (input file).n", call.=FALSE)
} 

########### Load data #############
# Piggy Everything
piggy_all<- read.table(args[1], header=TRUE)
  ##### Scoary input #####
scoary_sig<- read.csv(args[2], header = TRUE, sep=",", stringsAsFactors = FALSE) %>%  filter(Bonferroni_p < as.double(args[6]))
# Non-clonal metadata
metadata_nonclonal  <- read.csv(args[3], header=TRUE) %>%
  rename(Filename=args[4],Source.Host=args[5]) %>%
  select("Filename","Source.Host")

####################################################################### Collapsed #######################################################
# Combine scoary sig proteins (no duplicates)
sig_IGRs_c <-unique(c(scoary_sig$Gene ))

# Get first member in collapsed cluster as representative
sig_IGRs_c_rep <- as.data.frame(sig_IGRs_c) %>% 
  separate(sig_IGRs_c, into = c("V1", "Extra"), sep="--",extra = "merge") %>% 
  select(-Extra)

# Filter sig from piggy results
piggy_sig_c <- piggy_all %>% filter(Gene %in% sig_IGRs_c_rep$V1) 

# Transpose
piggy_sig_trans_c <- data.table::transpose(piggy_sig_c)
colnames(piggy_sig_trans_c) <- piggy_sig_c$Gene # set IGRs as column names
piggy_sig_trans_c$Assembly <- colnames(piggy_sig_c) # get assemblies to column
piggy_sig_trans_c = piggy_sig_trans_c[-1,] # remove header row

###### Combine IGRs with metadata and save #######
# non-clonal
piggy_sig_met_nonclonal_c <- inner_join(piggy_sig_trans_c, metadata_nonclonal, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly")
write.table(piggy_sig_met_nonclonal_c, file="igr_all.tsv", sep="\t")

# non-clonal & livestock
piggy_sig_met_nonclonal_bps_c <- inner_join(piggy_sig_trans_c, metadata_nonclonal, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly") %>%
  filter(Source.Host != "Human")
write.table(piggy_sig_met_nonclonal_bps_c, file="igr_bps.tsv", sep="\t")
