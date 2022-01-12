library(tidyverse)
library(data.table)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=6) {
  length(args)
  stop("Exactly 6 arguments must be supplied (input file).n", call.=FALSE)
} 

##########################  Load Data #################################
# panaroo everything
panaroo_all <- read.table(args[1], header=TRUE)
##### Scoary input #####
scoary_sig<- read.csv(args[2], header = TRUE, sep=",", stringsAsFactors = FALSE) %>%  filter(Bonferroni_p < as.double(args[6]))
# Non-clonal metadata
metadata_nonclonal  <- read.csv(args[3], header=TRUE) %>%
  rename(Filename=args[4],Source.Host=args[5]) %>%
  select("Filename","Source.Host")

####################################################################### Collapsed #######################################################
# Combine scoary sig proteins (no duplicates)
sig_PV_c <-unique(c(scoary_sig$Gene ))

# Get first member in collapsed cluster as representative
sig_PV_c_rep <- as.data.frame(sig_PV_c) %>% 
  separate(sig_PV_c, into = c("V1", "Extra"), sep="--",extra = "merge") %>% 
  select(-Extra)

# Filter sig from piggy results
panaroo_sig_c <- panaroo_all %>% filter(Gene %in% sig_PV_c_rep$V1)

# Transpose
panaroo_sig_trans_c <- data.table::transpose(panaroo_sig_c)
colnames(panaroo_sig_trans_c) <- panaroo_sig_c$Gene # set IGRs as column names
panaroo_sig_trans_c$Assembly <- colnames(panaroo_sig_c) # get assemblies to column
panaroo_sig_trans_c = panaroo_sig_trans_c[-1,] # remove header row

###### Combine IGRs with metadata and save #######
# non-clonal 
panaroo_sig_c_met_nonclonal <- inner_join(panaroo_sig_trans_c, metadata_nonclonal, by=c("Assembly" = "Filename")) %>%
  column_to_rownames(var="Assembly")
write.table(panaroo_sig_c_met_nonclonal, file="pv_all.tsv", sep="\t")

# non-clonal & livestock
panaroo_sig_c_met_nonclonal_bps <- panaroo_sig_c_met_nonclonal %>%
  filter(Source.Host != "Human")
write.table(panaroo_sig_c_met_nonclonal_bps, file="pv_bps.tsv", sep="\t")

# human vs non human
panaroo_sig_c_met_nonclonal_human <- panaroo_sig_c_met_nonclonal  %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))
write.table(panaroo_sig_c_met_nonclonal_human, file="pv_human.tsv", sep="\t")
