library("tidyverse")


metadata_all <- read.csv("in/predict_all/good_noext_metadata.csv")
panaroo_all <- read.table("in_old/predict_all/gene_presence_absence.Rtab", header=TRUE)  

panaroo_sig_c <- panaroo_all %>% filter(Gene %in% c("gsk","ybaL","ybaL~~~fsr~~~kefB","group_3075","group_2854","group_1910")) 
panaroo_sig_trans_c <- as.data.frame(t(panaroo_sig_c))
colnames(panaroo_sig_trans_c) <- panaroo_sig_c$Gene # set IGRs as column names
panaroo_sig_trans_c <- panaroo_sig_trans_c %>% mutate_if(is.character,as.numeric)
panaroo_sig_trans_c$Assembly <- colnames(panaroo_sig_c) # get assemblies to column
panaroo_sig_trans_c = panaroo_sig_trans_c[-1,] # remove header row

pv_all<- inner_join(panaroo_sig_trans_c, metadata_all, by=c("Assembly" = "Filename")) %>%
  select(!Collection.Year & !Region)

pv_counts <-pv_all %>% 
  pivot_longer(cols = !Source.Host & !Assembly) %>%
  count(Source.Host, name, value) %>% 
  group_by(Source.Host,name) %>% 
  mutate(perc = (n/sum(n))*100) 
  
pv_unique <- pv_all %>% select(!Assembly) %>% 
  count(group_1910, `ybaL~~~fsr~~~kefB`,group_2854,ybaL,gsk,group_3075, Source.Host) %>% 
  relocate(group_1910, `ybaL~~~fsr~~~kefB`,group_2854,ybaL,gsk,group_3075, Source.Host) %>%
  pivot_wider( names_from="Source.Host", values_from="n")
