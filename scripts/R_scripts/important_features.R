library(tidyverse)
library(C50)
library(tidyrules)
library(caret)
library(pagoo)

# load models of interest
pv_all_model <- readRDS("./in/models/pv_all_model.rds")
pv_bps_model <- readRDS("./in/models/pv_bps_model.rds")
pv_human_model <- readRDS("./in/models/pv_human_model.rds")

igr_all_model <- readRDS("./in/models/igr_all_model.rds")
igr_bps_model <- readRDS("./in/models/igr_bps_model.rds")
igr_human_model <- readRDS("./in/models/igr_human_model.rds")

# TODO LOAD SNP MODELS
#snp_abudance_all_model <- readRDS("./in/models/snp_abudance_all_model.rds")
#snp_abudance_bps_model <- readRDS("./in/models/snp_abudance_bps_model.rds")
#snp_abudance_human_model <- readRDS("./in/models/snp_abudance_human_model.rds")

importance_model_list <-  list( 
  pv_all_model=pv_all_model, pv_bps_model=pv_bps_model, pv_human_model=pv_human_model,
  igr_all_model=igr_all_model, igr_bps_model=igr_bps_model, igr_human_model=igr_human_model
  #snp_abudance_all_model=snp_abudance_all_model, snp_abudance_bps_model=snp_abudance_bps_model, snp_abudance_human_model=snp_abudance_human_model
  )

dir.create("./out/importance")
dir.create("./out/tuning")

i=1
important_features_all <- data.frame(model=character(), Feature=character(), importance=numeric())

for (model in importance_model_list) {
  # Plot hyper-parameter tuning plot
  png(filename=paste0("./out/tuning/", names(importance_model_list[i]), "_hyperparameter_plot.png"), res=300, units="cm", width = 15, height=15)
  print(plot(model))
  dev.off()
  
  # Plot distribution of importance scores
  important_features <- C5imp(model$finalModel) %>%
    rownames_to_column(var="Feature")
  important_features$model <- names(importance_model_list[i])
  important_features_all <- rbind(important_features, important_features_all)
    
  # Check if rule based or tree based and extract importance list
  if (model$bestTune[1,2] == "rules") {
    rules <- tidyRules(model$finalModel)
    write.table(rules, file=paste0("./out/importance/", names(importance_model_list[i]),"_rules.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
    
  } else if (model$bestTune[1,2]=="tree") {
    #png(file=paste0("./out/importance/", names(importance_model_list[i]), "_treeplot.png"), res=300, units="cm", width = 15, height=15)
    #plot(model$finalModel, type="s", main="Decision Tree")
    #dev.off()
  } else {
    warning("Something has gone wrong. Model is neither tree or rules")
    break 
  }
  i=i+1 
}
write.table(important_features_all %>% filter( human >= 10| bps >= 10 | all >= 10), file=paste0("./out/importance/", names(importance_model_list[i]),"_importance.tsv"), sep="\t", row.names = FALSE, quote = FALSE)

############################# pv #####

pg <- roary_2_pagoo(gene_presence_absence_csv = "in/predict_all/gene_presence_absence.csv")
cluster_info_pv <-as.data.frame(pg$clusters) %>%
  mutate(cluster=str_replace_all(cluster, "~","."))

pv_important_all <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="pv")
ggplot(pv_important_all, aes(x=Overall, fill=model_type)) +
  geom_density(alpha=0.4)+
  theme_bw() +
  labs(title="Importance of all PVs across the three model types", x="Importance ('Overall')", y="% of features", fill="Model Type")

pv_important_contributing <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="pv", Overall>0)
pv_important_contributing_plot <- ggplot(pv_important_contributing, aes(x=Overall, fill=model_type)) +
  geom_density(alpha=0.4)+
  theme_bw() +
  labs(title="Importance of important PVs across the three model types", x="Importance ('Overall>0')", y="% of features with Overall>0", fill="Model Type")
ggsave(pv_important_contributing_plot, file="out/importance/pv_important_contributing.png",  device = "png",dpi = 300, width = 15,height=10, units="cm" )

pv_important_best <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="pv" & Overall>13) %>%
  select(!model & !feature_type)%>%
  pivot_wider(names_from = "model_type", values_from = "Overall") %>%
  left_join(cluster_info_pv, by=c("Feature"="cluster"))
write.csv(pv_important_best, file = "out/importance/pv_important_top.csv", row.names = FALSE,na="" )


pv_all <- read.table("in/model_input/pv_all.tsv",header=TRUE, row.names = 1)
pv_important_best_hostcounts <- data.frame() 
for (feature in pv_important_best$Feature){
  pv_important_best_hostcounts_single <- pv_all %>% 
    count(Source.Host, !!as.name(feature)) %>% 
    group_by(Source.Host) %>% 
    mutate(perc = (n/sum(n))*100, PV= feature) %>%
    dplyr::rename(Presence = !!as.name(feature))
  pv_important_best_hostcounts <- rbind(pv_important_best_hostcounts, pv_important_best_hostcounts_single)
  
}
write.csv(pv_important_best_hostcounts, file = "out/importance/pv_important_best_hostcounts.csv", row.names = FALSE,na="" )

 

#################### PV Correlation


test <- pv_all %>% select(!Source.Host)%>%
  as.matrix %>%
  cor %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1)


gene_corr <- pv_all %>% select(!Source.Host)
gene_corr <- cor(gene_corr[ , colnames(gene_corr) != "ybaL...fsr...kefB"], gene_corr$ybaL...fsr...kefB) %>%
  as.data.frame() %>%  rownames_to_column()








################# igr ####
  
igr <- roary_2_pagoo(gene_presence_absence_csv = "in/predict_all/IGR_presence_absence.csv")
cluster_info_igr <-as.data.frame(igr$clusters) %>%
  mutate(cluster=str_replace_all(cluster, "~","."))
igr_important_all <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="igr")
ggplot(igr_important_all, aes(x=Overall, fill=model_type)) +
  geom_density(alpha=0.4)+
  theme_bw() +
  labs(title="Importance of all IGRs across the three model types", x="Importance ('Overall')", y="% of features", fill="Model Type")
igr_important_contributing <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="igr", Overall>0)
igr_important_contributing_plot <- ggplot(igr_important_contributing, aes(x=Overall, fill=model_type)) +
  geom_density(alpha=0.4)+
  theme_bw() +
  labs(title="Importance of important IGRs across the three model types", x="Importance ('Overall>0')", y="% of features with Overall>0", fill="Model Type")
ggsave(igr_important_contributing_plot, file="out/importance/igr_important_contributing_plot.png",  device = "png",dpi = 300, width = 15,height=10, units="cm" )

igr_important_best <- important_features_all %>% 
  mutate(model_type=str_extract(model, "all|bps|human"), feature_type=str_extract(model,"igr|pv")) %>%
  filter(feature_type=="igr" & Overall>13) %>%
  select(!model & !feature_type)%>%
  pivot_wider(names_from = "model_type", values_from = "Overall") %>%
  left_join(cluster_info_igr, by=c("Feature"="cluster"))

write.csv(igr_important_best, file = "out/importance/igr_important_best.csv", row.names = FALSE,na="" )


igr_all <- read.table("in/model_input/igr_all.tsv",header=TRUE, row.names = 1)
igr_important_best_hostcounts <- data.frame()
for (feature in igr_important_best$Feature){
  igr_important_best_hostcounts_single <- igr_all %>% 
    count(Source.Host, !!as.name(feature)) %>% 
    group_by(Source.Host) %>% 
    mutate(perc = (n/sum(n))*100, IGR= feature) %>%
    dplyr::rename(Presence = !!as.name(feature))
  igr_important_best_hostcounts <- rbind(igr_important_best_hostcounts, igr_important_best_hostcounts_single)
  
}
write.csv(igr_important_best_hostcounts, file = "out/importance/igr_important_best_hostcounts.csv", row.names = FALSE,na="" )

  

  
 
#test <- rules %>% filter(str_detect(LHS,"ccmB")) %>% mutate(rule=str_extract(LHS,"ccmB (=|>|<|>=|<|<=|<>) [01]")) %>% count(rule,RHS)

 





