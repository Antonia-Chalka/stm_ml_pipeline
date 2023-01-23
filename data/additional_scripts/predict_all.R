library(tidyverse)
library(randomForest)
library(caret)
library(pROC) # VESA - DONT NEED TO INSTALL THAT, FEEL FREE TO DELETE LINE
set.seed(100)


########################################################################## Load Data #####################################################
# Make list of models
metadata_all <- read.csv("in/predict_all/good_noext_metadata.csv") # VESA - CHANGE THIS TO  ./model_out/1.assembly_quality/good_noext_metadata.csv

##### AMR Gene #####
amr_all <- read.table("in/predict_all/amr_all.tsv", sep="\t", quote="",header=TRUE)%>%  # VESA - CHANGE TO ./model_out/4.model/model_input/amr_all.tsv
  filter(Element.type =="AMR")

amr_gene_all <- full_join(amr_all, metadata_all, by="Filename") %>% 
  count(Gene.symbol, Filename, Source.Host) %>% 
  pivot_wider(names_from = Gene.symbol, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")

amr_gene_all[is.na(amr_gene_all)] <- 0

names(amr_gene_all) <- gsub(x = names(amr_gene_all), pattern = "\\/", replacement = ".") 
names(amr_gene_all) <- gsub(x = names(amr_gene_all), pattern = "\\-", replacement = ".") 
names(amr_gene_all) <- gsub(x = names(amr_gene_all), pattern = "\\)", replacement = ".") 
names(amr_gene_all) <- gsub(x = names(amr_gene_all), pattern = "\\(", replacement = ".") 
names(amr_gene_all) <- gsub(x = names(amr_gene_all), pattern = "\\'", replacement = ".") 

amr_gene_bps <- amr_gene_all %>% filter(!Source.Host=="Human")

amr_gene_human <- amr_gene_all %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))

##### AMR Class #####
amr_class_all <- full_join(amr_all, metadata_all, by="Filename") %>% 
  count(Class, Filename, Source.Host) %>% 
  pivot_wider(names_from = Class, values_from = n,  values_fill = 0) %>%
  column_to_rownames(var="Filename")
names(amr_class_all) <- gsub(x = names(amr_class_all), pattern = "\\/", replacement = ".") 
names(amr_class_all) <- gsub(x = names(amr_class_all), pattern = "\\-", replacement = ".") 

amr_class_bps <- amr_class_all %>% filter(!Source.Host=="Human")

amr_class_human <- amr_class_all %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))

rm(amr_all)

##### IGR #####
piggy_all<- read.table("in/predict_all/IGR_presence_absence.Rtab", header=TRUE)   # VESA - CHANGE TO ./model_out/2.genomic_features/igr_out/piggy_out/IGR_presence_absence.Rtab
scoary_sig<- read.csv("in/predict_all/IGR_scoary_all.csv", header = TRUE, sep=",", stringsAsFactors = FALSE) %>%  filter(Bonferroni_p < as.double(0.05)) # VESA - CHANGE TO ./model_out/2.genomic_features/igr_out/scoary_igr/scoary_all.csv

sig_IGRs_c <-unique(c(scoary_sig$Gene )) # Combine scoary sig proteins (no duplicates)
sig_IGRs_c_rep <- as.data.frame(sig_IGRs_c) %>% # Get first member in collapsed cluster as representative
  separate(sig_IGRs_c, into = c("V1", "Extra"), sep="--",extra = "merge") %>% 
  select(-Extra)
piggy_sig_c <- piggy_all %>% filter(Gene %in% sig_IGRs_c_rep$V1) # Filter sig from piggy results
piggy_sig_trans_c <- data.table::transpose(piggy_sig_c)
colnames(piggy_sig_trans_c) <- piggy_sig_c$Gene # set IGRs as column names
piggy_sig_trans_c <- piggy_sig_trans_c %>% mutate_if(is.character,as.numeric)
piggy_sig_trans_c$Assembly <- colnames(piggy_sig_c) # get assemblies to column
piggy_sig_trans_c = piggy_sig_trans_c[-1,] # remove header row

igr_all <- inner_join(piggy_sig_trans_c, metadata_all, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly") 

igr_bps <- igr_all %>% filter(!Source.Host=="Human")

igr_human <- igr_all %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))

rm(piggy_all, scoary_sig, sig_IGRs_c, sig_IGRs_c_rep, piggy_sig_c, piggy_sig_trans_c)

##### PV #####
panaroo_all <- read.table("in/predict_all/gene_presence_absence.Rtab", header=TRUE)  # VESA - CHANGE TO ./model_out/2.genomic_features/pv_out/panaroo_out/gene_presence_absence.Rtab
scoary_sig<- read.csv("in/predict_all/PV_scoary_all.csv", header = TRUE, sep=",", stringsAsFactors = FALSE) %>%  filter(Bonferroni_p < as.double(1.1)) # VESA - CHANGE TO ./model_out/2.genomic_features/pv_out/scoary_pv/scoary_all.csv

sig_PV_c <-unique(c(scoary_sig$Gene )) # Combine scoary sig proteins (no duplicates)
sig_PV_c_rep <- as.data.frame(sig_PV_c) %>% # Get first member in collapsed cluster as representative
  separate(sig_PV_c, into = c("V1", "Extra"), sep="--",extra = "merge") %>% 
  select(-Extra)
panaroo_sig_c <- panaroo_all %>% filter(Gene %in% sig_PV_c_rep$V1) # Filter sig from piggy results
panaroo_sig_trans_c <- data.table::transpose(panaroo_sig_c)
colnames(panaroo_sig_trans_c) <- panaroo_sig_c$Gene # set IGRs as column names
panaroo_sig_trans_c <- panaroo_sig_trans_c %>% mutate_if(is.character,as.numeric)
panaroo_sig_trans_c$Assembly <- colnames(panaroo_sig_c) # get assemblies to column
panaroo_sig_trans_c = panaroo_sig_trans_c[-1,] # remove header row

pv_all<- inner_join(panaroo_sig_trans_c, metadata_all, by=c("Assembly" = "Filename")) %>%
  column_to_rownames(var="Assembly")
names(pv_all) <- gsub(x = names(pv_all), pattern = "\\~", replacement = ".") 

pv_bps <- pv_all %>% filter(!Source.Host=="Human")

pv_human <- pv_all %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))

rm(panaroo_all, scoary_sig, sig_PV_c, sig_PV_c_rep, panaroo_sig_c, panaroo_sig_trans_c)

##### sNP #####
snippy_all <- read.table("in/predict_all/core.tab", sep="\t", header=TRUE) # VESA - CHANGE TO ./model_out/2.genomic_features/snp_core_out/core.tab

snp_all_counts <-  snippy_all %>% 
  pivot_longer(-c("CHR","POS")) %>% 
  filter(name!="REF") 
snp_all_counts <- data.table::data.table(snp_all_counts)
snp_all_counts<- snp_all_counts[, .N, by=list(CHR, POS, value)] %>% 
  rename(n=N)

# calculate thresholsd (num of assemblies)
abudance_lower <- round(ncol(snippy_all)/100 *0.1)
abudance_upper <- round(ncol(snippy_all)/100 *99)

# filter for >1 & <99
snp_all_include <- snp_all_counts %>%
  filter(n > abudance_lower & n < abudance_upper)

# get position and chromosome (unique id) ( + remove snps with < 1 variant)
snp_include_list <- snp_all_include %>%  select(CHR,POS) %>%
  count(CHR,POS) %>%
  filter(n > 1) %>%
  select(-n)

# filter snps and get metadata & combine chr and pos
snp_filter <-  inner_join(snippy_all, snp_include_list, by = c("CHR" = "CHR", "POS" = "POS")) %>%
  unite("CHR_POS", CHR:POS, remove = TRUE) %>%
  select(-REF)

snp_filter_trans <- data.table::transpose(snp_filter)  
colnames(snp_filter_trans) <- snp_filter$CHR_POS # set snps as column names
snp_filter_trans$Assembly <- colnames(snp_filter) # get assemblies to column
snp_filter_trans = snp_filter_trans[-1,] # remove header row

snp_abudance_all <- inner_join(snp_filter_trans, metadata_all, by=c("Assembly" = "Filename"))  %>%
  column_to_rownames(var="Assembly")

snp_abudance_bps <- snp_abudance_all %>% filter(!Source.Host=="Human")

snp_abudance_human <- snp_abudance_all %>%
  mutate(Source.Host = if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host))

rm(snippy_all, snp_all_counts, abudance_lower, abudance_upper, snp_all_include, snp_include_list, snp_filter, snp_filter_trans)

####################################################################### Prediction #######################################################

# Make list of models
model_all <- list(amr_class_all = amr_class_all,
                  amr_gene_all = amr_gene_all,
                  igr_all = igr_all,
                  pv_all = pv_all,
                  snp_abudance_all = snp_abudance_all
)
model_bps <- list( amr_class_bps = amr_class_bps,
                   amr_gene_bps = amr_gene_bps,
                   pv_bps = pv_bps,
                   igr_bps = igr_bps,
                   snp_abudance_bps = snp_abudance_bps
)
model_human <- list ( amr_class_human = amr_class_human,
                      amr_gene_human = amr_gene_human,
                      pv_human = pv_human,
                      igr_human = igr_human,
                      snp_abudance_human = snp_abudance_human
)

predictions_df <- pv_all %>% select(Source.Host) %>% 
  mutate(Source.Type=if_else(!(Source.Host %in% c("Human")), "Livestock", Source.Host)) %>%
  rownames_to_column(var="Assembly")
thresholds_df <- data.frame(Host=character(),model=character(),optimal_threshold=numeric())
f1_all <- data.frame()
dir.create("./out/roc_plots")

#### Run all assemblies on all types of models ######

# All Host Models
i=1
for (model in model_all) {
  model_filename <- paste("./in/models/", names(model_all)[i], "_model.rds", sep="") #VESA -change first string ("./in/models/") to "./model_out/4.model/models/"
  rf_random <- readRDS(model_filename) 
  
  # Get F1
  predict_class <- predict(rf_random, newdata = model, type = "raw") 
  cm <- confusionMatrix(predict_class, as.factor(model$Source.Host), mode="everything")
  f1 <- as.data.frame(cm$byClass) %>% select(F1) %>% setNames((names(model_all)[i])) %>%
     t() %>% as.data.frame()
  f1_all <- bind_rows(f1_all,f1)
  
  # Get predictions by class
  predict_prob <- predict(rf_random, newdata = model, type="prob") 
  predictions_df %>%
    full_join(rownames_to_column(predict_prob, var="Assembly") )
  
  # Get highest scored prediction
  predict_prob$prediction <- names(predict_prob)[1:4][apply(predict_prob[,1:4], 1, which.max)]
  
  # VESA delete the code below until indicated
  # Get roc scores
  predict_prob$observed <- model$Source.Host
  roc.human <- roc(ifelse(predict_prob$observed=="Human", "Human", "non-human"), as.numeric(predict_prob$Human))
  roc.bovine <- roc(ifelse(predict_prob$observed=="Bovine", "Bovine", "non-bovine"), as.numeric(predict_prob$Bovine))
  roc.poultry <- roc(ifelse(predict_prob$observed=="Poultry", "Poultry", "non-Poultry"), as.numeric(predict_prob$Poultry))
  roc.swine <- roc(ifelse(predict_prob$observed=="Swine", "Swine", "non-Swine"), as.numeric(predict_prob$Swine))
  roc.list <- list(Human = roc.human, Bovine = roc.bovine, Poultry=roc.poultry, Swine=roc.swine)
  
  # Plot auc
  myauc <- paste("Human AUC=",round(auc(roc.human), digits=2), 
                 "\n Bovine AUC=",round(auc(roc.bovine), digits=2),
                 "\n Poultry AUC=",round(auc(roc.poultry), digits=2), 
                 "\n Swine AUC=", round(auc(roc.swine), digits=2), sep='')
  roc_plot <- ggroc(roc.list) +
    geom_abline(intercept = 1, slope = 1) +
    labs(x="Specificity", y="Sensitivity",title=paste0("ROC Curve for the ",names(model_all)[i], " model")) +
    annotate(geom="text", x=0.7, y=0.7, label=myauc) + 
    scale_colour_discrete("Host") +
    theme_bw()
  ggsave(roc_plot, filename=paste0("out/roc_plots/roc_",names(model_all)[i], ".png", sep=''), device = "png", width = 15, height = 15, units="cm")
  
  # Calculate thresholds for optimal prediction
  thresholds_df <- structure(rbind(thresholds_df,c("Human",
                                                   names(model_all)[i],
                                                   coords(roc.human, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  thresholds_df <- structure(rbind(thresholds_df,c("Poultry",
                                                   names(model_all)[i],
                                                   coords(roc.poultry, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  thresholds_df <- structure(rbind(thresholds_df,c("Bovine",
                                                   names(model_all)[i],
                                                   coords(roc.bovine, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  thresholds_df <- structure(rbind(thresholds_df,c("Swine",
                                                   names(model_all)[i],
                                                   coords(roc.swine, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  
  # Change column names to include model name
  names(predict_prob) <- paste0(names(model_all)[i], "_", names(predict_prob))
  
  predict_prob <- predict_prob %>%
    mutate (!!paste0(names(model_all)[i], "_model_prediction_threshold"):=case_when(
      !!(as.name(paste0(names(model_all)[i],"_Human")))>=(coords(roc.human, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Human",
      !!(as.name(paste0(names(model_all)[i],"_Poultry")))>=(coords(roc.poultry, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Poultry",
      !!(as.name(paste0(names(model_all)[i],"_Bovine")))>=(coords(roc.bovine, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Bovine",
      !!(as.name(paste0(names(model_all)[i],"_Swine")))>=(coords(roc.swine, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Swine",
      TRUE ~ 'NA'))
  
  # Add prediction columns and create true/false column about prediction accuracy
  predictions_df <- predictions_df %>%
    full_join(rownames_to_column(predict_prob, var="Assembly") )  %>%
    mutate( !!paste0(names(model_all)[i], "_model_accuracy_highest") := ifelse(!!(as.name(paste0(names(model_all)[i],"_prediction")))==Source.Host, "TRUE", "FALSE"),
            !!paste0(names(model_all)[i], "_model_accuracy_threshold") := ifelse(!!(as.name(paste0(names(model_all)[i],"_model_prediction_threshold")))==Source.Host, "TRUE", "FALSE")
    )
  # VESA stop deleting here
  i=i+1
}

# BPS models
i=1
for (model in model_bps) {
  model_filename <- paste("./in/models/", names(model_bps)[i], "_model.rds", sep="") #VESA -change first string ("./in/models/") to "./model_out/4.model/models/"
  
  rf_random <- readRDS(model_filename) 
  predict_prob <- predict(rf_random, newdata = model, type="prob") 
  
  # Get F1
  predict_class <- predict(rf_random, newdata = model, type = "raw") 
  cm <- confusionMatrix(predict_class, as.factor(model$Source.Host), mode="everything")
  f1 <- as.data.frame(cm$byClass) %>% select(F1) %>% setNames((names(model_all)[i])) %>%
    t() %>% as.data.frame()
  f1_all <- bind_rows(f1_all,f1)
  
  # Get highest scored prediction
  predict_prob$prediction <- names(predict_prob)[1:3][apply(predict_prob[,1:3], 1, which.max)]
  # VESA delete the code below until indicated
  # Get roc scores
  predict_prob$observed <- model$Source.Host
  roc.bovine <- roc(ifelse(predict_prob$observed=="Bovine", "Bovine", "non-bovine"), as.numeric(predict_prob$Bovine))
  roc.poultry <- roc(ifelse(predict_prob$observed=="Poultry", "Poultry", "non-Poultry"), as.numeric(predict_prob$Poultry))
  roc.swine <- roc(ifelse(predict_prob$observed=="Swine", "Swine", "non-Swine"), as.numeric(predict_prob$Swine))
  roc.list <- list(Bovine = roc.bovine, Poultry=roc.poultry, Swine=roc.swine)
  
  # Plot auc
  myauc <- paste("Bovine AUC=",round(auc(roc.bovine), digits=2),
                 "\n Poultry AUC=",round(auc(roc.poultry), digits=2), 
                 "\n Swine AUC=", round(auc(roc.swine), digits=2), sep='')
  roc_plot <- ggroc(roc.list) +
    geom_abline(intercept = 1, slope = 1) +
    labs(x="Specificity", y="Sensitivity", title=paste0("ROC Curve for the ",names(model_bps)[i], " model")) +
    annotate(geom="text", x=0.7, y=0.7, label=myauc) + 
    scale_colour_discrete("Host") +
    theme_bw()
  ggsave(roc_plot, filename=paste0("out/roc_plots/roc_",names(model_bps)[i], ".png", sep=''), device = "png",width = 15, height = 15, units="cm")
  
  # Calculate thresholds for optimal prediction
  thresholds_df <- structure(rbind(thresholds_df,c("Poultry",
                                                   names(model_bps)[i],
                                                   coords(roc.poultry, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  thresholds_df <- structure(rbind(thresholds_df,c("Bovine",
                                                   names(model_bps)[i],
                                                   coords(roc.bovine, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  thresholds_df <- structure(rbind(thresholds_df,c("Swine",
                                                   names(model_bps)[i],
                                                   coords(roc.swine, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))
  
  # Change column names to include model name
  names(predict_prob) <- paste0(names(model_bps)[i], "_", names(predict_prob))
  
  predict_prob <- predict_prob %>%
    mutate (!!paste0(names(model_bps)[i], "_model_prediction_threshold"):=case_when(
    !!(as.name(paste0(names(model_bps)[i],"_Poultry")))>=(coords(roc.poultry, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Poultry",
    !!(as.name(paste0(names(model_bps)[i],"_Bovine")))>=(coords(roc.bovine, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Bovine",
    !!(as.name(paste0(names(model_bps)[i],"_Swine")))>=(coords(roc.swine, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Swine",
    TRUE ~ 'NA'))
  
  # Add prediction columns and create true/false column about prediction accuracy
  predictions_df <- predictions_df %>%
    full_join(rownames_to_column(predict_prob, var="Assembly") ) %>%
    mutate( !!paste0(names(model_bps)[i], "_model_accuracy_highest") := ifelse(!!(as.name(paste0(names(model_bps)[i],"_prediction")))==Source.Host, "TRUE", "FALSE"),
            !!paste0(names(model_bps)[i], "_model_accuracy_threshold") := ifelse(!!(as.name(paste0(names(model_bps)[i],"_model_prediction_threshold")))==Source.Host, "TRUE", "FALSE")
    )
  # VESA stop deleting here
  i=i+1
}

# Human Scoring Models
i=1
for (model in model_human) {
  model_filename <- paste("./in/models/", names(model_human)[i], "_model.rds", sep="") #VESA -change first string ("./in/models/") to "./model_out/4.model/models/"
    #load(model_filename)
    rf_random <- readRDS(model_filename) 
    
    # Get F1
    predict_class <- predict(rf_random, newdata = model, type = "raw") 
    cm <- confusionMatrix(predict_class, as.factor(model$Source.Host), mode="everything")
    f1 <- cm$byClass %>% t() %>% as.data.frame() %>% select(F1) %>% setNames((names(model_all)[i])) %>%
      t() %>% as.data.frame() 
    f1_all <- bind_rows(f1_all,f1)
    
    predict_prob <- predict(rf_random, newdata = model, type="prob") 
    # Get highest scored prediction
    predict_prob$prediction <- names(predict_prob)[1:2][apply(predict_prob[,1:2], 1, which.max)]
    # VESA delete the code below until indicated
    roc.human <- roc(predict_prob$prediction, as.numeric(predict_prob$Human))
    roc.livestock <- roc(predict_prob$prediction, as.numeric(predict_prob$Livestock))
    
    # Plot auc 
    roc_plot <- ggroc(roc.human) +
      geom_abline(intercept = 1, slope = 1) +
      labs(x="Specificity", y="Sensitivity", title=paste0("ROC Curve for the ",names(model_bps)[i], " model")) +
      annotate(geom="text", x=0.7, y=0.7, label=paste("Human AUC=",round(auc(roc.human), digits=2), sep='')) + 
      scale_colour_discrete("Host") +
      theme_bw()
    ggsave(roc_plot, filename=paste0("out/roc_plots/roc_",names(model_human)[i], ".png", sep=''), device = "png",width = 15, height = 15, units="cm")
    
    # Calculate thresholds for optimal prediction
    thresholds_df <- structure(rbind(thresholds_df,c("Human",
                                                     names(model_human)[i],
                                                     coords(roc.human, "best", ret = "threshold", best.method="youden")[1,1])), .Names = names(thresholds_df))

    # Change column names to include model name
    names(predict_prob) <- paste0(names(model_human)[i], "_", names(predict_prob))
    
    predict_prob <- predict_prob %>%
      mutate( !!paste0(names(model_human)[i], "_model_prediction_threshold"):=case_when(
        !!(as.name(paste0(names(model_human)[i],"_Human")))>=(coords(roc.human, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Human",
        !!(as.name(paste0(names(model_human)[i],"_Livestock")))>=(coords(roc.livestock, "best", ret = "threshold", best.method="youden")[1,1]) ~ "Livestock",
        TRUE ~ 'NA')) 
      
    # Add prediction columns and create true/false column about prediction accuracy
    predictions_df <- predictions_df %>%
      full_join(rownames_to_column(predict_prob, var="Assembly")) %>%
      mutate( !!paste0(names(model_human)[i], "_model_accuracy_highest") := ifelse(!!(as.name(paste0(names(model_human)[i],"_prediction")))==Source.Type, "TRUE", "FALSE"),
              !!paste0(names(model_human)[i], "_model_accuracy_threshold") := ifelse(!!(as.name(paste0(names(model_human)[i],"_model_prediction_threshold")))==Source.Type, "TRUE", "FALSE")
      )
    # VESA stop deleting here
  i=i+1
}
write.table(predictions_df, file ="./out/predict_all.tsv", sep="\t", row.names = FALSE, quote = FALSE) # VESA change to your desird output folder
write.table(thresholds_df, file ="./out/thresholds.tsv", sep="\t", row.names = FALSE, quote = FALSE) # VESA change to your desird output folder
write.table(f1_all, file ="./out/f1_scores.tsv", sep="\t", row.names = TRUE, quote = FALSE)# VESA change to your desird output folder

 # VESA delete everything below 
#### Run human seqs through bps models ######
predictions_human_to_bps <- amr_class_human %>% filter(Source.Host=="Human") %>% select(Source.Host) %>% rownames_to_column(var="Assembly")
i=1
for (model in model_human) {
  model_filename <- paste("./in/models/", names(model_bps)[i], "_model.rds", sep="")
    #load(model_filename)
    rf_random <- readRDS(model_filename) 
    
    predict_prob <- predict.train(rf_random, newdata = model, type="prob") 
    
    # Get highest scored prediction
    predict_prob$prediction <- names(predict_prob)[1:3][apply(predict_prob[,1:3], 1, which.max)]
    # Change column names to include model name
    names(predict_prob) <- paste0(names(model_bps)[i], "_", names(predict_prob))
    # Add prediction columns and create true/false column about prediction accuracy
    predictions_human_to_bps <- predictions_human_to_bps %>%
      full_join(rownames_to_column(predict_prob, var="Assembly") )
  i=i+1
}
write.table(predictions_human_to_bps, file ="./out/predict_human_seqs_to_bps_models.tsv", sep="\t", row.names = FALSE, quote = FALSE)

#### Run bps seqs through human models ######
predictions_bps_to_human <- amr_class_bps %>% select(Source.Host) %>% rownames_to_column(var="Assembly")
i=1
for (model in model_bps) {
  model_filename <- paste("./in/models/", names(model_human)[i], "_model.rds", sep="")
  #load(model_filename)
  rf_random <- readRDS(model_filename) 
  
  predict_prob <- predict.train(rf_random, newdata = model, type="prob") 
  
  # Get highest scored prediction
  predict_prob$prediction <- names(predict_prob)[1:3][apply(predict_prob[,1:2], 1, which.max)]
  # Change column names to include model name
  names(predict_prob) <- paste0(names(model_human)[i], "_", names(predict_prob))
  # Add prediction columns and create true/false column about prediction accuracy
  predictions_bps_to_human <- predictions_bps_to_human %>%
    full_join(rownames_to_column(predict_prob, var="Assembly") )
  i=i+1
}
write.table(predictions_bps_to_human, file ="./out/predict_bps_seqs_to_human_models.tsv", sep="\t", row.names = FALSE, quote = FALSE)

