library(tidyverse)
library(ggtree)
library(phytools)
library(ggnewscale)
set.seed(100)

##### Load Data ################################################ 
iqtree <- treeio::read.newick("in/gubbins.node_labelled.final_tree.tre")
iqtree <- drop.tip(iqtree, c("Reference"))
iqtree.rooted <- midpoint.root(iqtree)

predictions <- read.table("./out/predict_all.tsv", sep="\t", quote = "", header = TRUE)
phylo_pred <- read.table("./out/phylogeny_predictions.tsv", sep="\t", quote = "", header = TRUE)
baps <- read.table("./out/fastbaps_clusters.tsv", sep="\t", quote = "", header = TRUE)

mlst <- read.table("./in/mlst_extra.tsv",sep="\t",quote="",header = TRUE) %>% select(FILE,ST)
mlst$FILE <- substr(mlst$FILE,1,nchar(mlst$FILE)-6)

metadata_all <- read.csv("./in/all_metadata.csv", header=TRUE) 
metadata_all$Filename <- substr(metadata_all$Filename,1,nchar(metadata_all$Filename)-6)
metadata_all_combined <- full_join(metadata_all,baps, by = c("Filename"="assembly"))
metadata_all_combined <- full_join(metadata_all_combined,mlst, by = c("Filename"="FILE"))

metadata_nonclonal <- read.csv("./in/good_nonclonal_metadata.csv", header=TRUE)
metadata_nonclonal_combined <- left_join(metadata_nonclonal,baps, by = c("Filename"="assembly"))
metadata_nonclonal_combined <- left_join(metadata_nonclonal_combined,mlst, by = c("Filename"="FILE"))


##### Dataset Counts Plots ################################################ 
ggplot(data= metadata_nonclonal_combined,aes(x=Source.Host, fill=ST)) +
  geom_bar() +
  labs(title="A. Nonclonal STm Dataset", y="Counts",x="Host") +
  theme_bw()
ggsave("out/nc_dataset_counts.png", device = "png",dpi = 300, width = 15,height=18, units="cm")
dev.off()

##### AUC vs Highest Score Accuracy Plots ################################################ 
predictions_comparison <- select(predictions, "Assembly", contains("accuracy") ) %>%
  pivot_longer(!Assembly, names_to="Feature", values_to="Accuracy" ) %>%
  count(Feature,Accuracy) %>%
  filter(!is.na(Accuracy)) %>%
  mutate(Accuracy_Type=str_extract(Feature, "threshold|highest"), Feature= str_remove(Feature, "_threshold|_highest")) %>%
  mutate(Model_Type=str_extract(Feature, "all|bps|human"), Feature= str_remove(Feature, "_all|_bps|_human")) %>%
  mutate(Feature= str_remove(Feature, "_model_accuracy")) 

ggplot(predictions_comparison, aes(fill=Accuracy,y=n, x=Accuracy_Type)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  labs("Comparison of prediction methods across all models", x="Method", y="Percentage")+
  facet_wrap(~ Feature + Model_Type)
ggsave("out/accuracy_methods.png", device = "png",dpi = 300,width = 25,height=25, units="cm")
dev.off()

##### Tree Plots ################################################ 
predictions_heatmap <- select(predictions, "Assembly", 
         PV.All = pv_all_model_accuracy_threshold, PV.BPS = pv_bps_model_accuracy_threshold, PV.Human = pv_human_model_accuracy_threshold, 
         IGR.All = igr_all_model_accuracy_threshold, IGR.BPS = igr_bps_model_accuracy_threshold, IGR.Human = igr_human_model_accuracy_threshold,
         SNP.All = snp_abudance_all_model_accuracy_threshold, SNP.BPS = snp_abudance_bps_model_accuracy_threshold, SNP.Human = snp_abudance_human_model_accuracy_threshold,
         AMR.Class.All = amr_class_all_model_accuracy_threshold, AMR.Class.BPS = amr_class_bps_model_accuracy_threshold, AMR.Class.Human = amr_class_human_model_accuracy_threshold,
         AMR.Gene.All = amr_gene_all_model_accuracy_threshold, AMR.Gene.BPS = amr_gene_bps_model_accuracy_threshold, AMR.Gene.Human = amr_gene_human_model_accuracy_threshold
  ) %>%
  column_to_rownames(var="Assembly")

year_heatmap <- metadata_all_combined %>% 
  select(Filename, Collection.Year) %>% 
  column_to_rownames(var="Filename") %>%
  mutate(Collection.Year = as.factor(Collection.Year)) 

mlst_heatmap <- metadata_all_combined %>% 
  select(Filename, ST) %>% 
  column_to_rownames(var="Filename") %>%
  mutate(ST = as.factor(ST)) 

host_heatmap <- metadata_all_combined %>% 
  select(Filename, Host = Source.Host) %>% 
  column_to_rownames(var="Filename") %>%
  mutate(Host = as.factor(Host))

baps_heatmap <- metadata_all_combined %>%
  select(Filename, Cluster=clusters) %>%
  column_to_rownames(var="Filename") %>%
  mutate(Cluster= as.factor(Cluster)) 

phylogeny_heatmap <- phylo_pred %>% 
  select(Assembly,Phylo.Prediction =Phylogeny.Accuracy) %>%
  column_to_rownames(var="Assembly") 

### Tree plot 
tree_plot <- ggtree(iqtree.rooted) +
  geom_treescale(x=1, y=4700, width=100)
tree_plot
baps_treeplot <- gheatmap(tree_plot, baps_heatmap , color=NULL, offset = 100, width = 0.08, font.size=4, colnames_position = "top", colnames_angle=60, colnames_offset_y = -0.2, hjust=0)+
  scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7"), name="Baps Cluster") + 
  coord_cartesian(clip = 'off') +
  ylim(NA, 5500)
baps_treeplot

host_treeplot <- baps_treeplot + new_scale_fill()
host_treeplot <-  gheatmap(host_treeplot, host_heatmap , color=NULL, offset = 200, width = 0.08, font.size=4, colnames_position = "top", colnames_angle=60, colnames_offset_y = -0.2, hjust=0) +
  scale_fill_manual( name="Host", values=c("darkgoldenrod1", "coral2","cadetblue1", "mediumorchid1")) +
  coord_cartesian(clip = 'off') +
  ylim(NA, 5500) 
host_treeplot

phylo_treeplot <- host_treeplot + new_scale_fill()
phylo_treeplot <- gheatmap(phylo_treeplot, phylogeny_heatmap , color=NULL, offset = 300, width = 0.08, font.size=4, colnames_position = "top", colnames_angle=60, colnames_offset_y = -0.2, hjust=0) +
  scale_fill_manual(values = c("#DA5724", "#74D944"), na.value="white", name="Phylo Predictions") +
  ylim(NA, 5500)
phylo_treeplot

predictions_treeplot <- phylo_treeplot + new_scale_fill()

predictions_treeplot <-  gheatmap(predictions_treeplot, predictions_heatmap , color=NULL, offset =400, width = 1, font.size=4, colnames_position = "top", colnames_angle=60, colnames_offset_y = -0.2, hjust=0) +
  scale_fill_manual(values = c("#DA5724", "#74D944"), na.value="white", name="Predictions") + 
  coord_cartesian(clip = 'off') +
  ylim(NA, 5300)
predictions_treeplot

ggsave(predictions_treeplot, device=png(), file="./out/prediction_treeplot.png", dpi=400, height = 12, width=10)
dev.off()

##### Stacked Predictions Plots ####################################################
thresholds_df <- read.table(file ="./out/thresholds.tsv", header=TRUE,sep="\t", quote = "") 
dir.create("./out/stacked_predictions")

stacked_predictions <- predictions %>% 
  select(Assembly, Source.Host, 
         igr_all_prediction,igr_bps_prediction,igr_human_prediction, 
         igr_all_Bovine,igr_all_Human,igr_all_Poultry,igr_all_Swine,
         igr_bps_Bovine,igr_bps_Poultry,igr_bps_Swine,
         igr_human_Human,igr_human_Livestock) %>%
  pivot_longer(cols = contains("prediction"), names_to = "model", values_to = "final_prediction")%>%
  pivot_longer(cols=matches("Bovine|Poultry|Human|Swine|Livestock"),names_to="host",values_to="igr_score")%>%
  filter((model=="igr_all_prediction" & str_detect(host,"all")) | 
         (model=="igr_bps_prediction" & str_detect(host,"bps")) | 
         (model=="igr_human_prediction" & str_detect(host,"human"))) %>%
  mutate(model = str_remove_all(model, "igr_|_prediction"),
         host = str_remove_all(host,"igr_all_|igr_bps_|igr_human_")) %>%
  filter(!(Source.Host=="Human" & model=="bps"))
write.table(stacked_predictions, file="out/stacked_predictions/stacked_predictions.tsv",row.names = FALSE,quote = FALSE,col.names = TRUE,sep="\t")

# All Bovine
stacked_predictions_bovine_all <- stacked_predictions %>% filter(Source.Host=="Bovine" & model=="all") %>%
  mutate(host=factor(host, levels=c("Human","Poultry","Swine","Bovine")))
stacked_predictions_bovine_all_plot <- ggplot(data=stacked_predictions_bovine_all, 
                                              aes(x=factor(Assembly,levels=Assembly[host=="Bovine"][order(igr_score[host=="Bovine"])]), 
                                                  y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Bovine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("coral2","cadetblue1", "mediumorchid1","darkgoldenrod1"))
ggsave(stacked_predictions_bovine_all_plot,device=png(), file="./out/stacked_predictions/all_bovine.png", dpi=300, height = 3, width=5)

#BPS Bovine
stacked_predictions_bovine_bps <- stacked_predictions %>% filter(Source.Host=="Bovine" & model=="bps") %>%
  mutate(host=factor(host, levels=c("Poultry","Swine","Bovine")))
stacked_predictions_bovine_bps_plot <- ggplot(data=stacked_predictions_bovine_bps, 
                                              aes(x=factor(Assembly,levels=Assembly[host=="Bovine"][order(igr_score[host=="Bovine"])]), 
                                                  y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_bps"& Host=="Bovine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("cadetblue1", "mediumorchid1","darkgoldenrod1"))
ggsave(stacked_predictions_bovine_bps_plot,device=png(), file="./out/stacked_predictions/bps_bovine.png", dpi=300, height = 3, width=5)

# All Poultry
stacked_predictions_poultry_all <- stacked_predictions %>% filter(Source.Host=="Poultry" & model=="all") %>%
  mutate(host=factor(host, levels=c("Human","Swine","Bovine","Poultry")))
stacked_predictions_poultry_all_plot <- ggplot(data=stacked_predictions_poultry_all, 
                                               aes(x=factor(Assembly,levels=Assembly[host=="Poultry"][order(igr_score[host=="Poultry"])]), 
                                                   y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Poultry")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("coral2", "mediumorchid1","darkgoldenrod1","cadetblue1"))
ggsave(stacked_predictions_poultry_all_plot,device=png(), file="./out/stacked_predictions/all_poultry.png", dpi=300, height = 3, width=5)

# BPS Poultry
stacked_predictions_poultry_bps <- stacked_predictions %>% filter(Source.Host=="Poultry" & model=="bps") %>%
  mutate(host=factor(host, levels=c("Swine","Bovine","Poultry")))
stacked_predictions_poultry_bps_plot <- ggplot(data=stacked_predictions_poultry_bps, 
                                               aes(x=factor(Assembly,levels=Assembly[host=="Poultry"][order(igr_score[host=="Poultry"])]), 
                                                   y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_bps"& Host=="Poultry")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c( "mediumorchid1","darkgoldenrod1","cadetblue1"))
ggsave(stacked_predictions_poultry_bps_plot,device=png(), file="./out/stacked_predictions/bps_poultry.png", dpi=300, height = 3, width=5)

# All Swine
stacked_predictions_swine_all <- stacked_predictions %>% filter(Source.Host=="Swine" & model=="all") %>%
  mutate(host=factor(host, levels=c("Human","Poultry","Bovine","Swine")))
stacked_predictions_swine_all_plot <- ggplot(data=stacked_predictions_swine_all, 
                                             aes(x=factor(Assembly,levels=Assembly[host=="Swine"][order(igr_score[host=="Swine"])]), 
                                                 y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Swine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("coral2","cadetblue1","darkgoldenrod1", "mediumorchid1"))
ggsave(stacked_predictions_swine_all_plot,device=png(), file="./out/stacked_predictions/all_swine.png", dpi=300, height = 3, width=5)

# BPS Swine
stacked_predictions_swine_bps <- stacked_predictions %>% filter(Source.Host=="Swine" & model=="bps") %>%
  mutate(host=factor(host, levels=c("Poultry","Bovine","Swine")))
stacked_predictions_swine_bps_plot <- ggplot(data=stacked_predictions_swine_bps, 
                                              aes(x=factor(Assembly,levels=Assembly[host=="Swine"][order(igr_score[host=="Swine"])]), 
                                                  y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_bps"& Host=="Swine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("cadetblue1","darkgoldenrod1", "mediumorchid1"))
ggsave(stacked_predictions_swine_bps_plot,device=png(), file="./out/stacked_predictions/bps_swine.png", dpi=300, height = 3, width=5)

# All Human
stacked_predictions_human_all <- stacked_predictions %>% filter(Source.Host=="Human" & model=="all") %>%
  mutate(host=factor(host, levels=c("Poultry","Swine","Bovine","Human")))
stacked_predictions_human_all_plot <- ggplot(data=stacked_predictions_human_all, 
                                             aes(x=factor(Assembly,levels=Assembly[host=="Human"][order(igr_score[host=="Human"])]), 
                                                 y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Human")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("cadetblue1", "mediumorchid1","darkgoldenrod1","coral2"))
ggsave(stacked_predictions_human_all_plot,device=png(), file="./out/stacked_predictions/all_human.png", dpi=300, height = 3, width=5)

# Human Human
stacked_predictions_human_human <- stacked_predictions %>% filter(Source.Host=="Human" & model=="human") %>%
  mutate(host=factor(host, levels=c("Livestock","Human")))
stacked_predictions_human_human_plot <- ggplot(data=stacked_predictions_human_human, 
                                             aes(x=factor(Assembly,levels=Assembly[host=="Human"][order(igr_score[host=="Human"])]), 
                                                 y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_human"& Host=="Human")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("darkolivegreen1","coral2"))
ggsave(stacked_predictions_human_human_plot,device=png(), file="./out/stacked_predictions/human_human.png", dpi=300, height = 3, width=5)

# Human Livestock
stacked_predictions_human_livestock<- stacked_predictions %>% filter(Source.Host!="Human" & model=="human") %>%
  mutate(host=factor(host, levels=c("Human","Livestock")))
stacked_predictions_human_livestock_plot <- ggplot(data=stacked_predictions_human_livestock, 
                                               aes(x=factor(Assembly,levels=Assembly[host=="Livestock"][order(igr_score[host=="Livestock"])]), 
                                                   y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1 - filter(thresholds_df, model=="igr_human"& Host=="Human")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("coral2","darkolivegreen1"))
ggsave(stacked_predictions_human_livestock_plot,device=png(), file="./out/stacked_predictions/human_livestock.png", dpi=300, height = 3, width=5)


##### Model Accuracy Plots ################################################ 
accuracy_all <- list.files(path = "./in/predictions/", pattern = "*prediction_overall.csv", full.names = TRUE) %>% 
  lapply(read.table, quote="", sep=",", header=TRUE) %>%
  bind_rows %>%
  mutate(type=str_extract(model, "all|bps|human"), model= str_remove(model, "_all|_bps|_human")) %>% 
  filter(rowname %in% c("Kappa")) %>%
  dplyr::rename(Features=model)

accuracy_all$Features[accuracy_all$Features=="amr_class"] <- "AMR Class"
accuracy_all$Features[accuracy_all$Features=="amr_gene"] <- "AMR Gene"
accuracy_all$Features[accuracy_all$Features=="igr"] <- "IGR"
accuracy_all$Features[accuracy_all$Features=="pv"] <- "PV"
accuracy_all$Features[accuracy_all$Features=="snp_abudance"] <- "SNP"

accuracy_overall_plot <- ggplot(accuracy_all,aes(x=type, y=V1, fill=Features)) +
  geom_bar( stat="identity", position= position_dodge()) +
  theme_bw()  + 
  labs(title = "Accuracy (Kappa)", subtitle = "Calculated based on training data (25%, 10-fold validation)", y="%", x="Model Type")
  
accuracy_overall_plot
ggsave(accuracy_overall_plot, file ="./out/accuracy_overall.png", device=png(), height=15, width=30, units="cm")
dev.off()

accuracy_class <- list.files(path = "./in/predictions/", pattern = "*prediction_class.csv", full.names = TRUE) %>% 
  lapply(read.table, quote="", sep=",", header=TRUE) %>%
  bind_rows %>%
  mutate(type=str_extract(model, "all|bps|human")) %>% 
  filter(rowname=="F1" & !host=="V1") 

accuracy_class_plot <- ggplot(accuracy_class, aes(x=model, y=score)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Accuracy of model within hosts", subtitle = "Calculated based on training data (25%, 10-fold validation)", y="F1 %", x="Model Name") +
  facet_wrap(. ~ host)

accuracy_class_plot
ggsave(accuracy_class_plot, file ="./out/accuracy_hosts.png", device=png())
dev.off()

#### ST313 Plots ############################################
predictions_313_all<- read.table(file = "./in/prediction_all.tsv", sep="\t", header = TRUE) %>%
  select(!prediction) %>%
  pivot_longer(cols = c("Human", "Bovine", "Swine","Poultry"), names_to = "Host", values_to = "Confidence") %>%
  mutate(Assembly = str_replace(Assembly, ".fasta", ""))

predictions_313_all_plot <- ggplot(data=predictions_313_all, aes(x=Assembly, y=Confidence, fill=Host)) +
  geom_bar(stat="identity") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~model, nrow=1)+
  labs(title="A. Confidence Scores of All-Host Models on ST313 Sequences", x= "Sequence", y="Confidence Score")
predictions_313_all_plot
ggsave(predictions_313_all_plot, file ="./out/predictions_313_all_plot.png", device=png(), height=30, width=30, units="cm")
dev.off()

predictions_313_bps <- read.table(file = "./in/prediction_bps.tsv", sep="\t", header = TRUE) %>%
  select(!prediction) %>%
  pivot_longer(cols = c("Bovine", "Swine","Poultry"), names_to = "Host", values_to = "Confidence") %>%
  mutate(Assembly = str_replace(Assembly, ".fasta", ""))
predictions_313_bps_plot <- ggplot(data=predictions_313_bps, aes(x=Assembly, y=Confidence, fill=Host)) +
  geom_bar(stat="identity") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~model, nrow=1)+
  labs(title="C. Confidence Scores of BPS Models on ST313 Sequences", x= "Sequence", y="Confidence Score")
predictions_313_bps_plot
ggsave(predictions_313_bps_plot, file ="./out/predictions_313_bps_plot.png", device=png(), height=30, width=30, units="cm")
dev.off()

predictions_313_human <- read.table(file = "./in/prediction_human.tsv", sep="\t", header = TRUE) %>%
  select(!prediction) %>%
  pivot_longer(cols = c("Human", "Livestock"), names_to = "Host", values_to = "Confidence") %>%
  mutate(Assembly = str_replace(Assembly, ".fasta", ""))
predictions_313_human_plot <- ggplot(data=predictions_313_human, aes(x=Assembly, y=Confidence, fill=Host)) +
  geom_bar(stat="identity") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~model, nrow=1)+
  labs(title="B. Confidence Scores of Human-Scoring Models on ST313 Sequences", x= "Sequence", y="Confidence Score")
predictions_313_human_plot
ggsave(predictions_313_human_plot, file ="./out/predictions_313_human_plot.png", device=png(), height=30, width=30, units="cm")
dev.off()

#### Increase in Accuracy From Phylo Models Plots ############################################
wrong_predictions <- phylo_pred %>% filter(Phylogeny.Accuracy==FALSE)

wrong_predictions_all <- left_join(wrong_predictions, predictions) %>% 
  select("Assembly", contains("accuracy"), -Phylogeny.Accuracy ) %>%
  pivot_longer(cols=contains("accuracy"), names_to = "Model", values_to = "Accuracy") %>%
  count(Model, Accuracy) %>%
  filter(!is.na(Accuracy)) %>%
  pivot_wider(names_from = "Accuracy",names_prefix = "Accuracy_", values_from = "n") %>%
  mutate(Total= Accuracy_FALSE + Accuracy_TRUE, Perc_Increase=(Accuracy_TRUE/Total)*100) %>%
  select(Perc_Increase, Model)

wrong_predictions_all_plot <- ggplot(wrong_predictions_all,aes(x=Perc_Increase, y=Model)) +
  geom_bar( stat="identity", position= position_dodge(), fill = "#FF6666") +
  theme_bw()  + 
  labs(title = "Increase in Accuracy Compared to Phylogenetic Host Prediction", x="Percentage Increase",)

wrong_predictions_all_plot
ggsave(wrong_predictions_all_plot, file ="./out/wrong_predictions_all.png", device=png(), height=15, width=15, units="cm")
dev.off()

#### Poultry Cluster Seqs on BPS Models ############################################
poultry_cluster <- metadata_all_combined %>% filter(clusters==7 & Source.Host=="Human")

predictions_humanbps <- read.table("./out/predict_human_seqs_to_bps_models.tsv", sep="\t", header=TRUE , quote = "")
poultry_cluster <- left_join(poultry_cluster, predictions_humanbps, by=c("Filename"="Assembly")) %>% 
  select("Filename",  contains("bps")) 

#### BPS Seqs on Human ############################################
predict_bps_seqs_to_human_models <- read.table(file = "./out/predict_bps_seqs_to_human_models.tsv", sep="\t", header = TRUE) %>% 
  select(Source.Host, contains("human_Human")) %>%
  pivot_longer(!Source.Host, names_to="Model", values_to="Score" ) %>%
  mutate(Model=str_remove(Model,"_human_Human"))

predict_bps_seqs_to_human_models$Model[predict_bps_seqs_to_human_models$Model=="amr_class"] <- "AMR Class"
predict_bps_seqs_to_human_models$Model[predict_bps_seqs_to_human_models$Model=="amr_gene"] <- "AMR Gene"
predict_bps_seqs_to_human_models$Model[predict_bps_seqs_to_human_models$Model=="igr"] <- "IGR"
predict_bps_seqs_to_human_models$Model[predict_bps_seqs_to_human_models$Model=="pv"] <- "PV"
predict_bps_seqs_to_human_models$Model[predict_bps_seqs_to_human_models$Model=="snp_abudance"] <- "SNP"

predict_bps_seqs_to_human_models_plot <- ggplot(predict_bps_seqs_to_human_models, aes(x=Model, y=Score, fill=Source.Host)) + 
  geom_boxplot() +
  theme_bw() +
  labs( fill="Host", y="Confidence Score") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_human")$optimal_threshold)

predict_bps_seqs_to_human_models_plot
ggsave(predict_bps_seqs_to_human_models_plot, file ="./out/predict_bps_seqs_to_human_models_plot.png", device=png(), height=13, width=18, units="cm")
dev.off()

predict_bps_seqs_to_human_models_context <- read.table(file = "./out/predict_bps_seqs_to_human_models.tsv", sep="\t", header = TRUE) %>% 
  select(Source.Host, Assembly, contains("human_Human")) %>%
  left_join( read.csv(file="./in/metadata_quast_all_merged_filtered_discarded_modified.csv") %>%
               select( Source.Context, Filename) %>%
               dplyr::rename(Assembly=Filename)) %>%
  select(Source.Host, Source.Context,contains("human_Human")) %>%
  pivot_longer(cols = ,contains("human_Human"), names_to="Model", values_to="Score" ) %>%
  mutate(Model=str_remove(Model,"_human_Human")) %>%
  ggplot(aes(x=Model, y=Score, fill=Source.Host)) + 
  geom_boxplot() +
  theme_bw() +
  labs(title="Livestock Sequences on Human Scoring Models", fill="Host") +
  facet_wrap(~Source.Context)
ggsave(predict_bps_seqs_to_human_models_context, file ="./out/predict_bps_seqs_to_human_models_context_plot.png", device=png(), height=30, width=30, units="cm")
dev.off()

predict_bps_seqs_to_human_models %>% filter(Score>0.5) %>% count(Source.Host, Model)

#### Human Seqs on BPS ############################################
predict_human_seqs_to_bps_models <- read.table(file = "./out/predict_human_seqs_to_bps_models.tsv", sep="\t", header = TRUE)%>% 
  filter(Source.Host=="Human")%>%
  select( contains("prediction")) %>%
  pivot_longer(everything(), names_to="Model", values_to="Prediction" ) %>%
  mutate(Model=str_remove(Model,"_bps_prediction")) %>%
  count(Prediction, Model)

predict_human_seqs_to_bps_models_plot <- ggplot(predict_human_seqs_to_bps_models, aes(x=Model, y=n, fill=Prediction)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  coord_flip()+
  labs(title="Human Sequences on BPS Models", fill="Predicted Host")
predict_human_seqs_to_bps_models_plot
ggsave(predict_human_seqs_to_bps_models_plot, file ="./out/predict_human_seqs_to_bps_models.png", device=png(), height=30, width=30, units="cm")
dev.off()

###### Stacked Predictions Plots ####################################################
thresholds_df <- read.table(file ="./out/thresholds.tsv", header=TRUE,sep="\t", quote = "") 
dir.create("./out/stacked_predictions")

stacked_predictions_human_to_bps <-  read.table(file = "./out/predict_human_seqs_to_bps_models.tsv", sep="\t", header = TRUE)%>% 
  filter(Source.Host=="Human")%>% 
  select(Assembly, Source.Host, 
         igr_bps_prediction,
         igr_bps_Bovine,igr_bps_Poultry,igr_bps_Swine)%>%
  pivot_longer(cols = contains("prediction"), names_to = "model", values_to = "final_prediction")%>%
  pivot_longer(cols=matches("Bovine|Poultry|Human|Swine|Livestock"),names_to="host",values_to="igr_score")%>%
  filter( model=="igr_bps_prediction" & str_detect(host,"bps"))%>%
  mutate(model = str_remove_all(model, "igr_|_prediction"),
         host = str_remove_all(host,"igr_bps_"))

write.table(stacked_predictions_human_to_bps, file="out/stacked_predictions/stacked_predictions_human_to_bps.tsv",row.names = FALSE,quote = FALSE,col.names = TRUE,sep="\t")


stacked_predictions_human_to_bps_bovine_plot <- ggplot(data=stacked_predictions_human_to_bps %>% filter(final_prediction=="Bovine") %>%
                                                         mutate(host=factor(host, levels=c("Poultry","Swine","Bovine"))), 
                                              aes(x=factor(Assembly,levels=Assembly[host=="Bovine"][order(igr_score[host=="Bovine"])]), 
                                                  y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Bovine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("cadetblue1", "mediumorchid1","darkgoldenrod1"))
ggsave(stacked_predictions_human_to_bps_bovine_plot,device=png(), file="./out/stacked_predictions/human_to_bps_bovine.png", dpi=300, height = 3, width=5)


stacked_predictions_human_to_bps_swine_plot <- ggplot(data=stacked_predictions_human_to_bps %>% filter(final_prediction=="Swine")%>%
                                                        mutate(host=factor(host, levels=c("Bovine","Poultry","Swine"))), 
                                                        aes(x=factor(Assembly,levels=Assembly[host=="Swine"][order(igr_score[host=="Swine"])]), 
                                                            y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Swine")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("darkgoldenrod1","cadetblue1", "mediumorchid1"))
stacked_predictions_human_to_bps_swine_plot
ggsave(stacked_predictions_human_to_bps_swine_plot,device=png(), file="./out/stacked_predictions/human_to_bps_swine.png", dpi=300, height = 3, width=5)

stacked_predictions_human_to_bps_poultry_plot <- ggplot(data=stacked_predictions_human_to_bps %>% filter(final_prediction=="Poultry")%>%
                                                        mutate(host=factor(host, levels=c("Bovine","Swine","Poultry"))), 
                                                      aes(x=factor(Assembly,levels=Assembly[host=="Poultry"][order(igr_score[host=="Poultry"])]), 
                                                          y=igr_score, fill=host, width=1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = filter(thresholds_df, model=="igr_all"& Host=="Poultry")$optimal_threshold) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Confidence Score") +
  scale_fill_manual( name="Host", values=c("darkgoldenrod1", "mediumorchid1","cadetblue1"))


stacked_predictions_human_to_bps_poultry_plot
ggsave(stacked_predictions_human_to_bps_poultry_plot,device=png(), file="./out/stacked_predictions/human_to_bps_poultry.png", dpi=300, height = 3, width=5)

# stacked_predictions_human_to_bps %>% group_by(final_prediction,host)%>% dplyr::summarise(Mean=mean(igr_score))


###### Human to BPS chronology ####

totals_year <- left_join(stacked_predictions_human_to_bps,metadata_all_combined,by=c("Assembly"="Filename")) %>% 
  count(Collection.Year) %>% 
  mutate(totals_year=n/3) %>% 
  select(!n)
new_human <-left_join(stacked_predictions_human_to_bps,metadata_all_combined,by=c("Assembly"="Filename")) %>% 
  filter(str_detect(Assembly, "^A")) %>% 
  group_by(final_prediction) %>% 
  count(Collection.Year) %>% 
  mutate(n=n/3) %>%
  rename(count_new=n)

old_human <-left_join(stacked_predictions_human_to_bps,metadata_all_combined,by=c("Assembly"="Filename")) %>% 
  filter(str_detect(Assembly, "^S")) %>% 
  group_by(final_prediction) %>% 
  count(Collection.Year) %>% 
  mutate(n=n/3) %>%
  rename(count_old=n) 


bps_to_human_time <-left_join(stacked_predictions_human_to_bps,metadata_all_combined,by=c("Assembly"="Filename")) %>% 
  count(Collection.Year,final_prediction) %>% 
  mutate(n=n/3) %>%
  full_join(totals_year, by="Collection.Year") %>% 
  mutate(perc_all= (n/totals_year)*100) %>%
  rename(count_all=n) %>%
  full_join(new_human, by=c("Collection.Year","final_prediction")) %>%
  mutate(perc_new= (count_new/totals_year)*100) %>%
  full_join(old_human, by=c("Collection.Year","final_prediction")) %>%
  mutate(perc_old= (count_old/totals_year)*100,
         final_prediction=factor(final_prediction, levels=c("Bovine","Poultry","Swine"))) 
  
write.table(bps_to_human_time, file="out/stacked_predictions/bps_to_human_time.tsv",row.names = FALSE,quote = FALSE,col.names = TRUE,sep="\t")


ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=perc_all, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  scale_fill_manual( name="BPS Model Prediction", values=c("darkgoldenrod1","cadetblue1", "mediumorchid1")) +
  theme_classic() +
  labs(x="Collection Year",y="% Counts", fill="IGR BPS Model Prediction") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_perc_all.png", dpi=300, height = 3, width=5)

ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=count_all, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  scale_fill_manual( name="BPS Prediction", values=c("darkgoldenrod1","cadetblue1", "mediumorchid1")) +
  theme_classic() +
  labs(x="Collection Year",y="Counts",fill="BPS Prediction") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_count_all.png", dpi=300, height = 3, width=5)

ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=count_new, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  scale_fill_manual( name="IGR BPS Model Prediction", values=c("darkgoldenrod1","cadetblue1", "mediumorchid1")) +
  theme_classic() +
  labs(x="Collection Year",y="Counts (New)") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_count_new.png", dpi=300, height = 3, width=5)

ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=perc_new, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  scale_fill_manual( name="IGR BPS Model Prediction", values=c("darkgoldenrod1","cadetblue1", "mediumorchid1")) +
  theme_classic() +
  labs(x="Collection Year",y="% Counts (New)") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_perc_new.png", dpi=300, height = 3, width=5)

ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=count_old, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  labs(x="Collection Year",y="Counts (Old)") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_count_old.png", dpi=300, height = 3, width=5)

ggplot(data=bps_to_human_time,
       aes(x=Collection.Year, 
           y=perc_old, fill=final_prediction)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  labs(x="Collection Year",y="% Counts (Old)") 
ggsave(device=png(), file="./out/stacked_predictions/bps_to_human_time_perc_old.png", dpi=300, height = 3, width=5)



##### Human seqs with 0 ybal cluster ####
human_ybal_0  <- pv_all %>% 
  filter(ybaL...fsr...kefB ==0, Source.Host=="Human") %>% 
  rownames_to_column(var="Assembly") %>%
  select("Assembly") %>%
  left_join(read.table(file = "./out/predict_human_seqs_to_bps_models.tsv", sep="\t", header = TRUE), by="Assembly") %>%
  select( contains("prediction")) %>%
  pivot_longer(everything(), names_to="Model", values_to="Prediction" ) %>%
  mutate(Model=str_remove(Model,"_bps_prediction")) %>%
  count(Prediction, Model)

human_group_3075_1  <- pv_all %>% 
  filter(group_3075 ==1, Source.Host=="Human") %>% 
  rownames_to_column(var="Assembly") %>%
  select("Assembly") %>%
  left_join(read.table(file = "./out/predict_human_seqs_to_bps_models.tsv", sep="\t", header = TRUE), by="Assembly") %>%
  select( contains("prediction")) %>%
  pivot_longer(everything(), names_to="Model", values_to="Prediction" ) %>%
  mutate(Model=str_remove(Model,"_bps_prediction")) %>%
  count(Prediction, Model)

#### Human Dataset Comparison ############################################

human_comparison <- filter(predictions, Source.Host =="Human") %>%
  select("Assembly", contains("accuracy_threshold") ) %>% 
  mutate(Dataset= if_else(str_starts(Assembly, "AA"), "New", "Old")) %>%
  select(!Assembly) %>%
  pivot_longer(!Dataset, names_to="Model", values_to="Accuracy" )%>%
  count(Dataset,Model,Accuracy) %>%
  filter(!is.na(Accuracy)) %>%
  group_by(Dataset, Model) %>%
  mutate(Model= str_remove(Model, "_accuracy_threshold"), freq = (n / sum(n))*100 ) 

#### Human PVs of interest ###############################################
# # PVs: ybal/kefb, 3075, 5859, 732
# 
# pv_heatmap <- pv_all %>% #pv_all from predict_all script
#   select(ybaL...fsr...kefB,group_3075, group_1923, group_4147, group_1077) %>%
#   mutate(ybaL...fsr...kefB = as.factor(ybaL...fsr...kefB),group_1923=as.factor(group_1923), group_4147=as.factor(group_4147), group_1077=as.factor(group_1077))
# 
# humanpv_treeplot <- phylo_treeplot + new_scale_fill()
# humanpv_treeplot <- gheatmap(humanpv_treeplot, pv_heatmap , color=NULL, offset = 400, width = 0.8, font.size=4, colnames_position = "top", colnames_angle=60, colnames_offset_y = -0.2, hjust=0) +
#   scale_fill_discrete( name="Presence")+
#   ylim(NA, 5500)
# humanpv_treeplot
# ggsave(humanpv_treeplot, device=png(), file="./out/humanpv_treeplot.png", dpi=400, height = 12, width=10)
# dev.off()
# 
# pv_cor <- pv_all %>% select(!Source.Host) %>%
#   as.matrix %>%
#   cor %>%
#   as.data.frame %>%
#   rownames_to_column(var = 'var1') %>%
#   gather(var2, value, -var1) %>%
#   filter(var1=="ybaL...fsr...kefB")
# 


