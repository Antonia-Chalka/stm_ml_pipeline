library(castor)
library(phytools)
library(tidyverse)
library(plotly)
library(ggtree)
library(ggnewscale)

# #Load tree
tree <- treeio::read.newick("./in/gubbins.node_labelled.final_tree.tre")
tree <- drop.tip(tree, c("Reference"))

# Metadata and Predicted Seqs
metadata_all <- read.csv("./in/all_metadata.csv", header=TRUE) 
metadata_all$Filename <- substr(metadata_all$Filename,1,nchar(metadata_all$Filename)-6)

phylogeny_predictions <- data.frame()

for(i in 1:length(tree$tip.label)){
  results <- find_nearest_tips(tree=tree, target_tips = i)
  index <- which.min(results$nearest_distance_per_tip[-i])
  
  if (index >= i) {
    index <- index + 1
  }
  phylogeny_prediction <- data.frame(Assembly = tree$tip.label[i],
                  Closest.Assembly = tree$tip.label[index],
                  Phylo.Distance = results$nearest_distance_per_tip[index],
                  Source.Host = (metadata_all %>% filter(Filename==tree$tip.label[i]) %>% select(Source.Host))[1,],
                  Closest.Assembly.Host = (metadata_all %>% filter(Filename==tree$tip.label[index]) %>% select(Source.Host))[1,]
                  ) %>%
    mutate( Phylogeny.Accuracy = ifelse( Closest.Assembly.Host == Source.Host, "TRUE", "FALSE"))
  
  phylogeny_predictions <- rbind(phylogeny_predictions, phylogeny_prediction)
}

write.table(phylogeny_predictions, file ="./out/phylogeny_predictions.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# Accuracy % 
phylogeny_predictions$Phylogeny.Accuracy <- as.factor(phylogeny_predictions$Phylogeny.Accuracy)
accuracy <- ((plyr::count(phylogeny_predictions$Phylogeny.Accuracy))[2,2]/length(phylogeny_predictions$Phylogeny.Accuracy))*100
accuracy

# accuracy per host
accuracy_host <- phylogeny_predictions %>% group_by(Source.Host) %>% count(Phylogeny.Accuracy) %>% 
  pivot_wider(names_from = Phylogeny.Accuracy, values_from=n,names_prefix = "Prediction_") %>%
  mutate(accuracy_per = Prediction_TRUE/(Prediction_FALSE + Prediction_TRUE)*100)
ggplot(accuracy_host, aes(x=Source.Host,y=accuracy_per)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(title="Phylogeny-based Prediction Accuracy By Host")
ggsave(filename = "./out/phylogeny_prediction_host.png",device = "png",dpi = 300)

# Distance vs Accuracy
p <- ggplot(phylogeny_predictions, aes(x=Phylogeny.Accuracy,y=Phylo.Distance)) +
  geom_boxplot() +
  theme_bw() +
  labs(title="Accuracy judged by closest neighbour",x="Prediction Outcome",y="Phylogenetic (Patristic) Distance ")
ggplotly(p)
htmlwidgets::saveWidget(as_widget(ggplotly(p)), "./out/phylogeny_prediction_accuracy.html")



