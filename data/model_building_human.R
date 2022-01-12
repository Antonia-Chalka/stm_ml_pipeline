library(tidyverse)
library(randomForest)
library(caret)
#library(doParallel)
#library(pROC)

set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=6) {
  length(args)
  stop("Exactly 6 arguments must be supplied (input file).n", call.=FALSE)
} 
##################################################### Load Input Data #####################################################
### AMR ###
# Class
amr_class_nonclonal <- read.table(args[1], header=TRUE, row.names = 1)
amr_class_nonclonal$Source.Host <- as.factor(amr_class_nonclonal$Source.Host)
# Gene
amr_gene_nonclonal <- read.table(args[2], header=TRUE, row.names = 1)
amr_gene_nonclonal$Source.Host <- as.factor(amr_gene_nonclonal$Source.Host)

### PV ###
pv_coll_nonclonal <- read.table(args[3], header=TRUE, sep="\t", row.names = 1)
pv_coll_nonclonal$Source.Host <- as.factor(pv_coll_nonclonal$Source.Host)

### IGR ###
igr_coll_nonclonal <- read.table(args[4], header = TRUE, sep="\t", row.names = 1)
igr_coll_nonclonal$Source.Host <- as.factor(igr_coll_nonclonal$Source.Host)

### SNP ###
snp_abudance_nonclonal <- read.table(args[5], header = TRUE, sep="\t", row.names = 1)
snp_abudance_nonclonal$Source.Host <- as.factor(snp_abudance_nonclonal$Source.Host)
# Set SNP factor levels to AGCT for model reliability 
snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"] <- lapply(snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"],factor,levels=c("C","G","T","A"))

##################################################### Make Models ##################################################### 
input_all <- list(amr_class_nonclonal= amr_class_nonclonal,
                  amr_gene_nonclonal = amr_gene_nonclonal,
                  pv_coll_nonclonal = pv_coll_nonclonal,
                  igr_coll_nonclonal = igr_coll_nonclonal,
                  snp_abudance_nonclonal = snp_abudance_nonclonal
)
# Initialise model performance stats dataframes #####
prediction_overall_all <- data.frame()
prediction_class_all <- data.frame()

###### All Host Models ####
i = 1
for (x in input_all) {
  #Initialise cluster
  # Set up parallel computing ####
  #cl <- makePSOCKcluster(args[6])
  #registerDoParallel(cl)
  
  # Hyperparameter testing
  control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", p=0.75, savePredictions="final", classProbs = TRUE)
  tunegrid <- expand.grid(.mtry = c(sqrt(ncol(x))))
  rf_random <- train(Source.Host~., data=x, method="rf", metric="Kappa", tuneGrid = tunegrid, trControl=control, importance=TRUE,ntree=500)
  
  #### Save model ####
  model_filename <- paste(names(input_all)[i], "_model_human.rds", sep="")
  saveRDS(rf_random,model_filename)
  
  ## Get Prediction scores ####
  prediction <- predict.train(rf_random, newdata = x)
  prediction_overall <- as.data.frame(as.matrix(confusionMatrix(prediction, x$Source.Host), what = "overall")) %>%
    rownames_to_column()
  prediction_overall$model=names(input_all)[i]
  prediction_overall_all <- rbind(prediction_overall_all,prediction_overall)
  prediction_class <- as.data.frame(as.matrix(confusionMatrix(prediction, x$Source.Host), what = "classes")) %>%
    rownames_to_column() %>%
    pivot_longer(cols = !rowname, names_to = "host", values_to="score")
  prediction_class$model=names(input_all)[i]
  prediction_class_all <- rbind(prediction_class_all, prediction_class)
  
  i = i + 1
  #stopCluster(cl)
}

# save prediction scores here ####
write.table(prediction_overall_all, file="prediction_overall_all_human.csv", quote=FALSE, sep=",")
write.table(prediction_class_all, file="prediction_class_all_human.csv", quote=FALSE, sep=",")

### Model scores plotting ####
accuracy_overall <- prediction_overall_all %>% filter(rowname %in% c("Accuracy","Kappa")) %>% ggplot() +
  geom_bar(aes(x=model, y=V1, fill=rowname), stat="identity", position=position_dodge()) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Accuracy of model", subtitle = "Calculated based on training data (25%, 10-fold validation)", y="%", x="Model Name")
ggsave(accuracy_overall, file ="accuracy_overal_human.png", device=png())

accuracy_class <- prediction_class_all %>% filter(rowname=="F1") %>% ggplot() +
  geom_bar(aes(x=model, y=score), stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Accuracy of model within hosts", subtitle = "Calculated based on training data (25%, 10-fold validation)", y="F1 %", x="Model Name") +
  facet_wrap(. ~ host)
ggsave(accuracy_class, file ="accuracy_hosts_human.png", device=png())
