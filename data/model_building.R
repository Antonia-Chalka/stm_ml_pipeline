library(tidyverse)
library(randomForest)
library(caret)

set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=10) {
  length(args)
  stop("Exactly 10 arguments must be supplied (input file).n", call.=FALSE)
} 
##################################################### Load Input Data #####################################################
### AMR ###
# Class
amr_class_nonclonal <- read.table(args[1], header=TRUE, row.names = 1)
amr_class_nonclonal$Source.Host <- as.factor(amr_class_nonclonal$Source.Host)
amr_class_nonclonal_bps <- read.table(args[2], header=TRUE, row.names = 1)
amr_class_nonclonal_bps$Source.Host <- as.factor(amr_class_nonclonal_bps$Source.Host)
# Gene
amr_gene_nonclonal <- read.table(args[3], header=TRUE, row.names = 1)
amr_gene_nonclonal$Source.Host <- as.factor(amr_gene_nonclonal$Source.Host)
amr_gene_nonclonal_bps <- read.table(args[4], header=TRUE, row.names = 1)
amr_gene_nonclonal_bps$Source.Host <- as.factor(amr_gene_nonclonal_bps$Source.Host)

### PV ###
pv_coll_nonclonal <- read.table(args[5], header=TRUE, sep="\t", row.names = 1)
pv_coll_nonclonal$Source.Host <- as.factor(pv_coll_nonclonal$Source.Host)
pv_coll_nonclonal_bps <- read.table(args[6], header=TRUE, sep="\t", row.names = 1)
pv_coll_nonclonal_bps$Source.Host <- as.factor(pv_coll_nonclonal_bps$Source.Host)

### IGR ###
igr_coll_nonclonal <- read.table(args[7], header = TRUE, sep="\t", row.names = 1)
igr_coll_nonclonal$Source.Host <- as.factor(igr_coll_nonclonal$Source.Host)
igr_coll_monoclonal_bps <- read.table(args[8], header = TRUE, sep="\t", row.names = 1)
igr_coll_monoclonal_bps$Source.Host <- as.factor(igr_coll_monoclonal_bps$Source.Host)

### SNP ###
snp_abudance_nonclonal <- read.table(args[9], header = TRUE, sep="\t", row.names = 1)
snp_abudance_nonclonal$Source.Host <- as.factor(snp_abudance_nonclonal$Source.Host)
# Set SNP factor levels to AGCT for model reliability 
snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"] <- lapply(snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"],factor,levels=c("C","G","T","A"))
snp_abudance_nonclonal_bps <- read.table(args[10], header = TRUE, sep="\t", row.names = 1)
snp_abudance_nonclonal_bps$Source.Host <- as.factor(snp_abudance_nonclonal_bps$Source.Host)
# Set SNP factor levels to AGCT for model reliability 
snp_abudance_nonclonal_bps[,names(snp_abudance_nonclonal_bps) != "Source.Host"] <- lapply(snp_abudance_nonclonal_bps[,names(snp_abudance_nonclonal_bps) != "Source.Host"],factor,levels=c("C","G","T","A"))

##################################################### Make Models ##################################################### 
input_all <- list(amr_class_nonclonal= amr_class_nonclonal,
                  amr_gene_nonclonal = amr_gene_nonclonal,
                  pv_coll_nonclonal = pv_coll_nonclonal,
                  igr_coll_nonclonal = igr_coll_nonclonal,
                  snp_abudance_nonclonal = snp_abudance_nonclonal
)
input_bps <- list(amr_class_nonclonal_bps = amr_class_nonclonal_bps,
                  amr_gene_nonclonal_bps = amr_gene_nonclonal_bps,
                  pv_coll_nonclonal_bps = pv_coll_nonclonal_bps,
                  igr_coll_nonclonal_bps = igr_coll_monoclonal_bps,
                  snp_abudance_nonclonal_bps = snp_abudance_nonclonal_bps
)

# Initialise model performance stats dataframes #####
prediction_overall_all <- data.frame()
prediction_class_all <- data.frame()

###### All Host Models ####
i = 1
for (x in input_all) {
  
  # Hyperparameter testing
  control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", p=0.75, savePredictions="final", classProbs = TRUE)
  tunegrid <- expand.grid(.mtry = c(sqrt(ncol(x))))
  rf_random <- train(Source.Host~., data=x, method="rf", metric="Kappa", tuneGrid = tunegrid, trControl=control, importance=TRUE,ntree=500)
  
  #### Save model ####
  model_filename <- paste(names(input_all)[i], "_model.rds", sep="")
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
}

######## BPS Models #########
i = 1
for (x in input_bps) {
  # Hyperparameter testing
  control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", p=0.75, savePredictions="final", classProbs = TRUE)
  tunegrid <- expand.grid(.mtry = c(sqrt(ncol(x))))
  rf_random <- train(Source.Host~., data=x, method="rf", metric="Kappa", tuneGrid = tunegrid, trControl=control, importance=TRUE,ntree=500)
  
  #### Save model ####
  model_filename <- paste( names(input_bps)[i], "_model.rds", sep="")
  saveRDS(rf_random, model_filename)
  
  ## Get Prediction scores ####
  prediction <- predict.train(rf_random, newdata = x)
  prediction_overall <- as.data.frame(as.matrix(confusionMatrix(prediction, x$Source.Host), what = "overall")) %>%
    rownames_to_column()
  prediction_overall$model=names(input_bps)[i]
  prediction_overall_all <- rbind(prediction_overall_all,prediction_overall)
  prediction_class <- as.data.frame(as.matrix(confusionMatrix(prediction, x$Source.Host), what = "classes")) %>%
    rownames_to_column() %>%
    pivot_longer(cols = !rowname, names_to = "host", values_to="score")
  prediction_class$model=names(input_bps)[i]
  prediction_class_all <- rbind(prediction_class_all, prediction_class)
  
  i = i + 1
}

# save prediction scores here ####
write.table(prediction_overall_all, file="prediction_overall_all.csv", quote=FALSE, sep=",")
write.table(prediction_class_all, file="prediction_class_all.csv", quote=FALSE, sep=",")
