library(tidyverse)
library(randomForest)
library(caret)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=5) {
  length(args)
  stop("Exactly 5 arguments must be supplied (input file).n", call.=FALSE)
} 

############################# Loading Data ####
amr_class <- read.table(args[1], header=TRUE, row.names = 1)
amr_gene <- read.table(args[2], header=TRUE, row.names = 1)
pv <- read.table(args[3], header=TRUE, sep="\t", row.names = 1)
igr <- read.table(args[4], header = TRUE, sep="\t", row.names = 1)
snp_abudance <- read.table(args[5], header = TRUE, sep="\t", row.names = 1, colClasses = "character")

# Set SNP factor levels to AGCT for model reliability 
snp_abudance[,names(snp_abudance) != "Source.Host"] <- lapply(snp_abudance[,names(snp_abudance) != "Source.Host"],factor,levels=c("C","G","T","A"))

# Make list of models
all_models <- list( amr_class= amr_class,
                    amr_gene = amr_gene,
                    pv = pv,
                    igr = igr,
                    snp_abudance = snp_abudance
)

############################# Predict for all-host models ####
predictions <- list()
i=1
for (model_test in all_models) {
  # Load Model
  model_filename <- paste("./models/", names(all_models)[i], "_all_model.rds", sep="")
  rf_random <- readRDS(model_filename) 
  
  # Impute Missing columns
  model_train <- rf_random$trainingData %>% select(!.outcome) %>% mutate(Source.Host="N") # Get training dataset
  col_diff <- setdiff(names(model_train),names(model_test)) # Check if there is a difference between training and testing datasets
  
  model_test[,col_diff] <- 0 # Add missing columns to testing dataset and set to 0
  
  # Predict
  predict_prob <- predict.train(rf_random, newdata = model_test, type="prob") 
  
  # Saving predictions
  predict_prob$prediction <- names(predict_prob)[1:4][apply(predict_prob[,1:4], 1, which.max)] # Get highest scored prediction
  predict_prob$model <-names(all_models[i])
  predict_prob <- rownames_to_column(predict_prob, var="Assembly")
  
  predictions[[i]] <- predict_prob 
  i=i+1
}

predictions_df <- bind_rows(predictions)
write.table(predictions_df, file ="prediction_all.tsv", sep="\t", row.names = FALSE, quote = FALSE)

############################# Predict for livestock-only models ####
predictions <- list()
i=1
for (model_test in all_models) {
  # Load Model
  model_filename <- paste("./models/", names(all_models)[i], "_bps_model.rds", sep="")
  rf_random <- readRDS(model_filename) 
  
  # Impute Missing columns
  model_train <- rf_random$trainingData %>% select(!.outcome) %>% mutate(Source.Host="N") # Get training dataset
  col_diff <- setdiff(names(model_train),names(model_test)) # Check if there is a difference between training and testing datasets
  
  model_test[,col_diff] <- 0 # Add missing columns to testing dataset and set to 0
  
  # Predict
  predict_prob <- predict.train(rf_random, newdata = model_test, type="prob") 
  
  # Saving predictions
  predict_prob$prediction <- names(predict_prob)[1:3][apply(predict_prob[,1:3], 1, which.max)] # Get highest scored prediction
  predict_prob$model <-names(all_models[i])
  predict_prob <- rownames_to_column(predict_prob, var="Assembly")
  
  predictions[[i]] <- predict_prob 
  i=i+1
}

predictions_df <- bind_rows(predictions)
write.table(predictions_df, file ="prediction_bps.tsv", sep="\t", row.names = FALSE, quote = FALSE)

############################# Predict for human-scoring models ####
predictions <- list()
i=1
for (model_test in all_models) {
  # Load Model
  model_filename <- paste("./models/", names(all_models)[i], "_human_model.rds", sep="")
  rf_random <- readRDS(model_filename) 
  
  # Impute Missing columns
  model_train <- rf_random$trainingData %>% select(!.outcome) %>% mutate(Source.Host="N") # Get training dataset
  col_diff <- setdiff(names(model_train),names(model_test)) # Check if there is a difference between training and testing datasets
  
  model_test[,col_diff] <- 0 # Add missing columns to testing dataset and set to 0
  
  # Predict
  predict_prob <- predict.train(rf_random, newdata = model_test, type="prob") 
  
  # Saving predictions
  predict_prob$prediction <- names(predict_prob)[1:2][apply(predict_prob[,1:2], 1, which.max)] # Get highest scored prediction
  predict_prob$model <-names(all_models[i])
  predict_prob <- rownames_to_column(predict_prob, var="Assembly")
  
  predictions[[i]] <- predict_prob 
  i=i+1
}

predictions_df <- bind_rows(predictions)
write.table(predictions_df, file ="prediction_human.tsv", sep="\t", row.names = FALSE, quote = FALSE)
