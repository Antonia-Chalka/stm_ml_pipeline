library(tidyverse)
library(randomForest)
library(caret)
library(missForest)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least x argumentS: if not, return an error
if (length(args)!=5) {
  length(args)
  stop("Exactly 5 arguments must be supplied (input file).n", call.=FALSE)
} 

### AMR #
# Class
amr_class_nonclonal <- read.table(args[1], header=TRUE, row.names = 1)
amr_class_nonclonal$Source.Host <- as.factor(amr_class_nonclonal$Source.Host)
# Gene
amr_gene_nonclonal <- read.table(args[2], header=TRUE, row.names = 1)
amr_gene_nonclonal$Source.Host <- as.factor(amr_gene_nonclonal$Source.Host)

### PV #
pv_coll_nonclonal <- read.table(args[3], header=TRUE, sep="\t", row.names = 1)NT
pv_coll_nonclonal$Source.Host <- as.factor(pv_coll_nonclonal$Source.Host)

### IGR #
igr_coll_nonclonal <- read.table(args[4], header = TRUE, sep="\t", row.names = 1) 
igr_coll_nonclonal$Source.Host <- as.factor(igr_coll_nonclonal$Source.Host)

### SNP #
snp_abudance_nonclonal <- read.table(args[5], header = TRUE, sep="\t", row.names = 1)
snp_abudance_nonclonal$Source.Host <- as.factor(snp_abudance_nonclonal$Source.Host)
# Set SNP factor levels to AGCT for model reliability 
snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"] <- lapply(snp_abudance_nonclonal[,names(snp_abudance_nonclonal) != "Source.Host"],factor,levels=c("C","G","T","A"))


##### Make list of models #
all_models <- list(amr_class_nonclonal= amr_class_nonclonal,
                   amr_gene_nonclonal = amr_gene_nonclonal,
                   pv_coll_nonclonal = pv_coll_nonclonal,
                   igr_coll_nonclonal = igr_coll_nonclonal,
                   snp_abudance_nonclonal = snp_abudance_nonclonal
)

# Get all assembly names and their observed host (from a dataset that includes all assemblies)
predictions_df <- pv_coll_nonclonal %>% select(Source.Host) %>% rownames_to_column(var="Assembly")

# For all models
i=1
for (model_test in all_models) {
  # Load Model
  model_filename <- paste("./models/", names(all_models)[i], "_model.rds", sep="")
  rf_random <- readRDS(model_filename) 
  
  # Impute Missing columns
  model_train <- rf_random$trainingData %>% select(!.outcome) %>% mutate(Source.Host="N") # Get training dataset
  pv_diff <- setdiff(names(model_train),names(model_test)) # Check if there is a difference between training and testing datasets
  model_test[,pv_diff] <- NA # Add missing columns to testing dataset and set to NA
  model_all <- rbind(model_test, model_train) # Combine training and testing dataset for imputation of missing values
  index <- 1:ncol(model_all) 
  model_all[ , index] <- lapply(model_all[ , index], as.factor) # Convert from numeric (0,1) to factor so that imputation is based on categorical values
  
  model_imp <- missForest(model_all, xtrue=model_train, verbose=TRUE) # Impute missing data
  
  model_test <- subset(model.imp$ximp, row.names(model_all) %in% row.names(model_test)) # Isolate testing data only
  index <- 1:ncol(model_test)
  model_test[ , index] <- lapply(model_test[ , index],function(x) as.numeric(as.character(x)))  # Convert from factor to numeric
  
  # Predict
  predict_prob <- predict.train(rf_random, newdata = model_test, type="prob") 
  
  # Saving predictions
  predict_prob$prediction <- names(predict_prob)[1:4][apply(predict_prob[,1:4], 1, which.max)] # Get highest scored prediction
  predictions_df <- predictions_df %>%
    full_join(rownames_to_column(predict_prob, var="Assembly") ) 
  i=i+1
}

write.table(predictions_df, file ="prediction_all_assemblies.tsv", sep="\t", row.names = FALSE, quote = FALSE)
