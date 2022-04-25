library(tidyverse)
library(randomForest)
library(caret)
set.seed(100)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
    length(args)
    stop("Exactly 1 arguments must be supplied (input file).n", call.=FALSE)
} 

# Load Input Data
input_data <- read.table(args[1], header=TRUE, row.names = 1)
input_data$Source.Host <- as.factor(input_data$Source.Host)
if (!is.numeric(input_data[3,3])) {
    input_data[,names(input_data) != "Source.Host"] <- lapply(input_data[,names(input_data) != "Source.Host"],factor,levels=c("C","G","T","A"))
}
data_name <- tools::file_path_sans_ext(basename(args[1]))

write.table(input_data, file="./test.tsv", sep="\t")


# Hyperparameter testing
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid", p=0.75, savePredictions="final", classProbs = TRUE)
tunegrid <- expand.grid(.mtry = c(sqrt(ncol(input_data))))

rf_random <- train(Source.Host~., data=input_data, method="rf", metric="Kappa", tuneGrid = tunegrid, trControl=control, importance=TRUE, ntree=500)

# Save model
model_filename <- paste(data_name, "_model.rds", sep="")
saveRDS(rf_random,model_filename)

# Get Prediction scores
prediction <- predict.train(rf_random, newdata = input_data)
prediction_overall <- as.data.frame(as.matrix(confusionMatrix(prediction, input_data$Source.Host), what = "overall")) %>%
    rownames_to_column()
prediction_overall$model=data_name
write.table(prediction_overall, file=paste(data_name, "_prediction_overall.csv", sep=""), quote=FALSE, sep=",")

prediction_class <- as.data.frame(as.matrix(confusionMatrix(prediction, input_data$Source.Host), what = "classes")) %>%
    rownames_to_column() %>%
    pivot_longer(cols = !rowname, names_to = "host", values_to="score")
prediction_class$model=data_name
write.table(prediction_class, file=paste(data_name, "_prediction_class.csv", sep=""), quote=FALSE, sep=",")
