library(tidyverse)
library(caret)
library(ncvreg)
library(parallel)
source("file_loading.R")
source("preprocessing_fxns.R")
source("rank_transform_fxns.R")
source("TSPs_fxns.R")
set.seed(100)

# Make clusters for paralleling
no_cores <- 2
cl <- makeCluster(no_cores)

num.studies <- length(studies.df) # Number of datasets
studies.df <- lapply(studies.df, geneSort)
commonGeneNames <- getGeneSymbols(studies.df[[1]])
for (i in 2:num.studies) {
  commonGeneNames <- intersect(commonGeneNames, getGeneSymbols(studies.df[[i]]))
}
numCommonGenes <- length(commonGeneNames) # Number of common genes

ranked.studies.df  <- lapply(studies.df, matchCommonGenes) %>% 
  lapply(rankTransform) 

# The labels for patients with a label
labels <- lapply(ranked.studies.df, function(x){ x$sampInfo$cluster.MT[!is.na(x$sampInfo$cluster.MT)] })

studies <- lapply(ranked.studies.df, extractData, 
                        common.gene.names = commonGeneNames) %>% 
  reorderGeneBySignificance(common.gene.names = commonGeneNames)

# This variable indicates how many genes of top differential expression are kept to make TSPs
# this will give us a maximum of (df.length)*(df.length - 1)/2 variables in TSPs
df.length <- ceiling(numCommonGenes/5)

# Append the classification to the training data
for (i in 1:num.studies) {
  studies[[i]] <- add_column(studies[[i]][ , 1:df.length], class = labels[[i]]) 
}


# We create an empty tibble first to hold the leanring results
learning.results <- tibble(training.set = rep('NA', length(studies)^2), 
                           testing.set = rep('NA', length(studies)^2),
                           strategy = rep('NA', length(studies)^2),
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))

# We also create a list to store the final learning models
final_models <- list()
final_preds <- list()

## Below is our loop for our penalized logistic regression model
for (i in 1:num.studies) {
  studies.min.1 <- studies[-i]
  studies.names.min.1 <- studies.names[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:df.length])
  truth <- studies[[i]]$class # The truth vector
  
  ## Run the rf on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    training_set <- studies.min.1[[j]]
    
    # Separate the feature matrix and the labels
    feature_matrix <- data.matrix(training_set[ ,1:df.length])
    class_labels <- training_set$class
    numeric_labels <- ifelse(class_labels == "basal", 1, -1)
    max_num_of_folds <- length(class_labels)
    
    # Cross-validation to find the best penalty parameter
    cvfit <- ncvreg::cv.ncvreg(feature_matrix, numeric_labels,
                               family = "binomial", 
                               cluster = cl
                               # penalty = "lasso", 
                               # nfolds = 30, 
                               # seed = 100, 
                               # max.iter = 1000000
                               ) 
    
    # Fit the regression model based on the best lambda
    fit <- ncvreg::ncvreg(feature_matrix, numeric_labels, 
                          family = "binomial", 
                          lambda = cvfit$lambda.min, 
                          penalty = "SCAD", 
                          seed = 100, 
                          max.iter = 1000000) 
    final_models[[5 * (i - 1) + j]] <- fit
    
    # Predict on the hold-out set
    pred <- predict(fit, X = data.matrix(testing_set), type = "class")
    final_preds[[5 * (i - 1) + j]] <- pred
    confusion.mat <- confusionMatrix(data = factor(ifelse(pred == 1, "basal", "classical"), levels = c("basal", "classical")),
                                     reference = truth)
    
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
    learning.results[5 * (i - 1) + j,] <- c(studies.names.min.1[j],
                                            studies.names[i],
                                            "Single",
                                            signif(accu, digits = 3),
                                            signif(sen, digits = 3),
                                            signif(spe, digits = 3))
    
  }
  
  ## Run the plr on combined datasets
  training_set <- studies.min.1[[1]]
  
  for (k in 1:(length(studies.min.1) - 1)) {
    training_set <- rbind(training_set, studies.min.1[[1 + k]])
  }
  
  training_set[ ,1:df.length] %>% apply(2, scale) # rescaling the columns of the combined dataset
  # Separate the feature matrix and the labels
  feature_matrix <- data.matrix(training_set[ ,1:df.length])
  class_labels <- training_set$class
  numeric_labels <- ifelse(class_labels == "basal", 1, -1)
  max_num_of_folds <- length(class_labels)
  
  # Cross-validation to find the best penalty parameter
  cvfit <- ncvreg::cv.ncvreg(feature_matrix, numeric_labels,
                             family = "binomial", 
                             cluster = cl
                             # penalty = "lasso", 
                             # nfolds = 30, 
                             # seed = 100, 
                             # max.iter = 1000000
  ) 
  
  # Fit the regression model based on the best lambda
  fit <- ncvreg::ncvreg(feature_matrix, numeric_labels, 
                        family = "binomial", 
                        lambda = cvfit$lambda.min, 
                        penalty = "SCAD", 
                        seed = 100, 
                        max.iter = 1000000) 
  final_models[[5 * i]] <- fit
  
  # Predict on the hold-out set
  pred <- predict(fit, X = data.matrix(testing_set), type = "class")
  final_preds[[5 * i]] <- pred
  confusion.mat <- confusionMatrix(data = factor(ifelse(pred == 1, "basal", "classical"), levels = c("basal", "classical")),
                                   reference = truth)
  
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  
  learning.results[5 * i,] <- c(paste("Combined minus", studies.names[i], sep = " "),
                                studies.names[i],
                                "Combined",
                                signif(accu, digits = 3),
                                signif(sen, digits = 3),
                                signif(spe, digits = 3))
}

# We compile a list of the outputs that contains the models as well as the result table
final_results <- list("models" = final_models, 
                      "results" = learning.results)
print(learning.results, n = length(studies)^2)
write.table(learning.results, file = "result_tables/PLR_RT_results.csv")
save(x = final_results, file = "models_and_predictions/PLR_RT_results.Rdata")
save(x = final_preds, file = "models_and_predictions/PLR_RT_predictions.Rdata")
save.image()

