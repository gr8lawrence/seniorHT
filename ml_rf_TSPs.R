library(tidyverse)
library(caret)
library(ranger) # package containing rf
source("file_loading.R")
source("preprocessing_fxns.R")
source("rank_transform_fxns.R")
source("TSPs_fxns.R")
set.seed(100)

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

expression.df <- lapply(ranked.studies.df, extractData, 
                        common.gene.names = commonGeneNames) %>% 
  reorderGeneBySignificance(common.gene.names = commonGeneNames)

# This variable indicates how many genes of top differential expression are kept to make TSPs
# this will give us a maximum of (df.length)*(df.length - 1)/2 variables in TSPs
numGenesToKeep <- 100

# Make the TSPs 
geneSubsetNames <- colnames(expression.df[[1]])[1:numGenesToKeep]
studies <- lapply(expression.df, transposeTibble, 
                  df.length = numGenesToKeep) %>% 
  makeTSPs(gene.subset.names = geneSubsetNames, 
           df.length = numGenesToKeep)

# This is the same as the total number of TSPs in the final training data
new.df.length <- dim(studies[[1]])[2]

# Append the classification to the training data
for (i in 1:num.studies) {
  studies[[i]] <- add_column(studies[[i]], class = labels[[i]]) 
}


# We create an empty tibble first to hold the leanring results
learning.results <- tibble(training.set = rep('NA', length(studies)^2), 
                           testing.set = rep('NA', length(studies)^2), 
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))

# We also create a list to store the final learning models
final_models <- list()
final_preds <- list()

## Below is our loop for random forest (rf)

for (i in 1:num.studies) {
  studies.min.1 <- studies[-i]
  studies.names.min.1 <- studies.names[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:new.df.length])
  truth <- studies[[i]]$class # The truth vector
  
  ## Run the rf on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    training_set <- studies.min.1[[j]]
    
    # We train the random forest model on the training set.
    # Because the high dimensinality of our data, to avoid crashes we turn on the memory saving function
    rf <- ranger::ranger(class~., data = training_set, probability = TRUE, save.memory = TRUE)
    
    # We store the models to the list
    final_models[[5 * (i - 1) + j]] <- rf
    
    # We output the relevant values to measure the results of learning by our rf
    # The predicted subtype is assigned based on probability (of being a "basal"); 
    # if the probability is greater than 0.5 we assign it to the "basal" class;
    # otherwise, we assign it to the "classical" class.
    
    pred <- predict(rf, testing_set)
    
    # We also store the predictions to the list for plotting ROC curves
    final_preds[[5 * (i - 1) + j]] <- pred
    
    pred_column <-  ifelse(pred$predictions[ ,1] > 0.5, "basal", "classical")  %>% 
      factor(levels = c("basal", "classical"))  
    
    confusion.mat <- confusionMatrix(data = pred_column,
                                     reference = truth)
    
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
    learning.results[5 * (i - 1) + j,] <-
      c(
        studies.names.min.1[j],
        studies.names[i],
        signif(accu, digits = 3),
        signif(sen, digits = 3),
        signif(spe, digits = 3)
      )
  }
  
  ## Run the rf on combined datasets
  training_set <- studies.min.1[[1]]
  
  for (k in 1:(length(studies.min.1) - 1)) {
    training_set <- rbind(training_set, studies.min.1[[1 + k]])
  }
  
  training_set[ ,1:new.df.length] %>% apply(2, scale) # rescaling the columns of the combined dataset
  
  rf <- ranger::ranger(class~., data = training_set, probability = TRUE, save.memory = TRUE)
  final_models[[5 * i]] <- rf
  pred <- predict(rf, testing_set)
  final_preds[[5 * i]] <- pred
  pred_column <- ifelse(pred$predictions[ ,1] > 0.5, "basal", "classical") %>%
    factor(levels = c("basal", "classical"))
  
  confusion.mat <- confusionMatrix(data = pred_column,
                                   reference = truth)
  
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  
  learning.results[5 * i,] <-
    c(
      paste("Combined minus", studies.names[i], sep = " "),
      studies.names[i],
      signif(accu, digits = 3),
      signif(sen, digits = 3),
      signif(spe, digits = 3)
    )
  
}

# We compile a list of the outputs that contains the models as well as the result table
final_results <- list("models" = final_models, "results" = learning.results)

print(learning.results, n = length(studies)^2)

write.table(learning.results, file = "/nas/longleaf/home/tianyi96/TSPs_svm_results.csv")

save(x = final_results, file = "TSPs_rf_results.Rdata")
save(x = final_preds, file = "TSPs_rf_predictions.Rdata")
save.image()

