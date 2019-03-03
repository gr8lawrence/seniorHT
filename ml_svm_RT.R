library(tidyverse)
library(caret)
library(e1071) # svm
source("file_loading.R")
source("preprocessing_fxns.R")
source("rank_transform_fxns.R")
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

studies <- lapply(ranked.studies.df, extractData, 
                  common.gene.names = commonGeneNames) %>% 
  reorderGeneBySignificance(common.gene.names = commonGeneNames)

# This variable indicates how many genes of top differential expression are kept
df.length <- ceiling(length(commonGeneNames)/5)

# Append the classification to the training data
for (i in 1:num.studies) {
  studies[[i]] <- add_column(studies[[i]][ , 1:df.length], class = labels[[i]]) 
}


# We create an empty tibble first to hold the leanring results
learning.results <- tibble(learning.set = rep('NA', length(studies)^2), 
                           validation.set = rep('NA', length(studies)^2), 
                           gamma = rep(0, length(studies)^2),
                           cost = rep(0, length(studies)^2),
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))
len <- length(studies)
final_models <- list()

for (i in 1:len) {
  studies.min.1 <- studies[-i]
  studies.names.min.1 <- studies.names[-i]
  validation.set <- studies[[i]][, 1:df.length]
  truth <- studies[[i]]$class # The truth vector
  
  ## Run the SVM on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    learning.set <- studies.min.1[[j]]
    
    # We go through parameter tuning at first to search for the optimal cost C and tuning parameter gamma for the radial basis kernel
    # Use a grid search to find the best (C, gamma) pair
    tune <- tune.svm(class ~ ., data = learning.set, 
                     gamma = 10^seq(-5, 3, 2), 
                     cost = 10^seq(-3, 5, 2)) 
    classical.basal.ratio = sum(learning.set$class == "classical")/sum(learning.set$class == "basal")
    supp.vec <- svm(class ~ ., data = learning.set, 
                    class.weights = c("basal" = classical.basal.ratio, "classical" = 0.1), 
                    gamma = tune$best.parameters$gamma, 
                    cost = tune$best.parameters$cost, 
                    probability = TRUE)
    
    # We store the models into the list
    final_models[[5 * (i - 1) + j]] <- supp.vec
    
    pred <- predict(supp.vec, validation.set, probability = TRUE)
    confusion.mat <- confusionMatrix(data = pred,
                                     reference = truth)
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
    learning.results[5 * (i - 1) + j,] <-
      c(
        studies.names.min.1[j],
        studies.names[i],
        tune$best.parameters$gamma,
        tune$best.parameters$cost,
        signif(accu, digits = 3),
        signif(sen, digits = 3),
        signif(spe, digits = 3)
      )
  }
  
  ## Run the SVM on combined datasets
  learning.set <- studies.min.1[[1]]
  for (k in 1:(length(studies.min.1) - 1)) {
    learning.set <- rbind(learning.set, studies.min.1[[1 + k]])
  }
  tune <- tune.svm(class ~ ., data = learning.set, 
                   gamma = 10^seq(-5, 3, 2), 
                   cost = 10^seq(-3, 5, 2)) 
  classical.basal.ratio = sum(learning.set$class == "classical")/sum(learning.set$class == "basal")
  supp.vec <- svm(class ~ ., data = learning.set, 
                  class.weights = c("basal" = classical.basal.ratio, "classical" = 1), 
                  gamma = tune$best.parameters$gamma,
                  cost = tune$best.parameters$cost,
                  probability = TRUE)
  
  final_models[[5 * i]] <- supp.vec
  
  pred <- predict(supp.vec, validation.set, probability = TRUE)
  confusion.mat <- confusionMatrix(data = pred,
                                   reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[5 * i,] <-
    c(
      paste("Combined minus", studies.names[i], sep = " "),
      studies.names[i],
      tune$best.parameters$gamma,
      tune$best.parameters$cost,
      signif(accu, digits = 3),
      signif(sen, digits = 3),
      signif(spe, digits = 3)
    )
  
}
final_results <- list("models" = final_models, 
                      "results" = learning.results)
print(learning.results, n = length(studies)^2)
write.table(learning.results, file = "result_tables/learning_results_SVM_RT.csv")
save(x = final_results, file = "models_and_predictions/SVM_RT_results.Rdata")
save.image()
