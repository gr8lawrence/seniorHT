library(tidyverse)
library(caret)
library(ranger) # fast implementation of random forest
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
  
  ## Run the RF (with 1500 trees) and SVM on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    learning.set <- studies.min.1[[j]]
    
    # Random Forest
    rf <- ranger::ranger(class ~ ., data = learning.set, num.trees = 1500, probability = TRUE)
    
    # We store the models into the list
    final_models[[5 * (i - 1) + j]] <- rf
    
    pred <- predict(rf, validation.set)
    pred_column <- factor(ifelse(pred$predictions[,1] > 0.5, "basal", "classical"),
                          levels = c("basal", "classical"))
    confusion.mat <- confusionMatrix(data = pred_column,
                                     reference = truth)
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
    learning.results[5 * (i - 1) + j,] <-
      c(
        studies.names.min.1[j],
        studies.names[i],
        signif(accu, digits = 4),
        signif(sen, digits = 4),
        signif(spe, digits = 4)
      )
  }
  
  ## Run the rf on combined datasets
  learning.set <- studies.min.1[[1]]
  for (k in 1:(length(studies.min.1) - 1)) {
    learning.set <- rbind(learning.set, studies.min.1[[1 + k]])
  }
  
  rf <- ranger::ranger(class ~ ., data = learning.set, num.trees = 1500, probability = TRUE)
  
  final_models[[5 * i]] <- rf
  
  pred <- predict(rf, validation.set)
  pred_column <- factor(ifelse(pred$predictions[,1] > 0.5, "basal", "classical"), 
                        levels = c("basal", "classical"))
  confusion.mat <- confusionMatrix(data = pred_column, 
                                   reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[5 * i,] <-
    c(
      paste("comb. minus", studies.names[i],
            sep = " "),
      studies.names[i],
      signif(accu, digits = 4),
      signif(sen, digits = 4),
      signif(spe, digits = 4)
    )

}

# We compile a list of the outputs that contains the models as well as the result table
final_results <- list("models" = final_models, "results" = learning.results)

print(learning.results, n = length(studies)^2)

write.table(learning.results, file = "/nas/longleaf/home/tianyi96/learning_results_RT_rf.csv")

save(x = final_results, file = "learning_result_RT_rf.Rdata")
save.image()
