library(tidyverse)
library(caret)
library(e1071) # package containing svm
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
numGenesToKeep <- 10

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
                           gamma = rep(0, length(studies)^2),
                           cost = rep(0, length(studies)^2),
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))
# We also create an empty list to hold the models
final_models <- list()
# Below is our loop for svm
for (i in 1:num.studies) {
  studies.min.1 <- studies[-i]
  studies.names.min.1 <- studies.names[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:new.df.length])
  # The truth vector
  truth <- studies[[i]]$class 
  ## Run the SVM on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    training_set <- studies.min.1[[j]]
    # Now because our data frame (tibble) of training data with TSPs is giant, when it is passed to the svm 
    # too many recursion may result and this will overflow the protection stack
    # so we separate the data frame into a data matrix and a classification vector
    dat_mat <- data.matrix(training_set[ ,1:new.df.length])
    class <- training_set$class
    # We go through parameter tuning at first to search for the optimal cost C and tuning parameter gamma for the radial basis kernel
    # Use a grid search to find the best (C, gamma) pair
    # Cross-validation generate 40 SVMs  
    tune <- tune.svm(dat_mat, class, 
                     gamma = 10^seq(-7, 3, 2), 
                     cost = 10^seq(-3, 7, 2), 
                     scale = FALSE) 
    # This line is for adjusting class weights as the classes are not balanced (have the same number of observations)
    # in the training data
    classical.basal.ratio = sum(class == "classical")/sum(class == "basal")
    # We put the best (C, gamma) values and adjusted class weights into our final svm model
    supp.vec <- svm(dat_mat, class, 
                    class.weights = c("basal" = classical.basal.ratio, "classical" = 0.1), 
                    gamma = tune$best.parameters$gamma, 
                    cost = tune$best.parameters$cost, 
                    probability = TRUE,
                    scale = FALSE)
    # We store the model to the list
    final_models[[5 * (i - 1) + j]] <- supp.vec
    # We output the relevant values to measure the results of learning by our svm
    pred <- predict(supp.vec, testing_set, probability = TRUE)
    confusion.mat <- confusionMatrix(data = pred,
                                     reference = truth)
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
     learning.results[5 * (i - 1) + j,] <- c(studies.names.min.1[j],
                                             studies.names[i],
                                             tune$best.parameters$gamma,
                                             tune$best.parameters$cost,
                                             signif(accu, digits = 3),
                                             signif(sen, digits = 3),
                                             signif(spe, digits = 3))
  }
  # Run the SVM on combined datasets
  training_set <- studies.min.1[[1]]
  for (k in 1:(length(studies.min.1) - 1)) {
    training_set <- rbind(training_set, studies.min.1[[1 + k]])
  }
  dat_mat <- data.matrix(training_set[ ,1:new.df.length])
  dat_mat %>% apply(2, scale) # rescaling the columns of the combined dataset
  class <- training_set$class
  tune <- tune.svm(dat_mat, class,
                   gamma = 10^seq(-7, 3, 2), 
                   cost = 10^seq(-3, 7, 2),
                   scale = FALSE) 
  classical.basal.ratio = sum(class == "classical")/sum(class == "basal")
  supp.vec <- svm(dat_mat, class,
                  class.weights = c("basal" = classical.basal.ratio, "classical" = 1), 
                  gamma = tune$best.parameters$gamma,
                  cost = tune$best.parameters$cost,
                  probability = TRUE,
                  scale = FALSE)
  final_models[[5 * i]] <- supp.vec
  pred <- predict(supp.vec, testing_set, probability = TRUE)
  confusion.mat <- confusionMatrix(data = pred,
                                   reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[5 * i,] <- c(paste("Combined minus", studies.names[i], sep = " "),
                                studies.names[i],
                                tune$best.parameters$gamma,
                                tune$best.parameters$cost,
                                signif(accu, digits = 3),
                                signif(sen, digits = 3),
                                signif(spe, digits = 3))
}
# We compile a list of the outputs that contains the models as well as the result table
final_results <- list("models" = final_models, 
                      "results" = learning.results)
# Print the learning results
print(learning.results, n = 2*length(studies)^2)
write.table(learning.results, file = "result_tables/SVM_TSPs_results.csv")
save(x = final_results, file = "models_and_predictions/SVM_TSPs_results.Rdata")
save.image()

