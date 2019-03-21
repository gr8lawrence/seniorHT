library(tidyverse)
library(caret)
library(RMTL) # package containing multi-task learning functions
source("file_loading.R")
source("preprocessing_fxns.R")
source("rank_transform_fxns.R")
set.seed(100)

## In this program we explore the RMTL package
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
learning.results <- tibble(testing.set = rep('NA', length(studies)),
                           lambda.1 = rep(0, length(studies)),
                           lambda.2 = rep(0, length(studies)),
                           accuracy = rep(0, length(studies)),
                           sensitivity = rep(0, length(studies)),
                           specificity = rep(0, length(studies)))

# We create a list to hold our final model and a list to hold the final prediction results (for ROC curves)
final_models <- list()
final_preds <- list()

## Below is our loop for multi-task learning (MTL)
# Because this is the multi-task learning, we can only perform the hold-one-out task.
# We need to group the feature matrices and ground truth vectors (classification vectors) into two lists
for (i in 1:num.studies) {
  training_set <- studies[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:df.length])
  testing_set_names <- studies.names[i]
  truth <- studies[[i]]$class # The ground truth vector

  # We make the lists of feature matrices and truth vectors (must be -1 or 1 for both catogories)
  data_matrix_list <- lapply(training_set, function(x) data.matrix(x[ ,1:df.length]))
  class_vec <- lapply(training_set, function(x) ifelse(x$class == "basal", 1, -1))

  # We start with cross-validation to find the optimal parameters for lambda_1 and lambda_2 (using loop) in the L21 method
  # We include lambda_2 = 0 (equivalent to penalized logistic regression) to show if there is any
  # performance increase by using multi-task learning
  log_2_lam_2_range <- seq(-5, 5, 1)
  lam2s <- tibble(lam2 = c(0, 2^log_2_lam_2_range),
                best.lam1 = rep(0, length(log_2_lam_2_range) + 1),
                error = rep(0, length(log_2_lam_2_range) + 1))
  k <- 0 #index of the row needs to be rcorded
  for (j in lam2s$lam2) {
    k <- k + 1
    cross_val <- cv.MTC_L21(data_matrix_list, class_vec, nfolds = 20,
                            lam1 = 2^seq(-5, 5, 1),
                            lam2 = ifelse(j == 0, 0, 2^j))
    lam2s$best.lam1[k] <- cross_val$lam1.min
    lam2s$error[k] <- min(cross_val$cvm)
  }

  ind_min_cve <- which(lam2s$error == min(lam2s$error))
  lam1_best <- lam2s$best.lam1[ind_min_cve][1]
  lam2_best <- lam2s$lam2[ind_min_cve][1]

  # With the selected lambda_1 and lambda_2 we fit the MTL L21 model;
  # Keep the result when lam2 equals 0. Also put lam2 = 0 into the loop.
  mtl_model <- MTC_L21(data_matrix_list, class_vec,
                  lam1 = lam1_best,
                  lam2 = lam2_best)
  final_models[[i]] <- mtl_model

  # We predict on the test set and output the measurements of learning performance
  # The test set needs to be a list
  pred <- predict(mtl_model, list(testing_set))
  final_preds[[i]] <- pred
  pred_values <- factor(ifelse(pred[[1]] > 0.5, "basal", "classical"), levels = c("basal", "classical"))
  confusion.mat <- confusionMatrix(data = pred_values, reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[i,] <- c(testing_set_names,
                            lam1_best,
                            lam2_best,
                            signif(accu, digits = 3),
                            signif(sen, digits = 3),
                            signif(spe, digits = 3))

}

## We want to compare the performance of MTL with that of the penalized logistic regression model (lambda_1 = lambda_2 = 0)
## We create an empty tibble first to hold the leanring results
learning.results.plr <- tibble(testing.set = rep('NA', length(studies)), 
                               lambda.1 = 0,
                               lambda.2 = 0,
                               accuracy = rep(0, length(studies)), 
                               sensitivity = rep(0, length(studies)), 
                               specificity = rep(0, length(studies)))
# We create a list to hold our final model and a list to hold the final prediction results (for ROC curves)
final_models_plr <- list()
final_preds_plr <- list()
## Below is our loop for multi-task learning (MTL)
# Because this is the multi-task learning, we can only perform the hold-one-out task.
# We need to group the feature matrices and ground truth vectors (classification vectors) into two lists
for (i in 1:num.studies) {
  training_set <- studies[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:df.length])
  testing_set_names <- studies.names[i]
  truth <- studies[[i]]$class # The ground truth vector
  # We make the lists of feature matrices and truth vectors (must be -1 or 1 for both catogories)
  data_matrix_list <- lapply(training_set, function(x) data.matrix(x[ ,1:df.length]))
  class_vec <- lapply(training_set, function(x) ifelse(x$class == "basal", 1, -1))
  # Withlambda_1 = lambda_2  = 0we fit the MTL L21 model;
  # Keep the result when lam2 equals 0. Also put lam2 = 0 into the loop.
  mtl_model <- MTC_L21(data_matrix_list, class_vec, 
                       lam1 = 0,
                       lam2 = 0)
  final_models_plr[[i]] <- mtl_model
  # We predict on the test set and output the measurements of learning performance
  # The test set needs to be a list
  pred <- predict(mtl_model, list(testing_set))
  final_preds_plr[[i]] <- pred
  pred_values <- factor(ifelse(pred[[1]] > 0.5, "basal", "classical"), levels = c("basal", "classical"))  
  confusion.mat <- confusionMatrix(data = pred_values, reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results.plr[i,] <- c(testing_set_names,
                                0,
                                0,
                                signif(accu, digits = 3),
                                signif(sen, digits = 3),
                                signif(spe, digits = 3))
  
}

# We create a list to hold all models and learning results
final_results <- list("models" = final_models,
                      "results" = learning.results)
final_results_plr <- list("models" = final_models_plr,
                          "results" = learning.results.plr)
print(learning.results, n = length(studies))
print(learning.results.plr, n = length(studies))
write.table(learning.results, file = "result_tables/MTL_RT_results.csv")
save(x = final_results, file = "models_and_predictions/MTL_RT_results.Rdata")
save(x = final_results_plr, file = "models_and_predictions/MTL_RT_results_no_gp.Rdata")
save(x = final_preds, file = "models_and_predictions/MTL_RT_predictions.Rdata")
save(x = final_preds_plr, file = "models_and_predictions/MTL_RT_predictions_no_gp.Rdata")
save.image()

