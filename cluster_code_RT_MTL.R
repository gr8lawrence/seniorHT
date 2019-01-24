library(tidyverse)
library(caret)
library(RMTL) # package containing multi-task learning functions
set.seed(100)

## In this section we explore the RMTL package

dataPath <-
  "/Users/gr8lawrence/Desktop/Senior Honors Thesis/datasets/"
  #"/nas/longleaf/home/tianyi96/datasets_used/" # change this line to your local dataset directory

studies.names <- c("Aguirre-Seq",
                   "Linehan-Seq",
                   "COMPASS",
                   "Moffitt Arrays",
                   "TCGA-PAAD") # Names of the studies in the same order as loading below 

# Load datasets
load(paste(dataPath, "Aguirre_seq_plus.RData", sep = ""))
load(paste(dataPath, "Linehan_Seq_plus.RData", sep = ""))
load(paste(dataPath, "COMPASS.2017_plus.RData", sep = ""))
load(paste(dataPath, "Moffitt_GEO_array_plus.RData", sep = ""))
load(paste(dataPath, "TCGA_PAAD_plus.RData", sep = ""))

# Incorporate all datasets into a single list object
studies.df <- list(Aguirre_seq_plus, 
                   Linehan_Seq_plus,
                   COMPASS.2017_plus,
                   Moffitt_GEO_array_plus,
                   TCGA_PAAD_plus)

num.studies <- length(studies.df) # Number of datasets

# Sort gene in alphanumeric order
geneSort <- function(d) {
  # d is a dataset
  index <- order(d$featInfo$SYMBOL)
  d$ex[,] <- d$ex[index,]
  d$featInfo[,] <- d$featInfo[index,]
  return(d)
}

studies.df <- lapply(studies.df, geneSort)

# Obtain the list of gene symbols from datasets
getGeneSymbols <- function(d) {
  symbols <- d$featInfo$SYMBOL
  return(symbols)
} 

# Find the common genes across all datasets
commonGeneNames <- getGeneSymbols(studies.df[[1]])
for (i in 2:num.studies) {
  commonGeneNames <- intersect(commonGeneNames, getGeneSymbols(studies.df[[i]]))
}

num.common.genes <- length(commonGeneNames) # Number of common genes

# Subset the datasets to ones of common genes
matchCommonGenes <- function(d) {
  index <- match(commonGeneNames, d$featInfo$SYMBOL)
  d$ex <- d$ex[index, ]
  d$featInfo <- d$featInfo[index, ]
  return(d)
}

studies.df <- lapply(studies.df, matchCommonGenes)

rankTransform <- function(d) {
  apply(d$ex, 2, rank)
  return(d)
}

ranked.studies.df <- lapply(studies.df, rankTransform)

# Extract ranked expression from each dataset, remove observations of unknown cancer classification, and write them into tibbles
extractData <- function(d) {
  df <- as_tibble(t(d$ex))
  colnames(df) <- make.names(commonGeneNames) # Get all gene names into correct name formats
  df <- df[!is.na(d$sampInfo$cluster.MT),] # Remove observations with unknown subtype
  return(add_column(df, class = d$sampInfo$cluster.MT[!is.na(d$sampInfo$cluster.MT)]))
}

expression.df <- lapply(ranked.studies.df, extractData)

# Separating expression data for different subtypeså
expression.basal.df <- lapply(expression.df, function(x) x[x$class == "basal", ])
expression.classical.df <- lapply(expression.df, function(x) x[x$class == "classical", ])

# Transforming dataset to matrices so we can conduct statistical tests (Wilcoxon rank sum test)
basal.matrices <- list()
classical.matrices <- list()

for (i in 1:num.studies) {
  basal.matrices[[i]] <- data.matrix(expression.basal.df[[i]][ ,1:num.common.genes])
  classical.matrices[[i]] <- data.matrix(expression.classical.df[[i]][ ,1:num.common.genes])
}

# Build up empty list to contain the p-values of single tests (in one dataset) for each gene
gene.pvals <- list()
for (i in 1:num.studies){
  gene.pvals[[i]] <- tibble(name = make.names(commonGeneNames),
                            pval = rep(0, num.common.genes))
}

for (i in 1:num.studies) {
  for (j in 1:num.common.genes) {
    basal.vector <- basal.matrices[[i]][ ,j]
    classical.vector <- classical.matrices[[i]][ ,j]
    gene.pvals[[i]]$pval[j] <- wilcox.test(basal.vector, classical.vector, exact = FALSE)[["p.value"]]
  }
}

# Compute each gene's sum of log-10 p-values
gene.sum.log.pvals <- tibble(name = make.names(commonGeneNames),
                             pval = rep(0, length(commonGeneNames)))

for (i in 1:num.common.genes) {
  gene.pvals.vector <- vector()
  for (j in 1:num.studies) {
    gene.pvals.vector[j] = gene.pvals[[j]]$pval[i]
  }
  log.10.gene.pvals <- log10(gene.pvals.vector)
  gene.sum.log.pvals$pval[i] <- sum(log.10.gene.pvals)
}

index <- order(gene.sum.log.pvals$pval) # Give the sorted indexes of genes

reorderGene <- function(df) {
  df.colnames <- colnames(df)
  df[ ,1:length(commonGeneNames)] <- data.matrix(df[ ,index])
  colnames(df) <- df.colnames[index]
  return(as_tibble(df))
}

expression.df <- lapply(expression.df, reorderGene)

df.length <- 100
  
#num.common.genes/5 # We keep 1/5 of the total number of genes in the training sets for RMTL

# Append the classification to the training data

studies <- list()

for (i in 1:num.studies) {
  studies[[i]] <- add_column(expression.df[[i]], class = ranked.studies.df[[i]]$sampInfo$cluster.MT[!is.na(ranked.studies.df[[i]]$sampInfo$cluster.MT)]) 
}

# We create an empty tibble first to hold the leanring results
learning.results <- tibble(testing.set = rep('NA', length(studies)), 
                           lambda.1 = rep(0, length(studies)),
                           lambda.2 = rep(0, length(studies)),
                           accuracy = rep(0, length(studies)), 
                           sensitivity = rep(0, length(studies)), 
                           specificity = rep(0, length(studies)))

## Below is our loop for multi-task learning (MTL)
# Because this is the multi-task learning, we can only perform the hold-one-out task.
# We need to group the feature matrices and ground truth vectors (classification vectors) into two lists

#for (i in 1:num.studies) {
  i = 1
  training_set <- studies[-i]
  training_set_names <- studies.names[-i]
  testing_set <- data.matrix(studies[[i]][ ,1:df.length])
  truth <- studies[[i]]$class # The ground truth vector
  
  # We make the lists of feature matrices and truth vectors
  dat_mat <- lapply(training_set, function(x) data.matrix(x[ ,1:df.length]))
  class_vec <- lapply(training_set, function(x) ifelse(x$class == "basal", 1, -1))
    
  # We start with cross-validation to find the optimal parameters for lambda_1 and lambda_2 (using loop) in the L21 method
  log_10_lam_2_range <- seq(-5, 3, 2)
  
  lam2s <- list(lam2 = c(0, 10^log_10_lam_2_range), 
                best.lam1 = rep(0, length(log_10_lam_2_range) + 1),
                error = rep(0, length(log_10_lam_2_range) + 1))
  
  for (j in log_10_lam_2_range) {
    cross_val <- cv.MTC_L21(dat_mat, class_vec, nfolds = 15,
                            lam1 = 10^seq(-5, 3, 2), 
                            lam2 = 10^j)
    lam2s$best.lam1[j] <- cross_val$lam1.min
    lam2s$error[j] <- min(cross_val$cvm)
  }
  
  ind_min_cve <- which(lam2s$error == min(lam2s$error))
  lam1_best <- lam2s$best.lam1[ind_min_cve][1]
  lam2_best <- lam2s$lam2[ind_min_cve][1]
    
  # With the selected lambda_1 and lambda_2 we fit the MTL L21 model;
  ### Keep the result when lam2 equals 0. Also put lam2 = 0 into the loop.

  mtl_model <- MTC_L21(dat_mat, class_vec, 
                  lam1 = lam1_best,
                  lam2 = lam2_best)
    
    # We predict on the test set and output the measurements of learning performance
    pred <- predict(mtl_model, testing_set)
    
    confusion.mat <- confusionMatrix(data = pred, reference = truth)
    
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

print(learning.results, n = length(studies))

write.table(learning.results, file = "/nas/longleaf/home/tianyi96/TSPs_svm_results.csv")

save(x = learning.results, file = "TSPs_svm_results.Rdata")
save.image()
unlink("TSPs_svm_results.Rdata")
