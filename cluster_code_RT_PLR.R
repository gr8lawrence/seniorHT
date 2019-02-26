library(tidyverse)
library(caret)
library(parallel)
library(ncvreg) 
set.seed(100)

### We run our data with the penalized logistic regression model

dataPath <- "/Users/gr8lawrence/Desktop/Senior Honors Thesis/datasets/"
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
df.length <- num.common.genes/5 # Keep 1/5 of the total number of genes in the training sets for RMTL

# Append the classification to the training data
studies <- list()
for (i in 1:num.studies) {
  studies[[i]] <- add_column(expression.df[[i]], class = ranked.studies.df[[i]]$sampInfo$cluster.MT[!is.na(ranked.studies.df[[i]]$sampInfo$cluster.MT)]) 
}

# We create an empty tibble first to hold the leanring results
learning.results <- tibble(training.set = rep('NA', length(studies)^2), 
                           testing.set = rep(0, length(studies)^2),
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))

# We create a list to hold our final models
final_models <- list()

## Below is our loop for multi-task learning (MTL)
# Reduce the runtime for cross-validation by using paralleing
no_cores <- 4
cl <- makeCluster(no_cores)

# We create an empty tibble first to hold the leanring results
learning.results <- tibble(training.set = rep('NA', length(studies)^2), 
                           testing.set = rep('NA', length(studies)^2), 
                           strategy = rep("Single",  length(studies)^2),
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))
len <- length(studies)
for (i in 1:len) {
  studies.min.1 <- studies[-i]
  studies.names.min.1 <- studies.names[-i]
  validation.set <- studies[[i]][, 1:df.length]
  truth <- studies[[i]]$class # The truth vector
  
  ## Run the penalized logistic regression on one dataset, and validate on studies[i]
  for (j in 1:length(studies.min.1)) {
    learning.set <- studies.min.1[[j]]
    
    # Cross-validation to find the best penalty parameter
    x <- learning.set[ ,1:df.length]
    y <- learning.set$class
    cvfit <- ncvreg::cv.ncvreg(X = x, y = y, family = "binomial", 
                               cluster = cl, nfolds = length(y), seed = 100, max.iter = 10000) 
    
    # Fit the regression model based on the best lambda
    fit <- ncvreg(x, y, family = "binomial", lambda = cvfit$lambda.min, seed = 100, max.iter = 10000) 
    pred <- predict(fit, X = data.matrix(validation.set))
    confusion.mat <- confusionMatrix(data = as.factor(ifelse(pred < 0, "basal", "classical")),
                                     reference = truth)
    accu <- confusion.mat$overall[["Accuracy"]]
    sen <- confusion.mat$byClass[["Sensitivity"]]
    spe <- confusion.mat$byClass[["Specificity"]]
    learning.results[5 * (i - 1) + j,] <- c(studies.names.min.1[j],
                                            studies.names[i],
                                            signif(accu, digits = 4),
                                            signif(sen, digits = 4),
                                            signif(spe, digits = 4))
  }
  
  ## Run the logistic regression model on combined datasets
  learning.set <- studies.min.1[[1]]
  for (k in 1:(length(studies.min.1) - 1)) {
    learning.set <- rbind(learning.set, studies.min.1[[1 + k]])
  }
  
  x <- learning.set[ ,1:df.length]
  y <- learning.set$class
  cvfit <- ncvreg::cv.ncvreg(X = x, y = y, family = "binomial", 
                             cluster = cl, nfolds = length(y), seed = 100, max.iter = 10000) 
  # Fit the regression model based on the best lambda
  fit <- ncvreg(x, y, family = "binomial", lambda = cvfit$lambda.min, seed = 100, max.iter = 10000) 
  pred <- predict(fit, X = data.matrix(validation.set))
  confusion.mat <- confusionMatrix(data = as.factor(ifelse(pred < 0, "basal", "classical")),
                                   reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[5 * i,] <- c(paste("Combined without ", studies.names.min.1[i]),
                                studies.names[i],
                                "Combined",
                                signif(accu, digits = 4),
                                signif(sen, digits = 4),
                                signif(spe, digits = 4))
}
learning.results$strategy <- factor(learning.results$strategy, levels = c("Single", "Combined"))
learning.results <- add_column(learning.results, transformation = "Rank")
learning.results$transformation <- factor(learning.results$transformation, levels = c("Rank", "TSPs"))

print(learning.results, n = 2*length(studies)^2)

# We create a list to hold all models and learning results
final_results <- list("models" = final_models, "results" = learning.results)
print(learning.results, n = length(studies)^2)
save(x = final_results, file = "/Users/gr8lawrence/Desktop/Senior Honors Thesis/PLR_results/PLR_RT_results.Rdata")
save.image()

