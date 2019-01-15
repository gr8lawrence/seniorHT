library(tidyverse)
library(caret)
library(ranger) # fast implementation of random forest
set.seed(100)

dataPath <- 
  "/Users/gr8lawrence/Desktop/Senior Honors Thesis/datasets/"
  # "/nas/longleaf/home/tianyi96/datasets_used/" # change this line to your local dataset directory

studies.names <-
  c("Aguirre-Seq",
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
studies.df <-
  list(
    Aguirre_seq_plus,
    Linehan_Seq_plus,
    COMPASS.2017_plus,
    Moffitt_GEO_array_plus,
    TCGA_PAAD_plus
  )

num.studies <- length(studies.df) # Number of datasets

# Sort gene in alphanumeric order
geneSort <- function(d) {
  # d is a dataset
  index <- order(d$featInfo$SYMBOL)
  d$ex[,] <- d$ex[index,]
  d$featInfo[,] <- d$featInfo[index,]
  return(d)
}

for (i in 1:num.studies) {
  studies.df[[i]] <- geneSort(studies.df[[i]])
}

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

for (i in 1:num.studies) {
  studies.df[[i]] <- matchCommonGenes(studies.df[[i]])
}

# Rank transform the common genes for each sample
rankTransform <- function(d) {
  for (i in 1:dim(d$ex)[2]) {
    d$ex[,i] <- rank(d$ex[,i])
  }
  return(d)
}

ranked.studies.df <- list()
for (i in 1:num.studies) {
  ranked.studies.df[[i]] <- rankTransform(studies.df[[i]])
}

# Extract ranked expression from each dataset, remove observations of unknown cancer classification, and write them into tibbles
extractData <- function(d) {
  df <- as_tibble(t(d$ex))
  colnames(df) <- make.names(commonGeneNames) # Get all gene names into correct name formats
  df <- df[!is.na(d$sampInfo$cluster.MT),] # Remove observations with unknown subtype
  return(add_column(df, class = d$sampInfo$cluster.MT[!is.na(d$sampInfo$cluster.MT)]))
}

expression.df <- list()
for (i in 1:num.studies) {
  expression.df[[i]] <- extractData(ranked.studies.df[[i]])
}

expression.basal.df <- list()
expression.classical.df <- list()
for (i in 1:num.studies) {
  expression.basal.df[[i]] <- expression.df[[i]][expression.df[[i]]$class == "basal", ]
  expression.classical.df[[i]] <- expression.df[[i]][expression.df[[i]]$class == "classical", ]
}

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

for (i in 1:num.studies) {
  expression.df[[i]] <- reorderGene(expression.df[[i]])  
}

df.length <- ceiling(length(commonGeneNames)/5)

learning.df <- list()
for (i in 1:num.studies) {
  learning.df[[i]] <- add_column(expression.df[[i]][ , 1:df.length], 
                                 class = ranked.studies.df[[i]]$sampInfo$cluster.MT[!is.na(ranked.studies.df[[i]]$sampInfo$cluster.MT)]) 
}

studies <- learning.df

# We create an empty tibble first to hold the leanring results
learning.results <- tibble(learning.set = rep('NA', length(studies)^2), 
                           validation.set = rep('NA', length(studies)^2), 
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))
len <- length(studies)

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
  pred <- predict(rf, validation.set)
  pred_column <- factor(ifelse(pred$predictions[,1] > 0.5, "basal", "classical"),
                            levels = c("basal", "classical"))
  confusion.mat <- confusionMatrix(data = , 
                                   reference = truth)
  accu <- confusion.mat$overall[["Accuracy"]]
  sen <- confusion.mat$byClass[["Sensitivity"]]
  spe <- confusion.mat$byClass[["Specificity"]]
  learning.results[10 * i,] <-
    c(
      paste("comb. minus", studies.names[i],
            sep = " "),
      studies.names[i],
      signif(accu, digits = 4),
      signif(sen, digits = 4),
      signif(spe, digits = 4)
    )

}

print(learning.results, n = length(studies)^2)

write.table(learning.results, file = "/nas/longleaf/home/tianyi96/learning_results_RT_rf.csv")

save(x = learning.results, file = "learning_result_RT_rf.Rdata")
save.image()
