library(tidyverse)
library(caret)
library(ranger) # package containing rf
set.seed(100)

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


# This variable indicates how many genes of top differential expression are kept to make TSPs
df.length <- 500 # We will only keep this number of genes, 
                # this will give us a maximum of (df.length)*(df.length - 1)/2 variables in TSPs

# Customary function to transpose a tibble while preserves the names
transpose_tibble <- function(tb) {
  data.mat <- t(data.matrix(tb))[1:df.length,]
  tr.tibble <- as.tibble(data.mat)
  colnames(tr.tibble) <- rownames(tb)
  rownames(tr.tibble) <- colnames(tb)[1:df.length]
  return(tr.tibble)
}

transposed.learning.df <- lapply(expression.df, transpose_tibble)

# Now make the TSPs for the gene pairs below

gene.subset.names <- colnames(expression.df[[1]])[1:df.length]
TSPs.mat <- matrix(0, nrow = df.length*(df.length - 1)/2, ncol = 2)

range <- df.length - 1
range.total = 0
while(range > 0) {
 index1 <- df.length - range # index of the first gene's name in the pair in the subset of gene names
 for (i in 1:range) {
    index <- range.total + i # row index of the pair in the TSPs matrix
    index2 <- df.length - range + i # index of the second gene's name in the pair
    TSPs.mat[index, ] <- c(gene.subset.names[index1], gene.subset.names[index2])
  }
  range.total = range.total + range 
  range = range - 1
}

## This part is mostly inherited from Dr. Rashid's code; minor changes were made to meet the needs of this script
# train_sub is your expression matrix, with rows labeled with the gene names, columns labeled with the sample names
# TSPs is a two column matrix, where the rows are TSPs.  
# column 1 is the 1st gene name in the TSP, column 2 is the 2nd gene name in the TSP

ind_fun = function(train_sub, TSPs){
  indmat = matrix(-1, ncol(train_sub), nrow(TSPs))
  for(i in 1:nrow(TSPs)){
    p1 = which(rownames(train_sub) == TSPs[i,1])
    p2 = which(rownames(train_sub) == TSPs[i,2])
    indmat[,i] = (train_sub[p1,] > train_sub[p2,])^2
  }
  indmat = as.tibble(indmat)
  colnames(indmat) = make.names(apply(TSPs, 1, paste, collapse = "/"))
  return(indmat)
}

# apply ind_fun to every member of the learning datasets
studies <- lapply(transposed.learning.df, ind_fun, TSPs = TSPs.mat) 

# Now determine which TSPs does not vary across different subtypes (always 0 or 1 in all observations).
all_1_inds <- lapply(studies, function(x) which(apply(x, 2, sum) == (dim(x)[1])))
all_0_inds <- lapply(studies, function(x) which(apply(x, 2, sum) == 0))

# Find out the unique indexes of those columns
unique_inds <-  unique(c(unlist(all_1_inds), unlist(all_0_inds)))

# remove those columns from the studies
studies <- lapply(studies, function(x) x[ ,-unique_inds])

# scale the data to mean 0 and variance 1 using the Standardize function
#studies <- lapply(studies, function(x) as.tibble(apply(x, 2, scale)))
studies <- lapply(studies, function(x) as.tibble(apply(x, 2, function (y) (y - mean(y))/sd(y))))

# This is the same as the total number of TSPs in the final training data
new.df.length <- dim(studies[[1]])[2]

# Append the classification to the training data
for (i in 1:num.studies) {
  studies[[i]] <- add_column(studies[[i]], class = ranked.studies.df[[i]]$sampInfo$cluster.MT[!is.na(ranked.studies.df[[i]]$sampInfo$cluster.MT)]) 
}

# We create an empty tibble first to hold the leanring results
learning.results <- tibble(training.set = rep('NA', length(studies)^2), 
                           testing.set = rep('NA', length(studies)^2), 
                           accuracy = rep(0, length(studies)^2), 
                           sensitivity = rep(0, length(studies)^2), 
                           specificity = rep(0, length(studies)^2))

# We also create a list to store the final learning models
final_models <- list()

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
  
  rf <- ranger::ranger(class~., data = training_set, probability = TRUE, memory.saving = TRUE)
  final_models[[5 * i]] <- rf
  pred <- predict(rf, testing_set)
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

print(learning.results, n = 2*length(studies)^2)

write.table(learning.results, file = "/nas/longleaf/home/tianyi96/TSPs_svm_results.csv")

save(x = final_results, file = "TSPs_rf_results.Rdata")
save.image()

