library(tidyverse)
library(caret)
library(parallel)
library(ranger) # fast implementation of random forest
library(e1071) # svm
set.seed(100)

# Initializing paralleling
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)

### this file processes the rna-seq data from 
  # 5 clinical trials

# Load data
dataPath <- "/Users/gr8lawrence/Desktop/Senior Honors Thesis/datasets/"
load(paste(dataPath, "Aguirre_seq.RData", sep = ""))
load(paste(dataPath, "COMPASS.RData", sep = ""))
load(paste(dataPath, "Linehan_Seq.RData", sep = ""))
load(paste(dataPath, "Moffitt_GEO_array.RData", sep = ""))
load(paste(dataPath, "TCGA_PAAD.RData", sep = ""))

# Reduce datasets to common genes
# adjust Aguirre's gene symbol to factors
# Aguirre_seq$featInfo$SYMBOL <- as.factor(Aguirre_seq$featInfo$SYMBOL)
# note: Aguirre is dropped for the moment

# Sort expression by gene name
geneSort <- function(d) { # d is the name of the dataset
  index <- order(d$featInfo$SYMBOL)
  d$ex[ , ] <- d$ex[index, ] 
  d$featInfo[ , ] <- d$featInfo[index, ]
  return(d)
}

COMPASS <- geneSort(COMPASS)
Linehan_Seq <- geneSort(Linehan_Seq)
Moffitt_GEO_array <- geneSort(Moffitt_GEO_array)
TCGA_PAAD <- geneSort(TCGA_PAAD)

# Get the list of gene symbols from datasets
getGeneSymbols <- function(d) {
  symbols <- d$featInfo$SYMBOL
  return(symbols)
}

commonNames1 <- intersect(getGeneSymbols(Linehan_Seq), getGeneSymbols(COMPASS))
commonNames2 <- intersect(getGeneSymbols(Moffitt_GEO_array), getGeneSymbols(TCGA_PAAD))
commonNames <- intersect(commonNames1, commonNames2)

# Subset the datasets to ones of common genes
matchCommonGenes <- function(d) {
  index <- match(commonNames, d$featInfo$SYMBOL)
  d$ex <- d$ex[index, ]
  d$featInfo <- d$featInfo[index, ]
  return(d)
}

Lin <- matchCommonGenes(Linehan_Seq) 
Com <- matchCommonGenes(COMPASS)
Mof <- matchCommonGenes(Moffitt_GEO_array)
Tcga <- matchCommonGenes(TCGA_PAAD)

# Rank transform the common genes for each sample

rankTransform <- function(d) {
  for (i in 1:dim(d$ex)[2]) {
    d$ex[,i] <- rank(d$ex[,i])
  }
  return(d)
}

rankedLin <- rankTransform(Lin)
rankedCom <- rankTransform(Com)
rankedMof <- rankTransform(Mof)
rankedTcga <- rankTransform(Tcga)

# Extract common genes here
geneNames <- rankedCom$featInfo$SYMBOL

# Use Wilcox rank sum test to select the highly expressed genes across datasets

wilcoxTest <- function(d) {
  df <- tibble(genename = geneNames, 
               pval = rep(0, dim(d$ex)[1]) ) # initiate a data frame
  for (i in 1:dim(d$ex)[1]) {
    pval <- as.double(wilcox.test(data.matrix(d$ex[i,]))[["p.value"]])
    df$pval[i] <- pval
  }
  return(df)
}

# Run the Random Forest model on the common genes
# Extract the data we are going to work on and tidy it into a single tibble

extractData <- function(d) {
  df <- as_tibble(t(d$ex))
  colnames(df) <- make.names(geneNames)
  return(add_column(df, category = d$sampInfo$cluster.MT))
}

#dfLin <- extractData(rankedLin)
dfCom <- extractData(rankedCom)
dfCom <- dfCom[!is.na(dfCom$category), ]
dfMof <- extractData(rankedMof)
dfMof <- dfMof[!is.na(dfMof$category), ]
dfTcga <- extractData(rankedTcga)
dfTcga <- dfTcga[!is.na(dfTcga$category), ]

# Combined datasets
dfComMof <- rbind(dfCom, dfMof)
dfComTcga <- rbind(dfCom, dfTcga)
dfMofTcga <- rbind(dfMof, dfTcga)

# The reference values of cancer type (truth) from three datasets
ref.Com <- rankedCom$sampInfo$cluster.MT[!is.na(rankedCom$sampInfo$cluster.MT)]
ref.Mof <- rankedMof$sampInfo$cluster.MT[!is.na(rankedMof$sampInfo$cluster.MT)]
ref.Tcga <- rankedTcga$sampInfo$cluster.MT[!is.na(rankedTcga$sampInfo$cluster.MT)]

# Run the random forest model based on a single training dataset
Com.rf <- ranger::ranger(category~., data = dfCom, mtry = 1500, num.trees = 100)
Mof.rf <- ranger::ranger(category~., data = dfMof, mtry = 1500, num.trees = 100)
Tcga.rf <- ranger::ranger(category~., data = dfTcga, mtry = 1500, num.trees = 100)

# prediciton using the random forest models trained on a single dataset
pred.rf.Com.Mof <- predict(Com.rf, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Com.Mof$predictions, reference = ref.Mof) # sen. = 1, spe. = 0.16

pred.rf.Com.Tcga <- predict(Com.rf, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Com.Tcga$predictions, reference = ref.Tcga) # sen. = 0.35, spe. = 0.99

pred.rf.Mof.Com <- predict(Mof.rf, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Mof.Com$predictions, reference = ref.Com) # sen. = 0.93, spe. = 0.33

pred.rf.Mof.Tcga <- predict(Mof.rf, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Mof.Tcga$predictions, reference = ref.Tcga) # sen. = 0.92, spe. = 0.25

pred.rf.Tcga.Com <- predict(Tcga.rf, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Tcga.Com$predictions, reference = ref.Com) # sen. = 0.93, spe. = 0.97

pred.rf.Tcga.Mof <- predict(Tcga.rf, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.Tcga.Mof$predictions, reference = ref.Mof) # sen. = 0.55, spe. = 0.82

# Run RF on combined datasets
ComMof.rf <- ranger::ranger(category~., data = dfComMof, mtry = 1500, num.trees = 100)
ComTcga.rf <- ranger::ranger(category~., data = dfComTcga, mtry = 1500, num.trees = 100)
MofTcga.rf <- ranger::ranger(category~., data = dfMofTcga, mtry = 1500, num.trees = 100)

# Prediction on the other datasets
pred.rf.ComMof.Tcga <- predict(ComMof.rf, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.ComMof.Tcga$predictions, reference = ref.Tcga) # sen. = 0.52, spe. = 0.99

pred.rf.ComTcga.Mof <- predict(ComTcga.rf, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.ComTcga.Mof$predictions, reference = ref.Mof) # sen. = 0.87, spe. = 0.45

pred.rf.MofTcga.Com <- predict(MofTcga.rf, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.rf.MofTcga.Com$predictions, reference = ref.Com) # sen. = 1, spe. = 0.81

# Run the supporting vector machine on a single training dataset
training.length <- length(geneNames)
Com.svm <- svm(category~., data = dfCom)
Mof.svm <- svm(category~., data = dfMof)
Tcga.svm <- svm(category~., data = dfTcga)

# Prediction on the other datasets using svm trained on Com
pred.Com.Com <- predict(Com.svm, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.Com.Com, reference = ref.Com) # sen. = 1, spe. = 1

pred.Com.Mof <- predict(Com.svm, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.Com.Mof, reference = ref.Mof) # sen. = 0, spe. = 1

pred.Com.Tcga <- predict(Com.svm, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.Com.Tcga, reference = ref.Tcga) # seb. = 0, spe. = 1

# ... svm trained on Mof
pred.Mof.Mof <- predict(Mof.svm, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.Mof.Mof, reference = ref.Mof)  # sen. = 0.944, spe. = 1

pred.Mof.Com <- predict(Mof.svm, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.Mof.Com, reference = ref.Com)  # sen. = 1, spe. = 0

pred.Mof.Tcga <- predict(Mof.svm, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.Mof.Tcga, reference = ref.Tcga)  # sen. = 1, spe. = 0

# ... svm trained on Tcga
pred.Tcga.Tcga <- predict(Tcga.svm, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.Tcga.Tcga, reference = ref.Tcga) # sen. = 1, spe. = 1

pred.Tcga.Com <- predict(Tcga.svm, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.Tcga.Com, reference = ref.Com) # sen. = 0, spe. = 1

pred.Tcga.Mof <- predict(Tcga.svm, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.Tcga.Mof, reference = ref.Mof) # sen. = 0, spe. = 1

# Run the svm on two combined datasets
ComMof.svm <- svm(category~., data = dfComMof)
ComTcga.svm <- svm(category~., data = dfComTcga)
MofTcga.svm <- svm(category~., data = dfMofTcga)

# svm confusion matrix
pred.ComMof.Tcga <- predict(ComMof.svm, dfTcga[ , 1:length(geneNames)])
confusionMatrix(data = pred.ComMof.Tcga, reference = ref.Tcga) # sen. = 1, spe. = 0

pred.ComTcga.Mof <- predict(ComTcga.svm, dfMof[ , 1:length(geneNames)])
confusionMatrix(data = pred.ComTcga.Mof, reference = ref.Mof) # sen. = 0, spe. = 1

pred.MofTcga.Com <- predict(MofTcga.svm, dfCom[ , 1:length(geneNames)])
confusionMatrix(data = pred.MofTcga.Com, reference = ref.Com) # sen. = 1, spe. = 0.61

