## Functions for data preprocessing
## Written by: Tianyi Liu  

# Sort gene in alphanumeric order
geneSort <- function(d.set) 
{
  index <- order(d.set$featInfo$SYMBOL)
  d.set$ex[,] <- d.set$ex[index,]
  d.set$featInfo[,] <- d.set$featInfo[index,]
  
  return(d.set)
}

# Obtain the list of gene symbols from datasets
getGeneSymbols <- function(d.set) 
{
  symbols <- d.set$featInfo$SYMBOL
  
  return(symbols)
} 

# Subset the datasets to ones of common genes
matchCommonGenes <- function(d) 
{
  index <- match(commonGeneNames, d$featInfo$SYMBOL)
  d$ex <- d$ex[index, ]
  d$featInfo <- d$featInfo[index, ]
  
  return(d)
}
