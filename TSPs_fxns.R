## Functions for making the TSPs
## Written by: Tianyi Liu  

# Customary function to transpose a tibble while preserves the names
transposeTibble <- function(tb, df.length) 
{
  data.mat <- t(data.matrix(tb))[1:df.length,]
  tr.tibble <- as.tibble(data.mat)
  colnames(tr.tibble) <- rownames(tb)
  rownames(tr.tibble) <- colnames(tb)[1:df.length]
  
  return(tr.tibble)
}

## This part is mostly inherited from Dr. Rashid's code; minor changes were made to meet the needs of this script
# train_sub is your expression matrix, with rows labeled with the gene names, columns labeled with the sample names
# TSPs is a two column matrix, where the rows are TSPs.  
# column 1 is the 1st gene name in the TSP, column 2 is the 2nd gene name in the TSP
ind_fun = function(train_sub, TSPs)
{
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

# Function to make a matrix for the TSPs: tr.tibble is a transposed tibble from above, gene.subset.names the names of the subset of genes
# used, and df.length the number of genes used to make the TSPs.
makeTSPs <- function(transposed.learning.df, gene.subset.names, df.length) 
{
  # Now make the TSPs for the gene pairs below
  TSPs.mat <- matrix(0, nrow = df.length*(df.length - 1)/2, ncol = 2)
  range <- df.length - 1
  range.total = 0
  while(range > 0) {
    # The index of the first gene's name in the pair in the subset of gene names
    index1 <- df.length - range 
    for (i in 1:range) {
      # Row index of the pair in the TSPs matrix
      index <- range.total + i 
      # The index of the second gene's name in the pair
      index2 <- df.length - range + i 
      TSPs.mat[index, ] <- c(gene.subset.names[index1], gene.subset.names[index2])
    }
    range.total = range.total + range 
    range = range - 1
  }
  # apply ind_fun to every member of the learning datasets
  studies <- lapply(transposed.learning.df, ind_fun, TSPs = TSPs.mat) 
  # Determine which TSPs does not vary across different subtypes (always 0 or 1 in all observations) as they do not provide useful information.
  all_1_inds <- lapply(studies, function(x){which(apply(x, 2, sum) == (dim(x)[1]))})
  all_0_inds <- lapply(studies, function(x){which(apply(x, 2, sum) == 0)})
  # Find out the unique indexes of those columns
  unique_inds <-  unique(c(unlist(all_1_inds), unlist(all_0_inds)))
  # remove those columns from the studies
  studies <- lapply(studies, function(x){x[ ,-unique_inds]})
  # scale the data to mean 0 and variance 1 using the Standardize function
    #studies <- lapply(studies, function(x) as.tibble(apply(x, 2, scale)))
  studies <- lapply(studies, function(x){as.tibble(apply(x, 2, function(y){(y - mean(y))/sd(y)}))})
  
  return(studies)
}