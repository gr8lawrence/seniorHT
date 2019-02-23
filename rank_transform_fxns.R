## Functions for performing rank transformation
## Written by: Tianyi Liu  

# Rank-transform the count matrix
rankTransform <- function(d.set) 
{
  apply(d.set$ex, 2, rank)
  
  return(d.set)
}

# Extract the rank-transformed count matrix, remove the un-labeled observations, and append the labels
# Output is a data frame
extractData <- function(d.set, common.gene.names) 
{
  df <- as_tibble(t(d.set$ex))
  # Correct gene name formats (due to syntactic requirments by R)
  colnames(df) <- make.names(common.gene.names) 
  # Remove observations with unknown subtype
  df <- df[!is.na(d.set$sampInfo$cluster.MT),] 
  
  return(add_column(df, class = d.set$sampInfo$cluster.MT[!is.na(d.set$sampInfo$cluster.MT)]))
}


# Re-order the extracted data frame (from the function above) by their significance of differential expression
reorderGeneBySignificance <- function(ranked.exp.df.list, common.gene.names)
{
  num.common.genes <- length(common.gene.names)
  # Separate observations by subtypes
  expression.basal.df.list <- lapply(ranked.exp.df.list, function(x) x[x$class == "basal", ] )
  expression.classical.df.list <- lapply(ranked.exp.df.list, function(x) x[x$class == "classical", ] )
  # Transform dataset to matrices (Wilcoxon rank sum test)
  basal.matrices <- list()
  classical.matrices <- list()
  for (i in 1:num.studies) {
    basal.matrices[[i]] <- data.matrix(expression.basal.df.list[[i]][, 1:num.common.genes])
    classical.matrices[[i]] <- data.matrix(expression.classical.df.list[[i]][, 1:num.common.genes])
  }
  # Build up empty list to contain the p-values of single tests (in one dataset) for each gene
  gene.pvals <- list()
  for (i in 1:num.studies){
    gene.pvals[[i]] <- tibble(name = make.names(common.gene.names),
                            pval = rep(0, num.common.genes))
  }
  # Conduct Wilcoxon Rank-sum test for each genes in each estudy
  for (i in 1:num.studies) {
    for (j in 1:num.common.genes) {
      basal.vector <- basal.matrices[[i]][, j]
      classical.vector <- classical.matrices[[i]][, j]
      gene.pvals[[i]]$pval[j] <- wilcox.test(basal.vector, classical.vector, exact = FALSE)[["p.value"]]
    }
  }
  # Compute each gene's sum of log-10 p-values across studies
  gene.sum.log.pvals <- tibble(name = make.names(common.gene.names),
                               pval = rep(0, length(common.gene.names)))
  for (i in 1:num.common.genes) {
    gene.pvals.vector <- vector()
    for (j in 1:num.studies) {
      gene.pvals.vector[j] = gene.pvals[[j]]$pval[i]
    }
    log.10.gene.pvals <- log10(gene.pvals.vector)
    gene.sum.log.pvals$pval[i] <- sum(log.10.gene.pvals)
  }
  # Give the sorted indices of genes
  index <- order(gene.sum.log.pvals$pval) 
  # Re-order genes based on sorted indices
  reorderGene <- function(df) 
  {
    df.colnames <- colnames(df)
    df[ ,1:length(common.gene.names)] <- data.matrix(df[, index])
    colnames(df) <- df.colnames[index]
    
    return(as_tibble(df))
  }
  # Re-order genes
  ranked.exp.df.list <- lapply(ranked.exp.df.list, reorderGene)
  
  return(ranked.exp.df.list)
}

