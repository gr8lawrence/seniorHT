## This file is for file loading
## Written by: Tianyi Liu


dataPath <- "datasets_used/"
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