# Objective: This script will make PCA graphs of 2024 phenotype 

# Import libraries 
library(tidyverse)


# Loading in the data
phenotypes_2024_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv"
phenotypes_2024 <- read.csv(phenotypes_2024_loc, header = TRUE, strip.white=TRUE)




# Finding the cross
phenotypes_2024$cross <- apply(phenotypes_2024[, c("Mother", "Father")], 1, function(x) {
  if (all(is.na(x)) || all(x == "")) {
    return("Parent")
  } else {
    return(paste(sort(x), collapse = "x"))
  }
})


colnames(phenotypes_2024)
phenotype_data <- subset(phenotypes_2024, select = c("Delta_CT_OG","ng.g"))


# Makes the data more normal 
pca_result <- prcomp(phenotype_data, scale. = TRUE)  # Scaling the data is important











