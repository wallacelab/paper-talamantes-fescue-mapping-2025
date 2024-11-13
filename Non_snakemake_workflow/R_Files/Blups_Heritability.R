# The objective of this script is to turn my phenotypes into BLUPs and analyze the 
# genotypes aswell to hopefully get sometype of heritability caluculation out of it. 

# Installing libraries
library(lme4)
library(vcfR)
library(sommer)
library(dplyr)
library(tidyverse)
library(AGHmatrix)

# Loading data
residual_avg_24_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_data_avg.txt"
residual_avg_24 <- read.table(residual_avg_24_loc, header = TRUE, strip.white=TRUE)

parent_data_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Meta_Data/Mother_Father_Data.csv"
parent_data <- read.table(parent_data_loc, header = TRUE, strip.white=TRUE, sep = ",")
parent_data$ID <- gsub("-", "_", parent_data$ID)

# Load the VCF file
vcf_data_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf"
vcf_data <- read.vcfR(vcf_data_loc)
geno_matrix <- extract.gt(vcf_data, element = "GT")  
geno_matrix <- t(geno_matrix)

# Labeling the crosses and the parents.
residual_avg_24 <- merge(residual_avg_24,parent_data, by = "ID",all.x = TRUE)
residual_avg_24$Cross <- apply(residual_avg_24[, c("Mother", "Father")], 1, function(x) {
  if (all(x == "" | is.na(x))) {
    "parent"
  } else {
    paste(sort(x), collapse = "x")
  }
})

#Split the data into all, 314x312, and 314x310
residuals_310x314 <- subset(residual_avg_24, Cross == "310x314" | ID == "310" | ID == "314")
residuals_312x314 <- subset(residual_avg_24, Cross == "312x314" | ID == "312" | ID == "314")
residuals_all <- residual_avg_24

####################################
# FUnction!
####################################

###### Fix genotype files to work in R ###### 
convert_genotypes <- function(geno) {
  geno <- gsub("0/0", "0", geno)
  geno <- gsub("0/1", "1", geno)
  geno <- gsub("1/1", "2", geno)
  return(as.numeric(geno))
}

###### Getting Heritabilities function ###### 
find_heritability <- function(pheno_data, geno_matrix, trait) {
  
  # Ensure both datasets have the same IDs
  common_IDs <- intersect(pheno_data$ID, rownames(geno_matrix))
  
  # Subset data to include only common IDs
  pheno_data <- pheno_data[pheno_data$ID %in% common_IDs, ]
  geno_matrix <- geno_matrix[common_IDs, ]
  
  # Check if row names match
  if (!all(rownames(geno_matrix) == pheno_data$ID)) {
    stop("Row names of geno_data do not match the ID in pheno_data.")
  }
  
  # Convert geno_matrix to data frame, and then convert genotypes
  geno_data <- as.data.frame(geno_matrix)
  geno_data <- geno_data %>% mutate_all(convert_genotypes)
  
  # Convert geno_data back to matrix for kinship matrix calculation
  geno_data <- as.matrix(geno_data)
  
  # Create kinship matrix using Gmatrix function
  kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")
  
  # Ensure the data is in correct format (pheno_data should have necessary columns)
  pheno_data$ID <- as.character(pheno_data$ID)
  
  # Example mixed model for heritability (adjust trait to be a column name)
  model <- mmer(fixed = as.formula(paste(trait, "~ 1")),
                random = ~ vsr(ID, Gu = kinship_matrix),
                rcov = ~ units,
                data = pheno_data)
  
  # Summarize model results
  model_summary <- summary(model)
  varcomp_summary <- model_summary$varcomp
  
  # Extract variance components for heritability calculation
  var_ID <- varcomp_summary$VarComp[1]  # Genetic variance
  var_residual <- varcomp_summary$VarComp[2]  # Residual variance
  
  # Calculate heritability
  heritability <- var_ID / (var_ID + var_residual)
  
  # Return heritability estimate
  return(heritability)
}

#####################################
# Calculating Heritability
#####################################

find_heritability(residuals_all,geno_matrix,"Alkaloids_Res_avg")














# Ensuring the genomatrix and the residuals have the same IDs
common_IDs <- intersect(residuals_all$ID, rownames(geno_matrix))
residuals_all <- residuals_all[residuals_all$ID %in% common_IDs, ]
geno_matrix <- geno_matrix[common_IDs, ]
identical(rownames(geno_matrix), residuals_all$ID)



geno_data <- as.data.frame(geno_matrix)
geno_data <- geno_data %>%
  mutate_all(convert_genotypes)
geno_data <- as.matrix(geno_data)
kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")

if (!all(rownames(geno_data) == residuals_all$ID)) {
  stop("Row names of geno_data do not match the ID in residuals_all.")
}
residuals_all$ID <- as.character(residuals_all$ID)

# Example mixed model for heritability
model <- mmer(fixed = Alkaloids_Res_avg ~ 1,
              random = ~ vsr(ID, Gu = kinship_matrix),
              rcov = ~ units,
              data = residuals_all)

summary(model)

# Extract heritability estimate
varcomp_summary <- summary(model)$varcomp
var_ID <- varcomp_summary$VarComp[1]  # Genetic variance
var_residual <- varcomp_summary$VarComp[2]  # Residual variance
heritability <- var_ID / (var_ID + var_residual)
heritability




























# Creating models of my phenotype
model_biomass <- lmer(Delta_CT_OG ~ Data_Set + Extractor + Extraction_Date + (1|Cross), data = pheno_24)
model_alkaloid <- lmer(ng.g ~ Alkaloid_Plate + (1|Cross), data = pheno_24)


# Extract BLUPs for Biomass
blup_biomass_OG <- ranef(model_biomass)$Cross
colnames(blup_biomass_OG)[1] <- "BiomassOG_Blup"

# Extract BLUPs for ng alkaloid/g plant
blup_alkaloid <- ranef(model_alkaloid)$Cross
colnames(blup_alkaloid)[1] <- "Alkaloid_Blup"
blup_Data <- merge(blup_biomass_OG, blup_alkaloid, by = "row.names")
colnames(blup_Data)[1] <- "ID"

# Save the Blup data to get it ready for tassel
write.table(blup_Data, file = "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Blups.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Calculating blups for my data does not really work. You kinda need multiple samples to calculate Blups














