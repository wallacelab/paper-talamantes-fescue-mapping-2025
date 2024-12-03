# Objective: Make phenotype data into residual data and then subset all data 
# into all three data sets, 23, 24, avg all

library(tidyverse)
library(multcompView)
library(gvlma)
library(ggpubr)
library(car)
library(reshape2)
library(gridExtra)
library(grid)  
library(vcfR)
library(AGHmatrix)
library(sommer)
library(data.table)


##############################################
# Step 1: Making every data set
##############################################
# Loading in the data
# The two data sets you need are all_Data_2024 and phenotype_Data_2023
phenotype_Data <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/All_Data_Filtered/phenotype_data.txt", header = TRUE)
phenotype_Data$Delta_CT[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "314" & phenotype_Data$Data_Set != "315x320"), ]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Data_Set == "315x320"),]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Extraction_Date == "03/16/23"),]

all_Data_2024 <- read.csv("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv", header = TRUE)

# Fixing naming conventions
phenotype_Data <- phenotype_Data %>%
  rename(Delta_CT_adj = Delta_CT)
phenotype_Data_2023 <- phenotype_Data
all_Data_2024 <- all_Data_2024 %>%
  rename(ID = Treatment)



########## Getting residual data set for 2023 and 2024 avraged #################

# Removing data thats not star cross
list_314x310 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/310x314_list.txt")  
list_314x312 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/312x314_list.txt")
new_rows <- data.frame(V1 = c(301, 302, 303, 304, 305, 306, 307, 308, 310, 312, 313, 314, 315, 316, 318, 319, 320))

#List of non star cross parents
no_star_parents <- tibble(V1 = c(301, 302, 303, 304, 305, 306, 307, 308, 313, 315, 316, 318, 319, 320))
list_star <- rbind(list_314x310, list_314x312, new_rows)
list_star$V1 <- sub("_dupped.bam", "", list_star$V1)
phenotype_Data <- phenotype_Data[phenotype_Data$ID %in% list_star$V1, ]

# Subsetting the data into separate chunks
allpehnotype_data_export_23 <- subset(phenotype_Data, select = c(ID,Delta_CT_adj,Delta_CT_OG,ng.g,Data_Set))
allpehnotype_data_export_23$Year <- "2023"
allpehnotype_data_export_24 <- subset(all_Data_2024, select = c(ID,Delta_CT_adj,Delta_CT_OG,ng.g,Data_Set))
allpehnotype_data_export_24$Year <- "2024"

head(allpehnotype_data_export_23)
head(allpehnotype_data_export_24)

#Recombining the data into one large table
phenotypes23_24 <- rbind(allpehnotype_data_export_23,allpehnotype_data_export_24)
phenotypes23_24$ID <- gsub("-", "_", phenotypes23_24$ID)
head(phenotypes23_24,15)

#Removing batch effects and leaving only residuals.
lm_model_alk <- lm(ng.g ~ Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Alkaloids_Res <- resid(lm_model_alk)

lm_model_CT_OG <- lm(Delta_CT_OG ~ Data_Set + Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Delta_CT_OG_Res <- resid(lm_model_CT_OG)

lm_model_CT_adj <- lm(Delta_CT_adj ~ Data_Set + Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Delta_CT_adj_Res <- resid(lm_model_CT_adj)

head(phenotypes23_24,15)

# Now that we have residuals we have to get avarages and export the data
Alkaloid_residuals_avaraged <- phenotypes23_24 %>%
  group_by(ID) %>%
  summarise(Alkaloids_Res_avg = mean(Alkaloids_Res, na.rm = TRUE)) %>%
  filter(!is.na(Alkaloids_Res_avg))

DeltaCT_adj_residuals_avaraged <- phenotypes23_24 %>%
  group_by(ID) %>%
  summarise(DeltaCT_adj_Res_avg = mean(Delta_CT_adj_Res, na.rm = TRUE)) %>%
  filter(!is.na(DeltaCT_adj_Res_avg))

DeltaCT_OG_residuals_avaraged <- phenotypes23_24 %>%
  group_by(ID) %>%
  summarise(DeltaCT_OG_Res_avg = mean(Delta_CT_OG_Res, na.rm = TRUE)) %>%
  filter(!is.na(DeltaCT_OG_Res_avg))

Residual_data_avg <- merge(Alkaloid_residuals_avaraged, 
                           DeltaCT_adj_residuals_avaraged, 
                           by = "ID")

Residual_data_avg <- merge(Residual_data_avg, 
                           DeltaCT_OG_residuals_avaraged, 
                           by = "ID")


# This is the 2023 and 2024 residual data avraged with outliars
write.table(Residual_data_avg, file = '/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_data_avg.txt', sep = '\t', row.names=FALSE)


################### Getting residual data for 2023 only ########################
phenotype_Data_23 <- phenotype_Data
phenotype_Data_23 <- subset(phenotype_Data_23, select = c(ID,Delta_CT_adj,Delta_CT_OG,ng.g,Data_Set))

# Removing batch effects and recombining the data 
lm_model_alk <- lm(ng.g ~ 1, data = phenotype_Data_23, na.action = na.exclude)
phenotype_Data_23$Alkaloids_Res <- resid(lm_model_alk)

lm_model_CT_OG <- lm(Delta_CT_OG ~ Data_Set, data = phenotype_Data_23, na.action = na.exclude)
phenotype_Data_23$Delta_CT_OG_Res <- resid(lm_model_CT_OG)

lm_model_CT_adj <- lm(Delta_CT_adj ~ Data_Set, data = phenotype_Data_23, na.action = na.exclude)
phenotype_Data_23$Delta_CT_adj_Res <- resid(lm_model_CT_adj)

head(phenotype_Data_23,15)
phenotype_Data_23 <- subset(phenotype_Data_23, select = c(ID,Alkaloids_Res,Delta_CT_OG_Res,Delta_CT_adj_Res))

# This is the 2023 Data with no out liar removal
phenotype_Data_23$ID <- gsub("-", "_", phenotype_Data_23$ID) 
write.table(phenotype_Data_23, file = '/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2023_Data/Resisduals_Starcross_2023.txt', sep = '\t', row.names=FALSE)


################### Getting residual data for 2024 only ########################

# This should be the 2024 phenotype data
allpehnotype_data_export_24 <- subset(allpehnotype_data_export_24, select = -c(Year))
head(allpehnotype_data_export_24)

# Removing batch effects and recombining the data 
lm_model_alk <- lm(ng.g ~ 1, data = allpehnotype_data_export_24, na.action = na.exclude)
allpehnotype_data_export_24$Alkaloids_Res <- resid(lm_model_alk)

lm_model_CT_OG <- lm(Delta_CT_OG ~ Data_Set, data = allpehnotype_data_export_24, na.action = na.exclude)
allpehnotype_data_export_24$Delta_CT_OG_Res <- resid(lm_model_CT_OG)

lm_model_CT_adj <- lm(Delta_CT_adj ~ Data_Set, data = allpehnotype_data_export_24, na.action = na.exclude)
allpehnotype_data_export_24$Delta_CT_adj_Res <- resid(lm_model_CT_adj)

allpehnotype_data_export_24 <- subset(allpehnotype_data_export_24, select = c(ID,Alkaloids_Res,Delta_CT_OG_Res,Delta_CT_adj_Res))
head(allpehnotype_data_export_24,15)

# This is the 2024 Residual Data with no out liar removal
allpehnotype_data_export_24$ID <- gsub("-", "_", allpehnotype_data_export_24$ID) 
write.table(allpehnotype_data_export_24, file = '/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Resisduals_Starcross_2024.txt', sep = '\t', row.names=FALSE)

##################################################################
# Outliar exploration and removal
##################################################################

# load in previously saved data
Residual_data_avg <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_data_avg.txt", header = TRUE)
Residual_Data_23 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2023_Data/Resisduals_Starcross_2023.txt", header = TRUE)
Residual_Data_24 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Resisduals_Starcross_2024.txt", header = TRUE)

#Making all column names the same
colnames(Residual_data_avg) <- gsub("_avg", "", colnames(Residual_data_avg))
Residual_data_avg <- Residual_data_avg %>%
  rename(Delta_CT_adj_Res = DeltaCT_adj_Res)
Residual_data_avg <- Residual_data_avg %>%
  rename(Delta_CT_OG_Res = DeltaCT_OG_Res)

#### Make Function to display all phenotypes in graphical format ##########

Pheno_Graph <- function(Residual_data, title){
  if (!is.character(title)) {
    stop("The last three arguments must be strings.")
  }
  
  # Graphs to look at the avraged residual data
  p1 <- ggplot(Residual_data, aes(x = Alkaloids_Res)) + 
    geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Alkaloid Residuals", x = "Value", y = "Frequency") +
    theme_bw()
  
  p2 <- ggplot(Residual_data, aes(x = Delta_CT_adj_Res)) + 
    geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
    labs(title = "Delta CT Adj Residuals", x = "Value", y = "Frequency") +
    theme_bw()
  
  p3 <- ggplot(Residual_data, aes(x = Delta_CT_OG_Res)) + 
    geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
    labs(title = "Delta CT OG Residuals", x = "Value", y = "Frequency") +
    theme_bw()
  
  title1 <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"))
  residual_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2, top = title1)
  
  return(residual_plot)
}
############################## End Function ####################################
# Exploring residual data for outliars
Pheno_Graph(Residual_data_avg,"Residual Data of 2023 and 2024 Avaraged" )
Pheno_Graph(Residual_Data_23,"Residual Data of 2023")
Pheno_Graph(Residual_Data_24,"Residual Data of 2024")

#################Function to remove outliars ##################################

remove_outliers <- function(data, stds) {
  # Select columns 2, 3, and 4
  cols_to_check <- data[, 2:4]
  
  # Calculate the mean and standard deviation for each column
  means <- apply(cols_to_check, 2, mean, na.rm = TRUE)
  sds <- apply(cols_to_check, 2, sd, na.rm = TRUE)
  
  # Create a logical vector indicating rows within the limits
  within_limits <- apply(cols_to_check, 1, function(row) {
    # Check for NA in the row
    if (any(is.na(row))) {
      return(FALSE) # Exclude rows with NA
    }
    # Check if all columns are within the specified range
    all(abs(row - means) <= stds * sds)
  })
  
  # Subset data to include only rows that fall within the limits
  filtered_data <- data[within_limits, ]
  # Gets rid of duplicates
  filtered_data <- filtered_data[!duplicated(filtered_data$ID), ]
  
  return(filtered_data)
}

############################## End Function ####################################

Residual_data_avg_outliars_rm <- remove_outliers(Residual_data_avg, 2.5)
Pheno_Graph(Residual_data_avg,"Residual Data of 2023 and 2024 Avaraged" )
Pheno_Graph(Residual_data_avg_outliars_rm,"Residual Data of 2023 and 2024 Avaraged" )
nrow(Residual_data_avg) - nrow(Residual_data_avg_outliars_rm)

Residual_Data_23_outliars_rm <- remove_outliers(Residual_Data_23, 2.5)
Pheno_Graph(Residual_Data_23,"Residual Data of 2023")
Pheno_Graph(Residual_Data_23_outliars_rm,"Residual Data of 2023")
nrow(Residual_Data_23) - nrow(Residual_Data_23_outliars_rm)

Residual_Data_24_outliars_rm <- remove_outliers(Residual_Data_24, 2.5)
Pheno_Graph(Residual_Data_24,"Residual Data of 2024")
Pheno_Graph(Residual_Data_24_outliars_rm,"Residual Data of 2024")
nrow(Residual_Data_24) - nrow(Residual_Data_24_outliars_rm)

write.table(Residual_data_avg_outliars_rm,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt")
write.table(Residual_Data_23_outliars_rm,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt")
write.table(Residual_Data_24_outliars_rm,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt")


################################################################################
# Now seperating data to make all crosses.
################################################################################
#Reload Residual Data sets
Residual_data_avg_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt", header = TRUE)
Residual_Data_23_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt", header = TRUE)
Residual_Data_24_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt", header = TRUE)

# Loading in lists
list_314x310 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/310x314_list.txt")  
list_314x312 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Lists/Parental_Lists/312x314_list.txt")

list_314x310$V1 <- gsub("_dupped\\.bam", "", list_314x310$V1)  # Remove "_dupped.bam"
list_314x310$V1 <- gsub("-", "_", list_314x310$V1)  # Remove "_dupped.bam"
list_314x312$V1 <- gsub("_dupped\\.bam", "", list_314x312$V1)  # Remove "_dupped.bam"
list_314x312$V1 <- gsub("-", "_", list_314x312$V1)  # Remove "_dupped.bam"

########### Function that subsets residual data frame by list of IDs ###########
residual_subsetter <- function(residual,IDlist){
  common_IDs <- intersect(residual$ID, IDlist$V1)
  residual_filtered <- residual %>% filter(ID %in% common_IDs)
  return(residual_filtered)
}

Residual_data_avg_outliars_rm_314x310 <- residual_subsetter(Residual_data_avg_outliars_rm,list_314x310 )
Residual_data_avg_outliars_rm_314x312 <- residual_subsetter(Residual_data_avg_outliars_rm,list_314x312 )

Residual_Data_23_outliars_rm_314x310 <- residual_subsetter(Residual_Data_23_outliars_rm,list_314x310 )
Residual_Data_23_outliars_rm_314x312 <- residual_subsetter(Residual_Data_23_outliars_rm,list_314x312 )

Residual_Data_24_outliars_rm_314x310 <- residual_subsetter(Residual_Data_24_outliars_rm,list_314x310 )
Residual_Data_24_outliars_rm_314x312 <- residual_subsetter(Residual_Data_24_outliars_rm,list_314x312 )

# All datasets are now complete.
write.table(Residual_data_avg_outliars_rm_314x310,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x310.txt")
write.table(Residual_data_avg_outliars_rm_314x312,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x312.txt")

write.table(Residual_Data_23_outliars_rm_314x310,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x310.txt")
write.table(Residual_Data_23_outliars_rm_314x312,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x312.txt")

write.table(Residual_Data_24_outliars_rm_314x310,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x310.txt")
write.table(Residual_Data_24_outliars_rm_314x312,"/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x312.txt")

################################################################################
# Calculating heritability
################################################################################
# Load in all 9 datasets
Residual_data_avg_outliars_rm_314x310 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x310.txt", header = TRUE)
Residual_data_avg_outliars_rm_314x312 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x312.txt", header = TRUE)
Residual_data_avg_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt", header = TRUE)

Residual_Data_23_outliars_rm_314x310 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x310.txt", header = TRUE)
Residual_Data_23_outliars_rm_314x312 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x312.txt", header = TRUE)
Residual_Data_23_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt", header = TRUE)

Residual_Data_24_outliars_rm_314x310 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x310.txt", header = TRUE)
Residual_Data_24_outliars_rm_314x312 <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x312.txt", header = TRUE)
Residual_Data_24_outliars_rm <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt", header = TRUE)

# Load the VCF file
vcf_data_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf"
vcf_data <- read.vcfR(vcf_data_loc)
geno_matrix <- extract.gt(vcf_data, element = "GT")  
geno_matrix <- t(geno_matrix)

################ Functions needed for heritability finding ##################### 
##### Function to convert genotype into something usable ####
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
  heritability <- round(heritability, 4)
  resultstable <- data.table(
    H = heritability,
    GeneticVar = var_ID,
    ResidualVar = var_residual,
    N = model_summary$groups[1]
  )
  return(resultstable)
}

############################## End Function #################################### 
########################### 0 heritability fixer ###############################
# Creates a function that detects 0 heritability datasets, deletes more outliars
# Sends it back to recalculate heritability.

refine_heritability <- function(pheno_data, geno_matrix, trait, max_removals = 5, Hcap = .65, Hmin = .03) {
  removals <- 0
  removed_individuals <- character()  # Initialize an empty vector to store removed IDs
  
  while (removals < max_removals) {
    results <- find_heritability(pheno_data, geno_matrix, trait)
    heritability <- results$H[1]  # Extract heritability
    
    if (heritability >= Hmin && heritability <= Hcap) {
      cat("Heritability is ",heritability," which is acceptable. Returning results.\n")
      # Add removed individuals to the results
      results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
      return(results)
      
    }
    
    # Remove the furthest point from the mean
    cat("Heritability is ",heritability," which is unacceptable\n")
    mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    furthest_index <- which.max(abs_diff)
    
    cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    
    removals <- removals + 1
  }
  # Add removed individuals to the results
  cat("Heritability could not get into acceptable range, returning last results table\n")
  results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
  return(results)
}
########################## End FUnction ########################################
################# graphing heritability as individuals get removed #############
# The function will remove a single individual, graph the heritability on the y
# N on the x. itll do this until ...

graph_heritability <- function(pheno_data, geno_matrix, trait, max_removals = 14, Hdif = .04) {
  removals <- 0
  removed_individuals <- character()  # Initialize an empty vector to store removed IDs
  heritability_old = 100
  N_list <- c()
  H_list <- c()
  
  
  while (removals < max_removals) {
    
    results <- find_heritability(pheno_data, geno_matrix, trait)
    heritability <- results$H[1]  # Extract heritability
    heritability_dif <- abs(heritability - heritability_old)
    # Append the current N and H values
    N_list <- c(N_list, nrow(pheno_data))
    H_list <- c(H_list, heritability)
    cat("######################## \n")
    cat("N_list new value ", N_list, "\n")
    cat("H_list new value ", H_list, "\n")
    
    # if (heritability_dif <= Hdif && heritability != 0) {
    #   cat("Heritability difference is ",heritability_dif," which is acceptable. Returning results.\n")
    #   # Add removed individuals to the results
    #   results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
    #   
    #   # Remove and calculate 2 more before finishing.
    #   cat("Removing one more before finishing \n")
    #   mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    #   abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    #   furthest_index <- which.max(abs_diff)
    #   cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    #   removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    #   pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    #   heritability_old <- results$H[1]  
    #   removals <- removals + 1
    #   results <- find_heritability(pheno_data, geno_matrix, trait)
    #   heritability <- results$H[1]  # Extract heritability
    #   heritability_dif <- abs(heritability - heritability_old)
    #   # Append the current N and H values
    #   N_list <- c(N_list, nrow(pheno_data))
    #   H_list <- c(H_list, heritability)
    #   
    #   # Graphing the final product
    #   data <- data.frame(N = N_list, H = H_list)
    #   H_plot <- ggplot(data, aes(x = N, y = H)) +
    #     geom_point(color = "blue", size = 3) +  # Add points
    #     geom_line(color = "red") +             # Add a line connecting the points
    #     labs(title = "N vs Heritability",
    #          x = "N",
    #          y = "Heritability") +
    #     theme_bw() 
    #   
    #   return(list(graph = H_plot, Hdata = data, removed = removed_individuals))
    # 
    # }
    
    # Remove the furthest point from the mean
    cat("Heritability difference is ",heritability_dif," which is unacceptable\n")
    mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    furthest_index <- which.max(abs_diff)
    
    cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    heritability_old <- results$H[1]      # Update old heritability
    removals <- removals + 1
  }
  
  # Add removed individuals to the results
  # Graphing the final product
  data <- data.frame(N = N_list, H = H_list)
  H_plot <- ggplot(data, aes(x = N, y = H)) +
    geom_point(color = "blue", size = 3) +  # Add points
    geom_line(color = "red") +             # Add a line connecting the points
    labs(title = "N vs Heritability",
         x = "N",
         y = "Heritability") +
    theme_bw() 
  
  return(list(graph = H_plot, Hdata = data, removed = removed_individuals))
}
########################## End FUnction ########################################
# Calculating Heritability from all data sets

# # Both years avaraged, heritabilities 
# stats_avg_310_alk <- refine_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res")
# H_avg_310_alk <- stats_avg_310_alk$H[1]
# 
# stats_avg_312_alk <- refine_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res")
# H_avg_312_alk <- stats_avg_312_alk$H[1]
# 
# stats_avg_star_alk <- refine_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Alkaloids_Res")
# H_avg_star_alk <- stats_avg_star_alk$H[1]
# 
# stats_avg_310_ct <- refine_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_310_ct <- stats_avg_310_ct$H[1]
# 
# stats_avg_312_ct <- refine_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_312_ct <- stats_avg_312_ct$H[1]
# 
# stats_avg_star_ct <- refine_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_star_ct <- stats_avg_star_ct$H[1]
# 
# 
# # Just the 2023 heritabilities
# stats_2023_310_alk <- refine_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res")
# H_2023_310_alk <- stats_2023_310_alk$H[1]
# 
# stats_2023_312_alk <- refine_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res")
# H_2023_312_alk <- stats_2023_312_alk$H[1]
# 
# stats_2023_star_alk <- refine_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Alkaloids_Res")
# H_2023_star_alk <- stats_2023_star_alk$H[1]
# 
# stats_2023_310_ct <- refine_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res", 15, .6)
# H_2023_310_ct <- stats_2023_310_ct$H[1]
# 
# stats_2023_312_ct <- refine_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2023_312_ct <- stats_2023_312_ct$H[1]
# 
# stats_2023_star_ct <- refine_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res", 15, .60)
# H_2023_star_ct <- stats_2023_star_ct$H[1]
# 
# # Just the 2024 Heritabilities
# stats_2024_310_alk <- refine_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res", 15, .60)
# H_2024_310_alk <- stats_2024_310_alk$H[1]
# 
# stats_2024_312_alk <- refine_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res", 20, .60)
# H_2024_312_alk <- stats_2024_312_alk$H[1]
# 
# stats_2024_star_alk <- refine_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Alkaloids_Res", 20, .60)
# H_2024_star_alk <- stats_2024_star_alk$H[1]
# 
# stats_2024_310_ct <- refine_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res", 5, .60)
# H_2024_310_ct <- stats_2024_310_ct$H[1]
# 
# stats_2024_312_ct <- refine_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2024_312_ct <- stats_2024_312_ct$H[1]
# 
# stats_2024_star_ct <- refine_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2024_star_ct <- stats_2024_star_ct$H[1]


########## Plotting Heritability platoe ################

plot_avg_310_alk <- graph_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm_314x310)/2, .04)
plot_avg_312_alk <- graph_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm_314x312)/2, .04)
plot_avg_star_alk <- graph_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm)/2)
plot_avg_310_ct <- graph_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm_314x310)/2)
plot_avg_312_ct <- graph_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm_314x312)/2)
plot_avg_star_ct <- graph_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm)/2)

plot_2023_310_alk <- graph_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm_314x310)/2)
plot_2023_312_alk <- graph_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm_314x312)/2)
plot_2023_star_alk <- graph_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm)/2)
plot_2023_310_ct <-graph_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm_314x310)/2)
plot_2023_312_ct <- graph_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm_314x312)/2)
plot_2023_star_ct <- graph_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm)/2)

plot_2024_310_alk <- graph_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm_314x310)/2)
plot_2024_312_alk <- graph_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm_314x312)/2)
plot_2024_star_alk <- graph_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm)/2)
plot_2024_310_ct <- graph_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_24_outliars_rm_314x310)/2)
plot_2024_312_ct <- graph_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res", nrow(Residual_Data_24_outliars_rm_314x312)/2)
plot_2024_star_ct <- graph_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_24_outliars_rm)/2)




# Plot for avg data
# Create overall titles
title_avg <- textGrob("Average Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_avg_310_alk$graph, plot_avg_312_alk$graph, plot_avg_star_alk$graph,
              plot_avg_310_ct$graph,plot_avg_312_ct$graph,plot_avg_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)

# Plot for 2023 data
# Create overall titles
title_avg <- textGrob("2023 Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_2023_310_alk$graph, plot_2023_312_alk$graph, plot_2023_star_alk$graph,
              plot_2023_310_ct$graph,plot_2023_312_ct$graph,plot_2023_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)


# Plot for 2024 data
# Create overall titles
title_avg <- textGrob("2024 Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_2024_310_alk$graph, plot_2024_312_alk$graph, plot_2024_star_alk$graph,
              plot_2024_310_ct$graph,plot_2024_312_ct$graph,plot_2024_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)


####################### Saving the data ########################################
output_dir <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Heritability_Outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to save plots and data tables
save_results <- function(result, name, output_dir) {
  # Save the plot (assuming the first element is the plot)
  plot_file <- file.path(output_dir, paste0(name, ".png"))
  ggsave(plot_file, result[[1]])
  
  # Save the data tables (assuming the second and third elements are data tables)
  if (length(result) > 1) {
    for (i in 2:length(result)) {
      data_file <- file.path(output_dir, paste0(name, "_data_", i - 1, ".csv"))
      write.csv(result[[i]], data_file, row.names = FALSE)
    }
  }
}

# Save all results
save_results(plot_avg_310_alk, "plot_avg_310_alk", output_dir)
save_results(plot_avg_312_alk, "plot_avg_312_alk", output_dir)
save_results(plot_avg_star_alk, "plot_avg_star_alk", output_dir)
save_results(plot_avg_310_ct, "plot_avg_310_ct", output_dir)
save_results(plot_avg_312_ct, "plot_avg_312_ct", output_dir)
save_results(plot_avg_star_ct, "plot_avg_star_ct", output_dir)

save_results(plot_2023_310_alk, "plot_2023_310_alk", output_dir)
save_results(plot_2023_312_alk, "plot_2023_312_alk", output_dir)
save_results(plot_2023_star_alk, "plot_2023_star_alk", output_dir)
save_results(plot_2023_310_ct, "plot_2023_310_ct", output_dir)
save_results(plot_2023_312_ct, "plot_2023_312_ct", output_dir)
save_results(plot_2023_star_ct, "plot_2023_star_ct", output_dir)

save_results(plot_2024_310_alk, "plot_2024_310_alk", output_dir)
save_results(plot_2024_312_alk, "plot_2024_312_alk", output_dir)
save_results(plot_2024_star_alk, "plot_2024_star_alk", output_dir)
save_results(plot_2024_310_ct, "plot_2024_310_ct", output_dir)
save_results(plot_2024_312_ct, "plot_2024_312_ct", output_dir)
save_results(plot_2024_star_ct, "plot_2024_star_ct", output_dir)






######################
# Exploring the stats of data sets that needed to remove more individuals
######################

# These are the data that had things removed 
stats_2023_310_ct
stats_2024_310_ct
stats_2024_312_alk
stats_2024_312_ct

############### Function to make graphs showing removed individuals ############

removed_viz <- function(stats, pheno_data, phenotype){
  removed_individuals <- data.table(
    ID = unlist(strsplit(stats$RemovedIndividuals, split = ", "))
  )
  pheno_data_2 <-  pheno_data[!(pheno_data$ID %in% removed_individuals$ID), ]

  
  if(phenotype == "Alkaloids_Res"){
    p1 <- ggplot(pheno_data_2, aes(x = Alkaloids_Res)) + 
      geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = "Alkaloid Residuals Extra Removed", x = "Value", y = "Frequency") +
      theme_bw()
    
    p2 <- ggplot(pheno_data, aes(x = Alkaloids_Res)) + 
      geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = "Alkaloid Residuals", x = "Value", y = "Frequency") +
      theme_bw()
    
    residual_plot <- grid.arrange(p1, p2, ncol = 2, nrow = 1)
    return(residual_plot)
    
  }else if (phenotype == "Delta_CT_adj_Res") {
    
    p1 <- ggplot(pheno_data_2, aes(x = Delta_CT_adj_Res)) + 
      geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
      labs(title = "Delta CT Adj Residuals Extra Removed", x = "Value", y = "Frequency") +
      theme_bw()
    
    p2 <- ggplot(pheno_data, aes(x = Delta_CT_adj_Res)) + 
      geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
      labs(title = "Delta Adj Adj Residuals", x = "Value", y = "Frequency") +
      theme_bw()
    
    residual_plot <- grid.arrange(p1, p2, ncol = 2, nrow = 1)
    return(residual_plot)
    
  } else{
    cat("Paramaters are probably wrong\n")
  }
}
############################# End function #####################################
comparison1 <- removed_viz(stats_2023_310_ct, Residual_Data_23_outliars_rm_314x310, phenotype = "Delta_CT_adj_Res")
comparison2 <- removed_viz(stats_2024_310_ct, Residual_Data_24_outliars_rm_314x310, phenotype = "Delta_CT_adj_Res")
comparison3 <- removed_viz(stats_2024_312_alk, Residual_Data_24_outliars_rm_314x312, phenotype = "Alkaloids_Res")
comparison4 <- removed_viz(stats_2024_312_ct, Residual_Data_24_outliars_rm_314x312, phenotype = "Delta_CT_adj_Res")


######################## This is function testing #############################

# This function works with the geno_matrix
geno <- gsub("0/0", "0", geno)
geno <- gsub("0/1", "1", geno)
geno <- gsub("1/1", "2", geno)
return(as.numeric(geno))



# Load the VCF file
vcf_data_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf"
vcf_data <- read.vcfR(vcf_data_loc)
geno_matrix <- extract.gt(vcf_data, element = "GT")  
geno_matrix <- t(geno_matrix)



Residual_Data_24_outliars_rm_314x312 #Bad
Residual_Data_24_outliars_rm_314x310 #Good
# Ensure both datasets have the same IDs
common_IDs <- intersect(Residual_Data_24_outliars_rm_314x310$ID, rownames(geno_matrix))

# Subset data to include only common IDs
Residual_Data_24_outliars_rm_314x310 <- Residual_Data_24_outliars_rm_314x310[Residual_Data_24_outliars_rm_314x310$ID %in% common_IDs, ]
geno_matrix <- geno_matrix[common_IDs, ]

# Check if row names match
if (!all(rownames(geno_matrix) == Residual_Data_24_outliars_rm_314x310$ID)) {
  stop("Row names of geno_data do not match the ID in Residual_Data_24_outliars_rm_314x310.")
}

# Convert geno_matrix to data frame, and then convert genotypes
geno_data <- as.data.frame(geno_matrix)
geno_data <- geno_data %>% mutate_all(convert_genotypes)

# Convert geno_data back to matrix for kinship matrix calculation
geno_data <- as.matrix(geno_data)

# Create kinship matrix using Gmatrix function
kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")

# Ensure the data is in correct format (Residual_Data_24_outliars_rm_314x310 should have necessary columns)
Residual_Data_24_outliars_rm_314x310$ID <- as.character(Residual_Data_24_outliars_rm_314x310$ID)

# Example mixed model for heritability (adjust trait to be a column name)
model <- mmer(fixed = as.formula(paste("Alkaloids_Res", "~ 1")),
              random = ~ vsr(ID, Gu = kinship_matrix),
              rcov = ~ units,
              data = Residual_Data_24_outliars_rm_314x310)

# Summarize model results
model_summary <- summary(model)
varcomp_summary <- model_summary$varcomp

# Extract variance components for heritability calculation
var_ID <- varcomp_summary$VarComp[1]  # Genetic variance
var_residual <- varcomp_summary$VarComp[2]  # Residual variance

# Calculate heritability
heritability <- var_ID / (var_ID + var_residual)
heritability



################### Investigating problem heritabilities########################
# 0 heritability
stats_2024_312_alk
summary(Residual_Data_24_outliars_rm_314x312$Alkaloids_Res)
ggplot(Residual_Data_24_outliars_rm_314x312, aes(x = Alkaloids_Res)) +
  geom_histogram(binwidth = 5000, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Data set with 0 Heritability",
       x = "Alkaloids_Res",
       y = "Frequency") +
  theme_minimal()

# AN expected heritability
stats_2024_310_alk
summary(Residual_Data_24_outliars_rm_314x310$Alkaloids_Res)
ggplot(Residual_Data_24_outliars_rm_314x310, aes(x = Alkaloids_Res)) +
  geom_histogram(binwidth = 5000, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Dataset with expected Heritability",
       x = "Alkaloids_Res",
       y = "Frequency") +
  theme_minimal()

# Difference in variance and sd 

# Bad over Good
var(Residual_Data_24_outliars_rm_314x312$Alkaloids_Res) / var(Residual_Data_24_outliars_rm_314x310$Alkaloids_Res)
sd(Residual_Data_24_outliars_rm_314x312$Alkaloids_Res) / sd(Residual_Data_24_outliars_rm_314x310$Alkaloids_Res)
heatmap(as.matrix(kinship_matrix)) 

### Notes
# Number of individuals does not seem to cause the 0 in gentic heritability
# There is no missing data (e.g. NAs)

