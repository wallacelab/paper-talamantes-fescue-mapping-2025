# Objective: Make phenotype data into residual data and then subset all data 
# into all three data sets, 23, 24, avg all
# Uses PCA plots on tassel data to identify and remove genetic outliars

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
datafolder = "/home/darrian/Documents/Mapping_and_QTL/Data"

phenotype_Data <- read.table(paste0(datafolder,"/Phenotype_Data/All_Data_Filtered/phenotype_data.txt"), header = TRUE)
phenotype_Data$Delta_CT[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "314" & phenotype_Data$Data_Set != "315x320"), ]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Data_Set == "315x320"),]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Extraction_Date == "03/16/23"),]

all_Data_2024 <- read.csv(paste0(datafolder,"/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv"), header = TRUE)

# Load in Cross Identifier file
cross_list_loc <- ("../../Data/Lists/Parent_Progeny_Lists/314_Star_Cross_Parents.txt")
cross_list <- read.table(cross_list_loc, header = FALSE)

# Fixing naming conventions
phenotype_Data <- phenotype_Data %>%
  rename(Delta_CT_adj = Delta_CT)
phenotype_Data_2023 <- phenotype_Data
all_Data_2024 <- all_Data_2024 %>%
  rename(ID = Treatment)

# Loading the genotype data
vcf_data_loc <- "../../Data/VCF/all_snps_filtered_2.recode.vcf"
vcf_data <- read.vcfR(vcf_data_loc)
geno_matrix <- extract.gt(vcf_data, element = "GT")  
geno_matrix <- t(geno_matrix)

########## Getting residual data set for 2023 and 2024 avraged #################

# Removing data thats not star cross
list_314x310 <- read.table(paste0(datafolder,"/Lists/Parental_Lists/310x314_list.txt"))  
list_314x312 <- read.table(paste0(datafolder,"/Lists/Parental_Lists/312x314_list.txt"))
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
write.table(Residual_data_avg, file = paste0(datafolder, '/Phenotype_Data/Residual_data_avg.txt'), sep = '\t', row.names=FALSE)


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
write.table(phenotype_Data_23, file = paste0(datafolder, '/Phenotype_Data/2023_Data/Resisduals_Starcross_2023.txt'), sep = '\t', row.names=FALSE)


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
write.table(allpehnotype_data_export_24, file = paste0(datafolder, '/Phenotype_Data/2024_Data/Resisduals_Starcross_2024.txt'), sep = '\t', row.names=FALSE)

################################################################################
# Outliar exploration and removal
################################################################################

# load in previously saved data
Residual_data_avg <- read.table(paste0(datafolder,"/Phenotype_Data/Residual_data_avg.txt"), header = TRUE)
Residual_Data_23 <- read.table(paste0(datafolder,"/Phenotype_Data/2023_Data/Resisduals_Starcross_2023.txt"), header = TRUE)
Residual_Data_24 <- read.table(paste0(datafolder,"/Phenotype_Data/2024_Data/Resisduals_Starcross_2024.txt"), header = TRUE)

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

################ Creating a list of Geneetic outliars ##########################
######################## Small function for genotype data ######################
convert_genotypes <- function(geno) {
  geno <- gsub("0/0", "0", geno)
  geno <- gsub("0/1", "1", geno)
  geno <- gsub("1/1", "2", geno)
  return(as.numeric(geno))
}
############################## End Function ####################################
# Now we are going to look at a PCA and remove genetic outliars
geno_data <- as.data.frame(geno_matrix)
geno_data <- geno_data %>% mutate_all(convert_genotypes)
# Convert geno_data back to matrix for kinship matrix calculation
geno_data <- as.matrix(geno_data)
kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")
pca_result <- prcomp(kinship_matrix)


# Perform PCA on the kinship matrix (or SNP matrix if needed)
pca_scores <- pca_result$x  # Principal component scores
pca_df <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2])
pca_df <-merge(pca_df, cross_list, by.x = "row.names", by.y = "V1" )
pca_df <- pca_df %>% rename(Cross = V2)


# Create the scatter plot of PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Cross), size = 3, alpha = 0.3) +  # Plot the points
  geom_text(
    data = subset(pca_df, grepl("parent", Cross, ignore.case = TRUE)),  # Filter for "parent"
    aes(label = Row.names),
    vjust = -1,  # Adjust text position
    size = 3,
    color = "black"
  ) +
  labs(title = paste0("PCA"), x = "PC1", y = "PC2") +
  theme_minimal()

#Subset the data to get rid of weird genetic outliars
pca_df2 <- subset(pca_df, !((Cross == "314x312" & PC1 < 0) | (Cross == "314x310" & PC1 < -4)))

good_geentics <- pca_df2$Row.names
write.table(good_geentics, file = paste0(datafolder, '/Lists/Good_Genetics.txt'), sep = '\t', row.names=FALSE)


ggplot(pca_df2, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Cross), size = 3, alpha = 0.3) +  # Plot the points
  geom_text(
    data = subset(pca_df2, grepl("parent", Cross, ignore.case = TRUE)),  # Filter for "parent"
    aes(label = Row.names),
    vjust = -1,  # Adjust text position
    size = 3,
    color = "black"
  ) +
  labs(title = paste0("PCA"), x = "PC1", y = "PC2") +
  theme_minimal()


################ End Make list of  Genetic outliars ############################

#################Function to remove outliars ###################################

remove_outliers <- function(data, stds, good_geentics) {
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
  filtered_data <- filtered_data[filtered_data$ID %in% good_geentics, ]
  return(filtered_data)
}

############################## End Function ####################################


Residual_data_avg_outliars_rm <- remove_outliers(Residual_data_avg, 2.5, good_geentics)
Pheno_Graph(Residual_data_avg,"Residual Data of 2023 and 2024 Avaraged" )
Pheno_Graph(Residual_data_avg_outliars_rm,"Residual Data of 2023 and 2024 Avaraged" )
nrow(Residual_data_avg) - nrow(Residual_data_avg_outliars_rm)

Residual_Data_23_outliars_rm <- remove_outliers(Residual_Data_23, 2.5, good_geentics)
Pheno_Graph(Residual_Data_23,"Residual Data of 2023")
Pheno_Graph(Residual_Data_23_outliars_rm,"Residual Data of 2023")
nrow(Residual_Data_23) - nrow(Residual_Data_23_outliars_rm)

Residual_Data_24_outliars_rm <- remove_outliers(Residual_Data_24, 2.5, good_geentics)
Pheno_Graph(Residual_Data_24,"Residual Data of 2024")
Pheno_Graph(Residual_Data_24_outliars_rm,"Residual Data of 2024")
nrow(Residual_Data_24) - nrow(Residual_Data_24_outliars_rm)

write.table(Residual_data_avg_outliars_rm,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt"))
write.table(Residual_Data_23_outliars_rm,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt"))
write.table(Residual_Data_24_outliars_rm,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt"))


################################################################################
# Now seperating data to make all crosses.
################################################################################
#Reload Residual Data sets
Residual_data_avg_outliars_rm <- read.table(paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt"), header = TRUE)
Residual_Data_23_outliars_rm <- read.table(paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt"), header = TRUE)
Residual_Data_24_outliars_rm <- read.table(paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt"), header = TRUE)

# Loading in lists
list_314x310 <- read.table(paste0(datafolder,"/Lists/Parental_Lists/310x314_list.txt"))  
list_314x312 <- read.table(paste0(datafolder,"/Lists/Parental_Lists/312x314_list.txt"))

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
write.table(Residual_data_avg_outliars_rm_314x310,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x310.txt"))
write.table(Residual_data_avg_outliars_rm_314x312,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x312.txt"))

write.table(Residual_Data_23_outliars_rm_314x310,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x310.txt"))
write.table(Residual_Data_23_outliars_rm_314x312,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x312.txt"))

write.table(Residual_Data_24_outliars_rm_314x310,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x310.txt"))
write.table(Residual_Data_24_outliars_rm_314x312,paste0(datafolder,"/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x312.txt"))








































