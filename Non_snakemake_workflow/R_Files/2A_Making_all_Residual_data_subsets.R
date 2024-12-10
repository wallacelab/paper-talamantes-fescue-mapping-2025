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

################################################################################
# Outliar exploration and removal
################################################################################

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








































