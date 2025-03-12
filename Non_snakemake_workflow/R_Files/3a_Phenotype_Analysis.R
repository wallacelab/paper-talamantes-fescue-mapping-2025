# Purpose: This takes the 23 and 24 phenotype data and analyses them for the paper

library(tidyverse)
library(multcompView)
#library(gvlma) # Dont think I used this either
#library(mbQTL) #does not work? WHat did I use this for?
library(ggpubr)
library(car)
library(reshape2)
library(gridExtra)




##############################################
# Removing all outliars
##############################################
# Loading in the data
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

# The two data sets you need are all_Data_2024 and phenotype_Data_2023
phenotype_Data <- read.table(paste0(data_folder,"/Phenotype_Data/All_Data_Filtered/phenotype_data.txt") , header = TRUE)
phenotype_Data$Delta_CT[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "314" & phenotype_Data$Data_Set != "315x320"), ]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Data_Set == "315x320"),]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Extraction_Date == "03/16/23"),]
good_genetics <- read.table(paste0(data_folder, "/Lists/Good_Genetics.txt"), header = TRUE, )
all_Data_2024 <- read.csv(paste0(data_folder,"/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv"), header = TRUE)

# Fixing naming conventions
phenotype_Data <- phenotype_Data %>%
  rename(Delta_CT_adj = Delta_CT)
phenotype_Data_2023 <- phenotype_Data
all_Data_2024 <- all_Data_2024 %>%
  rename(ID = Treatment)


# Removing data thats not star cross
list_314x310 <- read.table(paste0(data_folder,"/Lists/Parental_Lists/310x314_list.txt"))  
list_314x312 <- read.table(paste0(data_folder,"/Lists/Parental_Lists/312x314_list.txt"))
new_rows <- data.frame(V1 = c(301, 302, 303, 304, 305, 306, 307, 308, 310, 312, 313, 314, 315, 316, 318, 319, 320))
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

# Removing any rows that were genetic outliars
colnames(good_genetics)[1] <- "ID"
gg_phenotypes23_24 <- merge(phenotypes23_24, good_genetics, by = "ID")
head(gg_phenotypes23_24,15)




# Now that we have residuals we have to get avarages and export the data
avarage_data <- function(data_set){
  Alkaloid_residuals_avaraged <- data_set %>%
    group_by(ID) %>%
    summarise(Alkaloids_Res_avg = mean(Alkaloids_Res, na.rm = TRUE)) %>%
    filter(!is.na(Alkaloids_Res_avg))
  
  DeltaCT_adj_residuals_avaraged <- data_set %>%
    group_by(ID) %>%
    summarise(DeltaCT_adj_Res_avg = mean(Delta_CT_adj_Res, na.rm = TRUE)) %>%
    filter(!is.na(DeltaCT_adj_Res_avg))
  
  DeltaCT_OG_residuals_avaraged <- data_set %>%
    group_by(ID) %>%
    summarise(DeltaCT_OG_Res_avg = mean(Delta_CT_OG_Res, na.rm = TRUE)) %>%
    filter(!is.na(DeltaCT_OG_Res_avg))
  
  Residual_data_avg <- merge(Alkaloid_residuals_avaraged, 
                             DeltaCT_adj_residuals_avaraged, 
                             by = "ID")
  
  Residual_data_avg <- merge(Residual_data_avg, 
                             DeltaCT_OG_residuals_avaraged, 
                             by = "ID")
  return(Residual_data_avg)
}

Residual_data_avg <- avarage_data(phenotypes23_24)
gg_Residual_data_avg <- avarage_data(gg_phenotypes23_24)

# Remove phenotype outliars
pg_Residual_data_avg <- gg_Residual_data_avg  # Start with the original data frame
for (col_name in colnames(gg_Residual_data_avg)[-1]) {  # Exclude the first column (ID)
  # Calculate the mean and standard deviation for the column
  mean_value <- mean(gg_Residual_data_avg[[col_name]], na.rm = TRUE)
  sd_value <- sd(gg_Residual_data_avg[[col_name]], na.rm = TRUE)
  
  # Define the threshold for removal (mean ± 2.5 * standard deviation)
  threshold_lower <- mean_value - 2.5 * sd_value
  threshold_upper <- mean_value + 2.5 * sd_value
  
  # Remove rows where the value is outside of the defined range for this column
  pg_Residual_data_avg <- pg_Residual_data_avg[pg_Residual_data_avg[[col_name]] >= threshold_lower & pg_Residual_data_avg[[col_name]] <= threshold_upper, ]
}

# Remove not important parents
abscent_parents <- data.frame(ID = c("301", "302", "304", "305", "307", "308", "313", "315", "316", "318", "319", "320"))
pg_Residual_data_avg <- pg_Residual_data_avg[!pg_Residual_data_avg$ID %in% abscent_parents$ID,]


#Save the phenotype data sets
write.table(Residual_data_avg, file = paste0(data_folder,"/Phenotype_Data/Residual_data_avg.txt"), sep = '\t', row.names=FALSE)
write.table(pg_Residual_data_avg, file = paste0(data_folder,"/Phenotype_Data/pg_Residual_data_avg"), sep = '\t', row.names=FALSE)

# Graphs to look at the avraged residual data with outliars
p1 <- ggplot(Residual_data_avg, aes(x = Alkaloids_Res_avg)) + 
  geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Alkaloids", x = "Value", y = "Frequency") +
  theme_bw()

p2 <- ggplot(Residual_data_avg, aes(x = DeltaCT_adj_Res_avg)) + 
  geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Delta CT adj", x = "Value", y = "Frequency") +
  theme_bw()

p3 <- ggplot(Residual_data_avg, aes(x = DeltaCT_OG_Res_avg)) + 
  geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Delta CT OG", x = "Value", y = "Frequency") +
  theme_bw()

grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

# Graphs to look at the avaraged residual data with genetic and phenotypic outliars removed
p1 <- ggplot(pg_Residual_data_avg, aes(x = Alkaloids_Res_avg)) + 
  geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Alkaloids \n Outliars Removed", x = "Value", y = "Frequency") +
  theme_bw()

p2 <- ggplot(pg_Residual_data_avg, aes(x = DeltaCT_adj_Res_avg)) + 
  geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Delta CT adj \n Outliars Removed", x = "Value", y = "Frequency") +
  theme_bw()

p3 <- ggplot(pg_Residual_data_avg, aes(x = DeltaCT_OG_Res_avg)) + 
  geom_histogram(binwidth = .25, fill = "yellow", color = "black", alpha = 0.7) +
  labs(title = "Values for Residuals of Delta CT OG \n Outliars Removed", x = "Value", y = "Frequency") +
  theme_bw()

grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

############### Scatter plots #######################
head(pg_Residual_data_avg)
phenotypes_23 <- subset(phenotypes23_24, Year == 2023)  
phenotypes_24 <- subset(phenotypes23_24, Year == 2024)  

# linear model
model <- lm(Alkaloids_Res_avg ~ DeltaCT_adj_Res_avg, data = pg_Residual_data_avg)
rsq <- summary(model)$r.squared #This is low so I think I should use spearman.


# Function to calculate r squared and make scatter plot
scatterplot_phenos <- function(dataset, Alkaloidcol, DeltaCTcol, Title) {
  # Perform Spearman correlation
  spearmodel <- cor.test(dataset[[Alkaloidcol]], 
                         dataset[[DeltaCTcol]], 
                         method = "spearman", 
                         use = "complete.obs")  # Ignore NAs
  rho <- spearmodel$estimate
  rsq <- rho^2
  p <- spearmodel$p.value
  
  # Create the scatter plot
  plot1 <- ggplot(dataset, aes(x = .data[[DeltaCTcol]], y = .data[[Alkaloidcol]])) +
    geom_point() +  # Scatter plot points
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
    annotate("text", 
             x = min(dataset[[DeltaCTcol]], na.rm = TRUE), 
             y = max(dataset[[Alkaloidcol]], na.rm = TRUE), 
             label = paste("Spearman R-squared = ", round(rsq, 3)), 
             hjust = 0, vjust = -1, size = 5, color = "red") +
    annotate("text", 
             x = min(dataset[[DeltaCTcol]], na.rm = TRUE), 
             y = max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
               0.05 * (max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
                         min(dataset[[Alkaloidcol]], na.rm = TRUE)), 
             label = paste("P-value = ", format(p, digits = 3, scientific = TRUE)), 
             hjust = 0, vjust = -1, size = 5, color = "red") +
    labs(title = Title, 
         x = "Efficiency Adjusted CT Ratio", 
         y = "Residual Alkaloids") + 
    theme_bw() +
    theme(text = element_text(size = 20))
  
  return(plot1)
}


scatterplot_phenos(pg_Residual_data_avg,"Alkaloids_Res_avg","DeltaCT_adj_Res_avg", "Residual Alkaloids vs Efficiency Adjusted CT Ratio")
scatterplot_phenos(phenotypes_23,"Alkaloids_Res","Delta_CT_adj_Res", "2023")
scatterplot_phenos(phenotypes_24,"Alkaloids_Res","Delta_CT_adj_Res", "2024")

IDs <- pg_Residual_data_avg$ID
write.table(IDs, file = paste0(data_folder,"/Lists/3a_Geno_List_Outliars_Removed.txt"), row.names = FALSE)




##############################################
# Creating Residual Data for only the 2023 data set
##############################################
# Loading in the data
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

# The two data sets you need are all_Data_2024 and phenotype_Data_2023
phenotype_Data <- read.table(paste0(data_folder,"/Phenotype_Data/All_Data_Filtered/phenotype_data.txt") , header = TRUE)
phenotype_Data$Delta_CT[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "314" & phenotype_Data$Data_Set != "315x320"), ]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Data_Set == "315x320"),]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Extraction_Date == "03/16/23"),]

# Fixing naming conventions
phenotype_Data <- phenotype_Data %>%
  rename(Delta_CT_adj = Delta_CT)
phenotype_Data_2023 <- phenotype_Data


#Removing batch effects and leaving only residuals.
lm_model_alk <- lm(ng.g ~ Harvest_Date + Alkaloid_Plate, data = phenotype_Data_2023, na.action = na.exclude)
phenotype_Data_2023$Alkaloids_Res <- resid(lm_model_alk)

lm_model_CT_OG <- lm(Delta_CT_OG ~ Data_Set  + Standard  + Extraction_Date + Extractor + Harvest_Date, data = phenotype_Data_2023, na.action = na.exclude)
phenotype_Data_2023$Delta_CT_OG_Res <- resid(lm_model_CT_OG)

lm_model_CT_adj <- lm(Delta_CT_adj ~  Data_Set  + Standard  + Extraction_Date + Extractor + Harvest_Date, data = phenotype_Data_2023, na.action = na.exclude)
phenotype_Data_2023$Delta_CT_adj_Res <- resid(lm_model_CT_adj)

# Saves table with MetaData
write.table(phenotype_Data_2023, file = paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_data_and_metadata.txt"), sep = '\t', row.names=FALSE)

# Only keeping phenotype data
phenotype_Data_2023 <- phenotype_Data_2023[, c(1, 13:15)]
# Remove phenotype outliars
pg_phenotype_Data_2023 <- phenotype_Data_2023  # Start with the original data frame
for (col_name in colnames(phenotype_Data_2023)[-1]) {  # Exclude the first column (ID)
  # Calculate the mean and standard deviation for the column
  mean_value <- mean(phenotype_Data_2023[[col_name]], na.rm = TRUE)
  sd_value <- sd(phenotype_Data_2023[[col_name]], na.rm = TRUE)
  
  # Define the threshold for removal (mean ± 2.5 * standard deviation)
  threshold_lower <- mean_value - 2.5 * sd_value
  threshold_upper <- mean_value + 2.5 * sd_value
  
  # Remove rows where the value is outside of the defined range for this column
  pg_phenotype_Data_2023 <- pg_phenotype_Data_2023[pg_phenotype_Data_2023[[col_name]] >= threshold_lower & pg_phenotype_Data_2023[[col_name]] <= threshold_upper, ]
}
pg_phenotype_Data_2023 <- pg_phenotype_Data_2023[rowSums(is.na(pg_phenotype_Data_2023)) != ncol(pg_phenotype_Data_2023), ]
nrow(pg_phenotype_Data_2023)

#Save the phenotype data sets
write.table(pg_phenotype_Data_2023, file = paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_only.txt"), sep = '\t', row.names=FALSE)






































################################################################################
# Useless from beyond this point
################################################################################
###
### Subsetting data for other heritability
###

phenotypes23_24 <- rbind(allpehnotype_data_export_23,allpehnotype_data_export_24)
phenotypes23_24$ID <- gsub("-", "_", phenotypes23_24$ID)
head(phenotypes23_24,15)

# Here we have to subset the data in 8 different ways.

#Removing batch effects and leaving only residuals.

#I suggest making this into a function that we can use many times.
lm_model_alk <- lm(ng.g ~ Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Alkaloids_Res <- resid(lm_model_alk)

lm_model_CT_OG <- lm(Delta_CT_OG ~ Data_Set + Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Delta_CT_OG_Res <- resid(lm_model_CT_OG)

lm_model_CT_adj <- lm(Delta_CT_adj ~ Data_Set + Year, data = phenotypes23_24, na.action = na.exclude)
phenotypes23_24$Delta_CT_adj_Res <- resid(lm_model_CT_adj)

head(phenotypes23_24,15)


##########################################
# Graphs to explore the 2023 vs 2024 data 
#########################################






################################################################################
###### Making graphs to explore the pehnotypic data
################################################################################
tassel_2024_data <- read.table(file = paste0(data_folder,"/Phenotype_Data/2024_Data/tassel_2024_data.txt"), sep = '\t', header = TRUE)
all_Data_2024 <- read.csv(paste0(data_folder,"/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv") , header = TRUE)
parents_2024 <- all_Data_2024[, c("Treatment", "Mother", "Father")]
parents_2024 <- parents_2024 %>% rename(ID = Treatment)
merged_2024 <- merge(tassel_2024_data, parents_2024, by = "ID")
merged_2024$Parent_Combination <- apply(merged_2024[, c("Mother", "Father")], 1, function(x) paste(sort(x), collapse = " x "))

head(merged_2024)



ggplot(merged_2024, aes(x = Delta_CT_adj, fill = Parent_Combination)) +
  geom_histogram(position = "stack", binwidth = .5, color = "black") +
  labs(title = "Histogram of Delta_CT_adj Colored by Parent Combination", x = "Delta_CT_adj Values", y = "Count") +
  theme_minimal()

ggplot(merged_2024, aes(x = ng.g, fill = Parent_Combination)) +
  geom_histogram(position = "stack", binwidth = 5000, color = "black") +
  labs(title = "Histogram of Alkaloids Colored by Parent Combination", x = "Delta_CT_adj Values", y = "Count") +
  theme_minimal()

