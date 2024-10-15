# Objective: This script will make PCA graphs of 2024 phenotype 

# Import libraries 
library(tidyverse)
library(vcfR)
library(factoextra)
library(patchwork)
library(gridExtra)
library(grid)
library(vegan)
library(ape)
library(ggrepel)


# Loading in the data
phenotypes_2024_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv"
phenotypes_2024 <- read.csv(phenotypes_2024_loc, header = TRUE, strip.white=TRUE)

vcf_314x312_filtered <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x312_filtered.recode.vcf")
vcf_314x310_filtered <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_314x310_filtered.recode.vcf")
mcr50_filtered <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_Extra_filtered.recode.vcf")



vcf_314x312 <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x312.vcf")
vcf_314x310 <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/MCR50_snps_314x310.vcf")
mcr50 <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/Tall.fescue.MCR50.snps.vcf.gz")


################################################
# Doing PCA's
################################################
#Treatment to ID
phenotypes_2024 <- phenotypes_2024 %>% rename(ID = Treatment)
phenotypes_2024$ID <- gsub("-", "_", phenotypes_2024$ID)



####
# Creating a PCA function
####
PCA_plotter <- function(vcf_file, phenotypes, title){
  # Step 2: Extract the genotype matrix (individuals as columns)
  genotype_matrix <- extract.gt(vcf_file, as.numeric = TRUE)
  genotype_matrix[is.na(genotype_matrix)] <- 0  # Handle missing data
  
  # Step 3: Transpose the genotype matrix so that individuals are in rows
  genotype_df <- as.data.frame(t(genotype_matrix))
  
  # Assuming you have a phenotype data frame with IDs and crosses
  colnames(phenotypes)
  phenotypes$ID <- gsub("-", "_", phenotypes$ID)
  phenotypes$cross <- apply(phenotypes[, c("Mother", "Father")], 1, function(x) {
    if (all(is.na(x) | x == "")) {
      return("Parent")
    } else {
      return(paste(sort(x), collapse = "x"))
    }
  })  
  # Step 5: Set row names to match IDs from the phenotype data
  genotype_ids <- rownames(genotype_df)
  phenotypes_unique <- phenotypes %>%
    filter(ID %in% genotype_ids)
  phenotypes_unique <- unique(phenotypes_unique)
  
  # Ensure there are no extra phenotypes without corresponding genotype IDs
  unique_genotype_ids <- unique(rownames(genotype_df))
  unique_phenotype_ids <- unique(phenotypes_unique$ID)
  common_ids <- intersect(unique_phenotype_ids, unique_genotype_ids)
  # Filter both datasets to only include matching IDs
  
  genotype_df_filtered <- genotype_df[rownames(genotype_df) %in% common_ids, ]
  phenotypes_filtered <- phenotypes_unique[phenotypes_unique$ID %in% common_ids, ]
  rownames(genotype_df_filtered) <- phenotypes_filtered$ID
  
  # Step 6: Perform PCA on the transposed genotype data (individuals in rows)
  constant_columns <- sapply(genotype_df, function(x) length(unique(x)) == 1)
  genotype_df_filtered <- genotype_df[, !constant_columns]
  ncol(genotype_df_filtered)
  ncol(genotype_df)
  pca_result <- prcomp(genotype_df_filtered, center = TRUE, scale. = TRUE)
  
  # Step 7: Prepare PCA results for visualization
  pca_data <- as.data.frame(pca_result$x)
  
  # Add the 'cross' information to the PCA data
  pca_data$cross <- phenotypes$cross[match(rownames(pca_data), phenotypes$ID)]
  
  # Step 8: Create PCA plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cross)) +
    geom_point(size = 3) +
    ggtitle(title) +
    labs(x = "Principal Component 1", y = "Principal Component 2") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette = "Set1") 
  return(pca_plot)
}
#####
# End Function
#####

PCA_314x310 <- PCA_plotter(vcf_314x310, phenotypes_2024, "PCA 314x310")
PCA_314x310

PCA_314x312 <- PCA_plotter(vcf_314x312, phenotypes_2024, "PCA 314x312")
PCA_314x312

PCA_mcr50 <- PCA_plotter(mcr50, phenotypes_2024, "PCA MCR50")
PCA_mcr50


PCA_314x310_filtered <- PCA_plotter(vcf_314x310_filtered, phenotypes_2024, "PCA 314x310 Filtered")
PCA_314x310_filtered

PCA_314x312_filtered <- PCA_plotter(vcf_314x312_filtered, phenotypes_2024, "PCA 314x312 Filtered")
PCA_314x312_filtered

PCA_mcr50_filtered <- PCA_plotter(mcr50_filtered, phenotypes_2024, "PCA MCR50 Filtered")
PCA_mcr50_filtered


all_no_filters <- grid.arrange(
  arrangeGrob(PCA_mcr50, PCA_314x312, PCA_314x310, ncol = 3),
  left = textGrob("MCR50 No Filteres", rot = 90, gp = gpar(fontsize = 20))
)


all_filtered <- grid.arrange(
  arrangeGrob(PCA_mcr50_filtered, PCA_314x312_filtered, PCA_314x310_filtered, ncol = 3),
  left = textGrob("Darrian Filters", rot = 90, gp = gpar(fontsize = 20))
)

grid.arrange(
  arrangeGrob(all_no_filters, all_filtered, ncol = 1),
  nrow = 1
)
  





################################################################################
# Code to test the pca maker function
################################################################################



  # Step 2: Extract the genotype matrix (individuals as columns)
 genotype_matrix <- extract.gt(vcf_314x310_filtered, as.numeric = TRUE)
 genotype_df <- as.data.frame(t(genotype_matrix))
 dist(genotype_df)
 
 
 distance_matrix <- vegdist(genotype_matrix, method = "euclidean", na.rm = TRUE)
 pcoa_result <- pcoa(distance_matrix)
 
  
  # Assuming you have a phenotype data frame with IDs and crosses
  colnames(phenotypes_2024)
  phenotypes_2024$ID <- gsub("-", "_", phenotypes_2024$ID)
  phenotypes_2024$cross <- apply(phenotypes_2024[, c("Mother", "Father")], 1, function(x) {
    if (all(is.na(x) | x == "")) {
      return("Parent")
    } else {
      return(paste(sort(x), collapse = "x"))
    }
  })  
  # Step 5: Set row names to match IDs from the phenotype data
  genotype_ids <- rownames(genotype_df)
  phenotypes_unique <- phenotypes_2024 %>%
    filter(ID %in% genotype_ids)
  phenotypes_unique <- unique(phenotypes_unique)
  
  # Ensure there are no extra phenotypes without corresponding genotype IDs
  unique_genotype_ids <- unique(rownames(genotype_df))
  unique_phenotype_ids <- unique(phenotypes_unique$ID)
  common_ids <- intersect(unique_phenotype_ids, unique_genotype_ids)
  # Filter both datasets to only include matching IDs
  
  genotype_df_filtered <- genotype_df[rownames(genotype_df) %in% common_ids, ]
  phenotypes_filtered <- phenotypes_unique[phenotypes_unique$ID %in% common_ids, ]
  rownames(genotype_df_filtered) <- phenotypes_filtered$ID
  
  # Step 6: Perform PCA on the transposed genotype data (individuals in rows)
  constant_columns <- sapply(genotype_df, function(x) length(unique(x)) == 1)
  genotype_df_filtered <- genotype_df[, !constant_columns]
  ncol(genotype_df_filtered)
  ncol(genotype_df)
  pca_result <- prcomp(genotype_df_filtered, center = TRUE, scale. = TRUE)
  
  # Step 7: Prepare PCA results for visualization
  pca_data <- as.data.frame(pca_result$x)
  
  # Add the 'cross' information to the PCA data
  pca_data$cross <- phenotypes_2024$cross[match(rownames(pca_data), phenotypes_2024$ID)]
  
  # Step 8: Create PCA plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cross)) +
    geom_point(size = 3) +
    ggtitle(title) +
    labs(x = "Principal Component 1", y = "Principal Component 2") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette = "Set1") 


################################################################################
# ME being stupid and doing a PCA on the phenotypes
###############################################################################

# 
# 
# 
# # Finding the cross
# phenotypes_2024$cross <- apply(phenotypes_2024[, c("Mother", "Father")], 1, function(x) {
#   if (all(is.na(x)) || all(x == "")) {
#     return("Parent")
#   } else {
#     return(paste(sort(x), collapse = "x"))
#   }
# })
# 
# 
# colnames(phenotypes_2024)
# phenotype_data <- subset(phenotypes_2024, select = c("Delta_CT_OG","ng.g"))
# 
# 
# # Makes the data more normal 
# pca_result <- prcomp(phenotype_data, scale. = TRUE)  # Scaling the data is important
# pca_data <- as.data.frame(pca_result$x)
# 
# # Adding catagorical data back into the data
# pca_data$cross <- phenotypes_2024$cross
# pca_data$Data_Set <- phenotypes_2024$Data_Set
# pca_data$Alkaloid_Plate <- phenotypes_2024$Alkaloid_Plate
# 
# # PCA 
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Alkaloid_Plate)) +
#   geom_point() +
#   labs(title = "PCA of Numerical Data Colored by Category") +
#   theme_minimal()

################################################################################
# Graphs from tassel data
################################################################################
PCA_MCR50_Tassel_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/PCA_MCR50.txt"
PCA_MCR50_Tassel <- read.csv(PCA_MCR50_Tassel_loc, header = FALSE, strip.white=TRUE, skip = 2)

all_snps_filtered_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/PCA_all_Filtered.txt"
all_snps_filtered <- read.csv(all_snps_filtered_loc, header = FALSE, strip.white=TRUE, skip = 2)

PCA_MCR50_314x312_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/PCA_MCR50_314x312.txt"
PCA_MCR50_314x312 <- read.csv(PCA_MCR50_314x312_loc, header = FALSE, strip.white=TRUE, skip = 2)

PCA_filtered_314x312_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/PCA_all_314x312.txt"
PCA_filtered_314x312 <- read.csv(PCA_filtered_314x312_loc, header = FALSE, strip.white=TRUE, skip = 2)

PCA_all_filtered_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/PCA_all_Filtered.txt"
PCA_all_filtered <- read.csv(PCA_all_filtered_loc, header = FALSE, strip.white=TRUE, skip = 2)


PCA_tassel_data <- function(Tassel_PCA, title){
  Tassel_PCA <- Tassel_PCA[, -((ncol(Tassel_PCA)-4):ncol(Tassel_PCA))]
  colnames(Tassel_PCA) <- Tassel_PCA[1, ]
  Tassel_PCA <- Tassel_PCA[-1, ]
  Tassel_PCA$PC1 <- as.numeric(as.character(Tassel_PCA$PC1))
  Tassel_PCA$PC2 <- as.numeric(as.character(Tassel_PCA$PC2))
  
  # Making ggplot with data
  TasselPlot <- ggplot(Tassel_PCA, aes(x = PC1, y = PC2, color = Cross)) +
    geom_point() +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(data = subset(as.data.frame(Tassel_PCA), Tassel_PCA$Cross == "parent"),
                    aes(label = Taxa),
                    max.overlaps = Inf,   
                    nudge_x = 0.1, 
                    nudge_y = 0.1)
  return(TasselPlot)
}



PCA_tassel_data(PCA_MCR50_Tassel, "MCR50 Tassel")
PCA_tassel_data(PCA_all_filtered, "DRT Filtered All SNPs")

PCA_tassel_data(PCA_MCR50_314x312, "MCR50 314x312")
PCA_tassel_data(PCA_filtered_314x312, "DRT Filtered 314x312")


PCA_MCR50_Tassel <- PCA_MCR50_Tassel[, -((ncol(PCA_MCR50_Tassel)-4):ncol(PCA_MCR50_Tassel))]
colnames(PCA_MCR50_Tassel) <- PCA_MCR50_Tassel[1, ]
PCA_MCR50_Tassel <- PCA_MCR50_Tassel[-1, ]


# Rgl package scatter plot 3d
# Colorize them by the maternal parent since that is the one we know well. 
# Get metrics on filtered data sets, like % missing, MAF, avg depth, how many SNPs, heterogeneity







