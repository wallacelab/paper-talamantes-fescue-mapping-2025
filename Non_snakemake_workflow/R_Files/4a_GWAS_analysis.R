#Purpose: Take the output of tassel and graph it. 
library(qqman)
library(grid)
library(gridExtra)
library(ggpubr)
library(tidyverse)



# Load in Data 
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/2025_all_analysis/MLM_all_snps_no_outliars_residuals_stats.txt"), header = TRUE, sep = '\t')


########### 2023 Phenos with 3000 OG Genotypes Alkaloids  ######################
MLM_DRT_Filters_residuals_loc <- "/home/darrian/Documents/Mapping_and_QTL/Data/Tassel_Outputs/2023_only/MLM_all3000_genos_residual_phenos_stats.txt"
MLM_DRT_Filters_residuals <- read.table(MLM_DRT_Filters_residuals_loc, header = TRUE, sep = '\t')
# Changing the chr to numbers
generate_id <- function(x) {
  if (grepl("^FACHR\\d+[A-Za-z]\\d+$", x)) {
    # Extract Number, Letter, Number
    num_letter_num <- gsub("FACHR(\\d+)([A-Za-z])(\\d+)", "\\1-\\2-\\3", x)
    
    # Convert Letter to a Number (A=1, B=2, ..., Z=26, a=27, b=28, ..., z=52)
    parts <- unlist(strsplit(num_letter_num, "-"))
    num1 <- as.numeric(parts[1])
    letter_val <- match(parts[2], c(LETTERS, letters))  # Convert letter to number
    num2 <- as.numeric(parts[3])
    
    # Create a unique numeric ID
    unique_id <- num1 * 1000 + letter_val * 100 + num2  # Ensures uniqueness
    
    return(paste0("1.", unique_id))
    
  } else if (grepl("^FACHR\\d+[A-Za-z]$", x)) {
    # Case: FACHR#X (Number, Letter) → Assign "9" as default last digit
    num_letter <- gsub("FACHR(\\d+)([A-Za-z])", "\\1-\\2", x)
    
    parts <- unlist(strsplit(num_letter, "-"))
    num1 <- as.numeric(parts[1])
    letter_val <- match(parts[2], c(LETTERS, letters))  # Convert letter to number
    
    # Create a unique numeric ID
    unique_id <- num1 * 1000 + letter_val * 100 + 9  # Default last digit as 9
    
    return(paste0("1.", unique_id))
    
  } else if (grepl("^SCAFFOLD\\d+$", x)) {
    # Case: SCAFFOLDX (Number only)
    scaffold_num <- as.numeric(gsub("SCAFFOLD", "", x))
    return(paste0("2.", scaffold_num))
    
  } else {
    return(NA)  # If it doesn't match any pattern
  }
}

# Apply function to the column
MLM_DRT_Filters_residuals$unique_id <- sapply(MLM_DRT_Filters_residuals$Chr, generate_id)
MLM_DRT_Filters_residuals$Chr <- MLM_DRT_Filters_residuals$unique_id
MLM_DRT_Filters_residuals$unique_id <- NULL
MLM_DRT_Filters_residuals$Chr <- as.numeric(MLM_DRT_Filters_residuals$Chr)
#Chosing the phenotype
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[MLM_DRT_Filters_residuals$Trait == "Alkaloids_Res", ] #Delta_CT_adj_Res or Alkaloids_Res
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[-c(1), ]
# BOnferonii line
alpha <- 0.05
num_tests <- nrow(MLM_DRT_Filters_residuals)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_DRT_Filters_residuals$padj <- p.adjust(MLM_DRT_Filters_residuals$p, method = "BH")
FDR_threshold <- max(MLM_DRT_Filters_residuals$p[MLM_DRT_Filters_residuals$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
Alkaloids <- manhattan(MLM_DRT_Filters_residuals, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
                       ylim = c(0, 6), main = "MLM on First Genotypes with 2023 Alkaloid Concentration",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))



################ 2023 Phenos with 3000 OG Genotypes Biomass ####################
MLM_DRT_Filters_residuals_loc <- "/home/darrian/Documents/Mapping_and_QTL/Data/Tassel_Outputs/2023_only/MLM_all3000_genos_residual_phenos_stats.txt"
MLM_DRT_Filters_residuals <- read.table(MLM_DRT_Filters_residuals_loc, header = TRUE, sep = '\t')
# Changing the chr to numbers
generate_id <- function(x) {
  if (grepl("^FACHR\\d+[A-Za-z]\\d+$", x)) {
    # Extract Number, Letter, Number
    num_letter_num <- gsub("FACHR(\\d+)([A-Za-z])(\\d+)", "\\1-\\2-\\3", x)
    
    # Convert Letter to a Number (A=1, B=2, ..., Z=26, a=27, b=28, ..., z=52)
    parts <- unlist(strsplit(num_letter_num, "-"))
    num1 <- as.numeric(parts[1])
    letter_val <- match(parts[2], c(LETTERS, letters))  # Convert letter to number
    num2 <- as.numeric(parts[3])
    
    # Create a unique numeric ID
    unique_id <- num1 * 1000 + letter_val * 100 + num2  # Ensures uniqueness
    
    return(paste0("1.", unique_id))
    
  } else if (grepl("^FACHR\\d+[A-Za-z]$", x)) {
    # Case: FACHR#X (Number, Letter) → Assign "9" as default last digit
    num_letter <- gsub("FACHR(\\d+)([A-Za-z])", "\\1-\\2", x)
    
    parts <- unlist(strsplit(num_letter, "-"))
    num1 <- as.numeric(parts[1])
    letter_val <- match(parts[2], c(LETTERS, letters))  # Convert letter to number
    
    # Create a unique numeric ID
    unique_id <- num1 * 1000 + letter_val * 100 + 9  # Default last digit as 9
    
    return(paste0("1.", unique_id))
    
  } else if (grepl("^SCAFFOLD\\d+$", x)) {
    # Case: SCAFFOLDX (Number only)
    scaffold_num <- as.numeric(gsub("SCAFFOLD", "", x))
    return(paste0("2.", scaffold_num))
    
  } else {
    return(NA)  # If it doesn't match any pattern
  }
}

# Apply function to the column
MLM_DRT_Filters_residuals$unique_id <- sapply(MLM_DRT_Filters_residuals$Chr, generate_id)
MLM_DRT_Filters_residuals$Og_Chr <- MLM_DRT_Filters_residuals$Chr
MLM_DRT_Filters_residuals$Chr <- MLM_DRT_Filters_residuals$unique_id
MLM_DRT_Filters_residuals$Chr <- as.numeric(MLM_DRT_Filters_residuals$Chr)
MLM_DRT_Filters_residuals$unique_id <- NULL
#Chosing the phenotype
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[MLM_DRT_Filters_residuals$Trait == "Delta_CT_adj_Res", ] #Delta_CT_adj_Res or Alkaloids_Res
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[-c(1), ]
# BOnferonii line
alpha <- 0.05
num_tests <- nrow(MLM_DRT_Filters_residuals)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_DRT_Filters_residuals$padj <- p.adjust(MLM_DRT_Filters_residuals$p, method = "BH")
FDR_threshold <- max(MLM_DRT_Filters_residuals$p[MLM_DRT_Filters_residuals$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
Alkaloids <- manhattan(MLM_DRT_Filters_residuals, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
                       ylim = c(0, 10), main = "MLM on First Genotypes with 2023 Fungal to Plant Ratio",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))

sig_snps <- subset(MLM_DRT_Filters_residuals, padj <= .05) 
write.table(sig_snps, paste0(data_folder , "/Significant_Snps/MLM_3000_genos_biomass.tsv"), row.names = FALSE)









################################################################################
# Redoing half sib phenotypes in greater deapth 
################################################################################


################ Half Sibs 2023 Alkaloids with new genos ##########################


# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_half_sib_2023.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Alkaloids_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
  ylim = c(0, 6), main = "MLM on 2023 Alkaloid Half Sibs \n and Freedom Marker Genos",
  genomewideline = -log10(bonferroni_threshold),
  suggestiveline = -log10(FDR_threshold),
  col = c("blue", "red", "darkgrey", "purple"))

################ Half Sibs 2024 Alkaloids with new genos ##########################
# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_stats_half_sib_2024.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Alkaloids_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          ylim = c(0, 6), main = "MLM on 2024 Alkaloid Half Sibs \n and Freedom Marker Genos",
          genomewideline = -log10(bonferroni_threshold),
          suggestiveline = -log10(FDR_threshold),
          col = c("blue", "red", "darkgrey", "purple"))


################ Half Sibs Averaged Alkaloids with new genos ##########################
# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_stats_half_sib_average.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Alkaloids_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          ylim = c(0, 6), main = "MLM on Averaged Alkaloid Half Sibs \n and Freedom Marker Genos",
          genomewideline = -log10(bonferroni_threshold),
          suggestiveline = -log10(FDR_threshold),
          col = c("blue", "red", "darkgrey", "purple"))



################ Half Sibs 2023 Fungal Biomass with new genos ##########################


# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_half_sib_2023.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Delta_CT_adj_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          ylim = c(0, 6), main = "MLM on 2023 Fungal Biomass Half Sibs \n and Freedom Marker Genos",
          genomewideline = -log10(bonferroni_threshold),
          suggestiveline = -log10(FDR_threshold),
          col = c("blue", "red", "darkgrey", "purple"))

################ Half Sibs 2024 Fungal Biomass with new genos ##########################
# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_stats_half_sib_2024.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Delta_CT_adj_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          ylim = c(0, 6), main = "MLM on 2024 Fungal Biomass Half Sibs \n and Freedom Marker Genos",
          genomewideline = -log10(bonferroni_threshold),
          suggestiveline = -log10(FDR_threshold),
          col = c("blue", "red", "darkgrey", "purple"))


################ Half Sibs Averaged Fungal Biomass with new genos ##########################
# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/Half_Sibs/mlm_new_genos_stats_half_sib_average.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Delta_CT_adj_Res", ] #CHANGE THIS
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]

#Making of manhattan plots
alpha <- 0.05
num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          ylim = c(0, 6), main = "MLM on Averaged Fungal Biomass Half Sibs \n and Freedom Marker Genos",
          genomewideline = -log10(bonferroni_threshold),
          suggestiveline = -log10(FDR_threshold),
          col = c("blue", "red", "darkgrey", "purple"))




























################################################################################
# Depcricated (used less filtered phenotypes)
################################################################################


# ### Avarage Alkaloids
# # Data prep
# MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
# MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
# MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "DeltaCT_adj_Res_avg", ]
# MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]
# 
# #Making of manhattan plots
# alpha <- 0.05
# num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
# bonferroni_threshold <- alpha / num_tests
# # Making FDR threshold
# MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
# FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
# if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
#   FDR_threshold <- 1e-5
# }
# print(FDR_threshold)
# CT_ratio <- manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
#                       ylim = c(0, 6), main = "MLM on CT Ratio Residuals Avaraged From 2023 and 2024",
#                       genomewideline = -log10(bonferroni_threshold),
#                       suggestiveline = -log10(FDR_threshold),
#                       col = c("blue", "red", "darkgrey", "purple"))
# 
# ### Avarage Adjusted Delta CT
# 
# # Data prep
# MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/2025_all_analysis/MLM_all_snps_no_outliars_residuals_stats.txt"), header = TRUE, sep = '\t')
# MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
# MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
# MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Alkaloids_Res_avg", ]
# MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[-c(1), ]
# 
# #Making of manhattan plots
# alpha <- 0.05
# num_tests <- nrow(MLM_all_snps_no_outliars)  # Number of SNPs or tests
# bonferroni_threshold <- alpha / num_tests
# # Making FDR threshold
# MLM_all_snps_no_outliars$padj <- p.adjust(MLM_all_snps_no_outliars$p, method = "BH")
# FDR_threshold <- max(MLM_all_snps_no_outliars$p[MLM_all_snps_no_outliars$padj <= alpha], na.rm = TRUE)
# if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
#   FDR_threshold <- 1e-5
# }
# print(FDR_threshold)
# Alkaloids <- manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
#                        ylim = c(0, 6), main = "MLM on Alkaloid Level Residuals Avaraged From 2023 and 2024",
#                        genomewideline = -log10(bonferroni_threshold),
#                        suggestiveline = -log10(FDR_threshold),
#                        col = c("blue", "red", "darkgrey", "purple"))
# 
# 
