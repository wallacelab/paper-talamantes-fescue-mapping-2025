#Purpose: Take the output of tassel and graph it. 
library(qqman)
library(grid)
library(gridExtra)
library(ggpubr)
library(tidyverse)



# Load in Data 
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/2025_all_analysis/MLM_all_snps_no_outliars_residuals_stats.txt"), header = TRUE, sep = '\t')












################################################################################
# Make Manhattan Plots
################################################################################

### Avarage Alkaloids
# Data prep
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "DeltaCT_adj_Res_avg", ]
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
CT_ratio <- manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
                       ylim = c(0, 6), main = "MLM on CT Ratio Residuals Avaraged From 2023 and 2024",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))

### Avarage Adjusted Delta CT

# Data prep
MLM_all_snps_no_outliars <- read.table(paste0(data_folder, "/Tassel_Outputs/2025_all_analysis/MLM_all_snps_no_outliars_residuals_stats.txt"), header = TRUE, sep = '\t')
MLM_all_snps_no_outliars$Chr<-gsub("SCAFFOLD_","",as.character(MLM_all_snps_no_outliars$Chr))
MLM_all_snps_no_outliars$Chr<- as.numeric(MLM_all_snps_no_outliars$Chr)
MLM_all_snps_no_outliars <- MLM_all_snps_no_outliars[MLM_all_snps_no_outliars$Trait == "Alkaloids_Res_avg", ]
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
Alkaloids <- manhattan(MLM_all_snps_no_outliars, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
                       ylim = c(0, 6), main = "MLM on Alkaloid Level Residuals Avaraged From 2023 and 2024",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))







########### Example 
MLM_DRT_Filters_residuals_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/MLM_DRT_filters_Phenotype_residuals_avaraged.txt"
MLM_DRT_Filters_residuals <- read.table(MLM_DRT_Filters_residuals_loc, header = TRUE, sep = '\t')
MLM_DRT_Filters_residuals$Chr<-gsub("SCAFFOLD_","",as.character(MLM_DRT_Filters_residuals$Chr))
MLM_DRT_Filters_residuals$Chr<- as.numeric(MLM_DRT_Filters_residuals$Chr)
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[MLM_DRT_Filters_residuals$Trait == "Alkaloids_Res_avg", ]
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
                       ylim = c(0, 6), main = "MLM on Alkaloid Level Residuals Avaraged From 2023 and 2024",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))

