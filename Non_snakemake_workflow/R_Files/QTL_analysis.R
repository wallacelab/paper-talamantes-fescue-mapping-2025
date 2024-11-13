# Objective: This script will do QTL analysis on residual 2023 and 2024 data. It uses DRT filtered genotypes
# Author: Darrian Talamantes


# Library insertion
library(qtl)
library(vcfR)
library(tibble)



# Loading files

# Turning the vcf file into a format RQTL can work with
vcf_data <- read.vcfR("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/VCF/all_snps_filtered_2.recode.vcf")
gt_matrix <- extract.gt(vcf_data, element = "GT", as.numeric = TRUE)
gt_matrix <- as.data.frame(t(gt_matrix))
gt_matrix <- rownames_to_column(gt_matrix, var = "ID")

# Load in phenotype data
phenotypes <- read.table("/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Phenotypes_23_24_tassel.txt", skip = 2,header = TRUE)
avg_phenotypes <- subset(phenotypes, select = c("ID", "Delta_CT_adj_avg", "Delta_CT_OG_avg", "ng.g_avg"))

# Combine the data
combined_data <- merge(avg_phenotypes, gt_matrix, by = "ID")
write.csv(combined_data, "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/QTL_Files/combined_data.txt", row.names = FALSE)

#Turn into rqtl format
genotypes <- read.cross("csv", file = "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/QTL_Files/combined_data.txt", genotypes=c("AA", "AB", "BB"))


# Do HK regression
results <- scanone(combined_data, method = "hk")  # 'hk' for Haley-Knott regression
plot(results)
perm_results <- scanone(combined_data, method = "hk", n.perm = 1000)
summary(perm_results)






