# THis file is created to check heritabilities for the half-sib 2024 data that
# failed to get heritability in tassel.


################################################################################
# heritability function
################################################################################

# Load necessary library
library(rrBLUP)

# Define Data Folder
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

# Define the function
calculate_heritability <- function(hapmap_file, phenotype_file, data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data") {
  
  #----------------------------
  # 1. Load your Hapmap file
  #----------------------------
  hapmap <- read.delim(hapmap_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  
  #----------------------------
  # 2. Extract just the genotypes
  #----------------------------
  # Remove first 11 columns (metadata), keep sample genotypes
  geno_data <- hapmap[, -(1:11)]  
  rownames(geno_data) <- hapmap$`rs#`  # set SNP names as rownames
  
  #----------------------------
  # 3. Convert alleles to numeric (0/1/2)
  #----------------------------
  # Define conversion function
  convert_geno <- function(x) {
    if (x %in% c("A", "C", "G", "T")) return(0)  # homozygous major
    if (grepl("/", x)) return(1)                # heterozygous (like A/T)
    if (x %in% c("N", "./.")) return(NA)         # missing data
    return(2)  # if somehow minor homozygous
  }
  
  # Apply to each cell
  geno_numeric <- as.data.frame(apply(geno_data, c(1,2), convert_geno))
  
  # Transpose to have samples as rows
  geno_numeric <- as.data.frame(t(geno_numeric))
  
  # Fix the rownames
  rownames(geno_numeric) <- sub("^X", "", rownames(geno_numeric))
  
  #----------------------------
  # 4. Match with phenotype
  #----------------------------
  phenotype  <- read.table(phenotype_file, sep = '\t', header = TRUE)
  
  # Make sure order matches (IMPORTANT!)
  geno_numeric <- geno_numeric[match(phenotype$ID, rownames(geno_numeric)), ]
  
  #----------------------------
  # 5. Run rrBLUP!
  #----------------------------
  kinship = A.mat(geno_numeric)
  
  model_alk_w = kin.blup(data=phenotype, geno="ID", pheno="Alkaloids_Res", K=kinship)
  model_mass_w = kin.blup(data=phenotype, geno="ID", pheno="Delta_CT_adj_Res", K=kinship)
  
  heritability_alk_w <- model_alk_w$Vg / (model_alk_w$Vg + model_alk_w$Ve)
  heritability_mass_w <- model_mass_w$Vg / (model_mass_w$Vg + model_mass_w$Ve)
  
  # Return heritability results
  return(list(
    alkaloid_heritability = heritability_alk_w,
    biomass_heritability = heritability_mass_w
  ))
}

# Example usage of the function
hapmap_file <- paste0(data_folder,"/hapmap_files/half_sib_Genos_hapmap.hmp") # Change to your hapmap file path
phenotype_file <- paste0(data_folder,"/Phenotype_Data/3a_data_2024_noOut.txt") # Change to your phenotype file path

result <- calculate_heritability(hapmap_file, phenotype_file)

# Print results
print("Alkaloid heritability:")
print(result$alkaloid_heritability)

print("Biomass heritability:")
print(result$biomass_heritability)



