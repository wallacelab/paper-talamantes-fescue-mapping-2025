# The objective of this script is to turn my phenotypes into BLUPs and analyze the 
# genotypes aswell to hopefully get sometype of heritability caluculation out of it. 

# Installing libraries
library(lme4)

# Loading data
pheno_24_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/2024_Data/Final_2024_phenotype_data.csv"
pheno_24 <- read.csv(pheno_24_loc, header = TRUE, strip.white=TRUE)

# Fiding the cross 
pheno_24$Cross <- apply(pheno_24[, c("Mother", "Father")], 1, function(x) paste(sort(x), collapse = "x"))

# Creating models of my ohenotpyes
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














