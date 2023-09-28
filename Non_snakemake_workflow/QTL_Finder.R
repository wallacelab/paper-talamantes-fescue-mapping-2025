library(tidyverse)

# Load in sample
snps_loc = "/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/Done/Variants_progeny.bcf"
snps <- read.table(snps_loc, header = TRUE)

#Get counts for sample
counts <- sample %>%
  count(X1)

# Remove lower depths
counts <- counts[-c(1,2,3),]


#Display on histogram 
ggplot(data=counts, aes(x=X1, y=n)) + 
  geom_bar(stat = "identity") +
  xlim(0,100) +
  ylim(0,20000)
  
