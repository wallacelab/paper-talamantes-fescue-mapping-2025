library(tidyverse)

# Load in sample
sample_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Mapped_Reads/320-5-9_depth.txt"
sample <- read.table(sample_loc, header = TRUE)

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
  
