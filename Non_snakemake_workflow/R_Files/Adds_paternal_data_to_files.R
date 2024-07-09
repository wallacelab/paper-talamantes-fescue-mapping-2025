# The point of this Rscript is to combine all the meta data to one common file.
library(ggplot2)

# Loading in data
all_residual_data_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/residual_data.txt"
father_data_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Meta_Data/Mother_Father_Data.csv" 
live_list_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/Live_List.csv"
save_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/residual_data_father_included.txt"
final_List_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/residual_data_father_included_2.txt"


all_residual_data <- read.table(all_residual_data_loc, header = TRUE, strip.white=TRUE)
father_data <- read.csv(father_data_loc, header = TRUE, strip.white=TRUE)
live_list <- read.csv(live_list_loc, header = TRUE, strip.white=TRUE)

father_data <- subset(father_data, select = -c(Mother))
merged_df <- merge(all_residual_data, father_data, by = "ID", all = TRUE)
merged_df <- merge(live_list, merged_df, by = "ID", all = TRUE)
merged_df$Status[is.na(merged_df$Status)] <- "dead"

write.table(merged_df, file = save_loc, sep = "\t", row.names = FALSE)

# Extract a certain cross
single_cross <- merged_df[merged_df$Status == "Alive" & (merged_df$Maternal_Parent %in% c("319", "320")) & (merged_df$Father %in% c("319", "320")), ]
write.csv(single_cross, "/home/drt06/Downloads/Plants2Keep/320x319")


### Visualize the best crosses
final_List <- read.csv(final_List_loc, header = TRUE, strip.white=TRUE)
final_List <- subset(final_List, Status != "dead")
# final_List <- subset(final_List, Concentration != "#N/A")
final_List$Combo <- paste(final_List$Maternal_Parent, final_List$Father, sep = "-")
cross_count <- data.frame(table(final_List$Combo))
cross_count_subseted <- cross_count[cross_count$Freq >= 10, ]

ggplot(cross_count_subseted,aes( x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Alive Plants", x = "Crosses", y = "Count") +
  theme(axis.text.x = element_text(size = 12)) +
  geom_text(aes(label = Freq), vjust = -0.5)


