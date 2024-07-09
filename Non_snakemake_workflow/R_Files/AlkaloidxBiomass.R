# The point of this Rscript is to combine all the meta data to one common file.
library(ggplot2)

# Loading in data
final_List_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/residual_data_father_included_2.txt"



### Visualize the best crosses
final_List <- read.csv(final_List_loc, header = TRUE, strip.white=TRUE)
final_List <- subset(final_List, Status != "dead")


# linear model of biomass efficiency adjusted x alkaloid
model1 <- lm(filtered_CT_model.residuals ~ alkaloid_model.residuals, data = final_List)
model2 <- lm(alkaloid_model.residuals ~ filtered_CT_model.residuals, data = final_List)

summary(model1)
summary(model2)

ggplot(final_List, aes(x = filtered_CT_model.residuals, y = alkaloid_model.residuals)) +
  geom_point() + 
  labs(x = "Adjusted Biomass Risidual Values", y = "Alkaloid Risiduals", title = "Alklaoid vs Biomass")  # Add labels and title




final_List$filtered_CT_model.residuals
final_List$alkaloid_model.residuals



CT_model_raw <- lm(DeltaCT_Raw ~ Harvest_Date + Standard + Extraction_Date + Extractor + Data_Set, data = Delta_CT_Data_raw)





ggplot(cross_count_subseted,aes( x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Alive Plants", x = "Crosses", y = "Count") +
  theme(axis.text.x = element_text(size = 12)) +
  geom_text(aes(label = Freq), vjust = -0.5)


