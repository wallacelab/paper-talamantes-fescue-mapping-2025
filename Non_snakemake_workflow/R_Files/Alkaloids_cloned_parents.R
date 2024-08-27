# This script will look at the alkaloid data for the cloned tested parents
library(ggplot2)
library(tidyverse)

# Loading in data
parent_alkaloid_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Phenotype_Data/Alkaloids_parent_testing.csv"



### Visualize the best crosses
parent_alkaloid <- read.csv(parent_alkaloid_loc, header = FALSE, strip.white=TRUE)

################################################################################
# cleaning up the data
################################################################################
parent_alkaloid$V4 <- as.character(parent_alkaloid$V4)
parent_alkaloid <- parent_alkaloid %>% rename(Replicates = V2)
parent_alkaloid <- parent_alkaloid %>% rename(Replication_Status = V6)


create_anova_table <-function(parent){
df_subset <- subset(parent_alkaloid, V4 == parent, select = c(Replicates, V4, V5))
anova_model <- aov(V5 ~ Replicates, data = df_subset)
anova_summary <- summary(anova_model) # f stat is The ratio of between-group variance to within-group variance
f_value <- anova_summary[[1]]["group", "F value"]
p_value <- anova_summary[[1]]["group", "Pr(>F)"]

result_table <- data.frame(
  Value = parent,
  F_value = f_value,
  P_value = p_value
)
return(result_table)
}


parent_list <- unique(parent_alkaloid$V4)

# For loop to cycle through the parent_list and run the create_anova_table function
result_list <- list()  # Initialize an empty list to store results
for (value in parent_list) {
  result <- create_anova_table(value)  # Run the function
  result_list[[as.character(value)]] <- result  # Store the result in a named list
}

################################################################################
# Making Graphs
################################################################################
ggplot(parent_alkaloid, aes(x = V4, y = V5, color = Replication_Status)) +
  scale_color_manual(values = c("violet","navy","seagreen","goldenrod", "black")) + 
  geom_point(size = 4, alpha = .6, position = position_jitter(width = 0.1)) +  # Use geom_point for scatter plot
  labs(title = "Alkaloid Levels in Parent Replicates", 
       x = "Parents", 
       y = "ng/g ergot alkaloids") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5)
  )


