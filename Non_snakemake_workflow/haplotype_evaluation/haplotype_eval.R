# Objective: This script should load in a vcf file. It will evaluate it and 
# determine where the parental haplotypes are.
library(tidyverse)
library(reshape2) 

# Load in data
prog_315_4_43_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/315-4-43"
parent_comparison_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/315x320"
prog_315_4_43 <- read.csv(prog_315_4_43_loc, header = TRUE, strip.white=TRUE)
parent_comparison <- read.csv(parent_comparison_loc, header = TRUE, strip.white=TRUE)

#Fixing the number labels
prog_315_4_43$X <- prog_315_4_43$X +1
prog_315_4_43$X <- prog_315_4_43$X +1

################################################################################
### Functions ###
################################################################################
#function to detect highest ratio in a row. Ratio must be .04 than second best ratio.
get_high_ratio <- function(row) 
  {
  row_list <- as.list(row)
  numeric_list <- unlist(row_list[sapply(row_list, is.numeric)])
  max_ratio <- sort(numeric_list, decreasing = TRUE)[1]
  second_max_ratio <- sort(numeric_list, decreasing = TRUE)[2]

  if (max_ratio - second_max_ratio > .04){
    max_ratio_index <- which.max(row)
    return(max_ratio_index)
    }
    else {
      return(0) 
    }
}

# use get_high_ratio function to append row with 0s and 1s to a larger data frame
Finding_highest_only <- function(data_set) 
{
  nrow(data_set) # Getting length of dataframe
  for (i in  1:nrow(data_set)){
    zero_row <- integer(4)
    max_index <- get_high_ratio(data_set[i,])
    if (max_index > 0){
      zero_row[max_index] <- 1 
    }
    new_row <- as.data.frame(t(zero_row))
    if (!exists("my_dataframe")) {
      my_dataframe <- new_row
    } else {
      # Append the new row to the existing data frame
      my_dataframe <- rbind(my_dataframe, new_row)
    }
  }
  return(my_dataframe)
}

# This function takes a data set from one parent and the entire data set.
# It cleans the data then puts it into Finding_highest_only and cleans its output
clean_data_Finding_highest_only <- function(parent_data_set,full_data)
{
  parent_data_setX <- parent_data_set[, -c(1,2)]
  highestonly_data_set <- Finding_highest_only(parent_data_setX)
  # Add in x axis numbers
  highestonly_data_set <- cbind(full_data$X, full_data$Chrom, highestonly_data_set)
  colnames(highestonly_data_set)[1] <- "X"
  colnames(highestonly_data_set)[2] <- "Chrom"
  return(highestonly_data_set)
}
################################################################################
### Use functions ###
################################################################################
#split data frame by A and B for progeny to parent comparisons
prog_315_4_43A <- prog_315_4_43[, c(1:2,3:6)]
prog_315_4_43B <- prog_315_4_43[, c(1:2,7:10)]

highestonly_A <- clean_data_Finding_highest_only(prog_315_4_43A, prog_315_4_43 )
highestonly_B <- clean_data_Finding_highest_only(prog_315_4_43B, prog_315_4_43 )
highestonly_A$Chrom <- as.factor(highestonly_A$Chrom)

# Split data frame by A and B for parent to parent comparisons
parent_comparisonA <- parent_comparison[, c(1:2,3:6)]
parent_comparisonB <- parent_comparison[, c(1:2,7:10)]

highestonly_A <- clean_data_Finding_highest_only(parent_comparisonA, parent_comparison )
highestonly_B <- clean_data_Finding_highest_only(parent_comparisonB, parent_comparison )
highestonly_A$Chrom <- as.factor(highestonly_A$Chrom)

################################################################################
### Plot Making ###
################################################################################
# These next sections will allow me to visualize each chromosome
# Get the unique groups
unique_groups <- unique(highestonly_A$Chrom)
# Calculate the x-axis range for each group
group_ranges <- tapply(highestonly_A$X, highestonly_A$Chrom, range)
# Create a data frame for rectangles
rect_df <- data.frame(
  xmin = sapply(group_ranges, function(x) if(length(x) >= 1) x[1] else NA),
  xmax = sapply(group_ranges, function(x) if(length(x) >= 2) x[2] else NA),
  ymin = -Inf,
  ymax = Inf,
  fill = unique_groups)
### Showing multiple chromosomes and displaying the highest match with dots
ggplot(data = highestonly_A, aes(x=X)) +
  geom_rect(data = rect_df, inherit.aes=FALSE, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.5, show.legend = FALSE) +
  geom_point(aes( y=V1), color="cyan", size = 1) +
  geom_point(aes( y=V2), color="green", size = 1) +
  geom_point(aes( y=V3), color="red", size = 1) +
  geom_point(aes( y=V4), color="purple", size = 1) +
  scale_y_continuous(limits = c(.1, 1.1)) +
  scale_x_continuous(limits = c(0,480)) +
  ggtitle("Haplotype") +
  xlab("SNP Number")

###### Line Graphs #####
prog_315_4_43Ac1 <- subset(prog_315_4_43A, subset = Chrom=="ptg000002l")
prog_315_4_43Ac1 <- melt(prog_315_4_43Ac1, id = c("X", "Chrom"))
# Line graph of all haplotype percent matches 
ggplot(data = prog_315_4_43Ac1, aes(x=X, y=value, group=variable, color=variable)) +
  geom_line(size=1) +
  ggtitle("Chromosome 2 ProgenyA to Parents") +
  xlab("Window Number") +
  ylab("Percent")

parent_comparisonAc1 <- subset(parent_comparisonA, subset = Chrom=="ptg000002l")
parent_comparisonAc1 <- melt(parent_comparisonAc1, id = c("X", "Chrom"))
# Line graph of all haplotype percent matches 
ggplot(data = parent_comparisonAc1, aes(x=X, y=value, group=variable, color=variable)) +
  geom_line(size=1) +
  ggtitle("Chromosome 2 Parent1A to Parents") +
  xlab("Window Number") +
  ylab("Percent")

parent_comparisonBc1 <- subset(parent_comparisonB, subset = Chrom=="ptg000002l")
parent_comparisonBc1 <- melt(parent_comparisonBc1, id = c("X", "Chrom"))
# Line graph of all haplotype percent matches 
ggplot(data = parent_comparisonBc1, aes(x=X, y=value, group=variable, color=variable)) +
  geom_line(size=1) +
  ggtitle("Chromosome 2 Parent2A to Parents") +
  xlab("Window Number") +
  ylab("Percent")
  
  
  

