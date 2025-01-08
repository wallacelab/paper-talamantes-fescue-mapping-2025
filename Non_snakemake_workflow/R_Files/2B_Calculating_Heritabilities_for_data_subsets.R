# Objective: This Script goes after Making_all_Residual_data_subsets.R
# It will find the heritabilities remove the data point furthest from the mean
# and find heritability again. Itll graph this.

library(tidyverse)
library(multcompView)
library(gvlma)
library(ggpubr)
library(car)
library(reshape2)
library(gridExtra)
library(grid)  
library(vcfR)
library(AGHmatrix)
library(sommer)
library(data.table)
library(pheatmap)



################################################################################
# Calculating heritability
################################################################################
base_path <- "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/"

# Load in all 9 datasets
Residual_data_avg_outliars_rm_314x310 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x310.txt", header = TRUE)
Residual_data_avg_outliars_rm_314x312 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm_314x312.txt", header = TRUE)
Residual_data_avg_outliars_rm <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_data_avg_outliars_rm.txt", header = TRUE)

Residual_Data_23_outliars_rm_314x310 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x310.txt", header = TRUE)
Residual_Data_23_outliars_rm_314x312 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm_314x312.txt", header = TRUE)
Residual_Data_23_outliars_rm <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_23_outliars_rm.txt", header = TRUE)

Residual_Data_24_outliars_rm_314x310 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x310.txt", header = TRUE)
Residual_Data_24_outliars_rm_314x312 <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm_314x312.txt", header = TRUE)
Residual_Data_24_outliars_rm <- read.table("../../Data/Phenotype_Data/Residual_Data/Residual_Data_24_outliars_rm.txt", header = TRUE)

# Load the VCF file
vcf_data_loc <- "../../Data/VCF/all_snps_filtered_2.recode.vcf"
vcf_data <- read.vcfR(vcf_data_loc)
geno_matrix <- extract.gt(vcf_data, element = "GT")  
geno_matrix <- t(geno_matrix)

# Load in Cross Identifier file
cross_list_loc <- ("../../Data/Lists/Parent_Progeny_Lists/314_Star_Cross_Parents.txt")
cross_list <- read.table(cross_list_loc, header = FALSE)


############### Done Loading in Data ###########################################
################ Functions needed for heritability finding ##################### 
##### Function to convert genotype into something usable ####
convert_genotypes <- function(geno) {
  geno <- gsub("0/0", "0", geno)
  geno <- gsub("0/1", "1", geno)
  geno <- gsub("1/1", "2", geno)
  return(as.numeric(geno))
}

###### Getting Heritabilities function ###### 
find_heritability <- function(pheno_data, geno_matrix, trait) {
  
  # Ensure both datasets have the same IDs
  common_IDs <- intersect(pheno_data$ID, rownames(geno_matrix))
  
  # Subset data to include only common IDs
  pheno_data <- pheno_data[pheno_data$ID %in% common_IDs, ]
  geno_matrix <- geno_matrix[common_IDs, ]
  
  # Check if row names match
  if (!all(rownames(geno_matrix) == pheno_data$ID)) {
    stop("Row names of geno_data do not match the ID in pheno_data.")
  }
  
  # Convert geno_matrix to data frame, and then convert genotypes
  geno_data <- as.data.frame(geno_matrix)
  geno_data <- geno_data %>% mutate_all(convert_genotypes)
  
  # Convert geno_data back to matrix for kinship matrix calculation
  geno_data <- as.matrix(geno_data)
  
  # Create kinship matrix using Gmatrix function
  kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")
  
  # Ensure the data is in correct format (pheno_data should have necessary columns)
  pheno_data$ID <- as.character(pheno_data$ID)
  
  # Example mixed model for heritability (adjust trait to be a column name)
  model <- mmer(fixed = as.formula(paste(trait, "~ 1")),
                random = ~ vsr(ID, Gu = kinship_matrix),
                rcov = ~ units,
                data = pheno_data)
  
  # Summarize model results
  model_summary <- summary(model)
  varcomp_summary <- model_summary$varcomp
  cat("\n ###### check for convergance ######\n")
  if (model$convergence) {
    print("Model converged successfully.")
  } else {
    print("Model failed to converge.")
  }
  # Extract variance components for heritability calculation
  var_ID <- varcomp_summary$VarComp[1]  # Genetic variance
  var_residual <- varcomp_summary$VarComp[2]  # Residual variance
  
  # Calculate heritability
  heritability <- var_ID / (var_ID + var_residual)
  
  # Return heritability estimate
  heritability <- round(heritability, 4)
  resultstable <- data.table(
    H = heritability,
    GeneticVar = var_ID,
    ResidualVar = var_residual,
    N = model_summary$groups[1]
  )
  return(resultstable)
}

############################## End Function #################################### 
########################### 0 heritability fixer ###############################
# Creates a function that detects 0 heritability datasets, deletes more outliars
# Sends it back to recalculate heritability.

refine_heritability <- function(pheno_data, geno_matrix, trait, max_removals = 5, Hcap = .65, Hmin = .03) {
  removals <- 0
  removed_individuals <- character()  # Initialize an empty vector to store removed IDs
  
  while (removals < max_removals) {
    results <- find_heritability(pheno_data, geno_matrix, trait)
    heritability <- results$H[1]  # Extract heritability
    
    if (heritability >= Hmin && heritability <= Hcap) {
      cat("Heritability is ",heritability," which is acceptable. Returning results.\n")
      # Add removed individuals to the results
      results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
      return(results)
      
    }
    
    # Remove the furthest point from the mean
    cat("Heritability is ",heritability," which is unacceptable\n")
    mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    furthest_index <- which.max(abs_diff)
    
    cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    
    removals <- removals + 1
  }
  # Add removed individuals to the results
  cat("Heritability could not get into acceptable range, returning last results table\n")
  results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
  return(results)
}
########################## End FUnction ########################################
################# graphing heritability as individuals get removed #############
# The function will remove a single individual, graph the heritability on the y
# N on the x. itll do this until ...

graph_heritability <- function(pheno_data, geno_matrix, trait, max_removals = 14, Hdif = .04) {
  removals <- 0
  removed_individuals <- character()  # Initialize an empty vector to store removed IDs
  heritability_old = 100
  N_list <- c()
  H_list <- c()
  
  
  while (removals < max_removals) {
    
    results <- find_heritability(pheno_data, geno_matrix, trait)
    heritability <- results$H[1]  # Extract heritability
    heritability_dif <- abs(heritability - heritability_old)
    # Append the current N and H values
    N_list <- c(N_list, nrow(pheno_data))
    H_list <- c(H_list, heritability)
    cat("\n ######################## \n")
    cat("N_list new value ", N_list, "\n")
    cat("H_list new value ", H_list, "\n")
    
    # if (heritability_dif <= Hdif && heritability != 0) {
    #   cat("Heritability difference is ",heritability_dif," which is acceptable. Returning results.\n")
    #   # Add removed individuals to the results
    #   results$RemovedIndividuals <- paste(removed_individuals, collapse = ", ")
    #   
    #   # Remove and calculate 2 more before finishing.
    #   cat("Removing one more before finishing \n")
    #   mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    #   abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    #   furthest_index <- which.max(abs_diff)
    #   cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    #   removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    #   pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    #   heritability_old <- results$H[1]  
    #   removals <- removals + 1
    #   results <- find_heritability(pheno_data, geno_matrix, trait)
    #   heritability <- results$H[1]  # Extract heritability
    #   heritability_dif <- abs(heritability - heritability_old)
    #   # Append the current N and H values
    #   N_list <- c(N_list, nrow(pheno_data))
    #   H_list <- c(H_list, heritability)
    #   
    #   # Graphing the final product
    #   data <- data.frame(N = N_list, H = H_list)
    #   H_plot <- ggplot(data, aes(x = N, y = H)) +
    #     geom_point(color = "blue", size = 3) +  # Add points
    #     geom_line(color = "red") +             # Add a line connecting the points
    #     labs(title = "N vs Heritability",
    #          x = "N",
    #          y = "Heritability") +
    #     theme_bw() 
    #   
    #   return(list(graph = H_plot, Hdata = data, removed = removed_individuals))
    # 
    # }
    
    # Remove the furthest point from the mean
    cat("Heritability difference is ",heritability_dif," which is unacceptable\n")
    mean_trait <- mean(pheno_data[[trait]], na.rm = TRUE)
    abs_diff <- abs(pheno_data[[trait]] - mean_trait)
    furthest_index <- which.max(abs_diff)
    
    cat("Removing index:", pheno_data$ID[furthest_index], "with value:", pheno_data[[trait]][furthest_index], "\n")
    removed_individuals <- c(removed_individuals, pheno_data$ID[furthest_index])  # Store removed individual
    pheno_data <- pheno_data[-furthest_index, ]  # Remove the individual
    heritability_old <- results$H[1]      # Update old heritability
    removals <- removals + 1
  }
  
  # Add removed individuals to the results
  # Graphing the final product
  data <- data.frame(N = N_list, H = H_list)
  H_plot <- ggplot(data, aes(x = N, y = H)) +
    geom_point(color = "blue", size = 3) +  # Add points
    geom_line(color = "red") +             # Add a line connecting the points
    labs(title = "N vs Heritability",
         x = "N",
         y = "Heritability") +
    theme_bw() 
  
  return(list(graph = H_plot, Hdata = data, removed = removed_individuals))
}
########################## End FUnction ########################################
# Calculating Heritability from all data sets

# # Both years avaraged, heritabilities 
# stats_avg_310_alk <- refine_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res")
# H_avg_310_alk <- stats_avg_310_alk$H[1]
# 
# stats_avg_312_alk <- refine_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res")
# H_avg_312_alk <- stats_avg_312_alk$H[1]
# 
# stats_avg_star_alk <- refine_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Alkaloids_Res")
# H_avg_star_alk <- stats_avg_star_alk$H[1]
# 
# stats_avg_310_ct <- refine_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_310_ct <- stats_avg_310_ct$H[1]
# 
# stats_avg_312_ct <- refine_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_312_ct <- stats_avg_312_ct$H[1]
# 
# stats_avg_star_ct <- refine_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res")
# H_avg_star_ct <- stats_avg_star_ct$H[1]
# 
# 
# # Just the 2023 heritabilities
# stats_2023_310_alk <- refine_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res")
# H_2023_310_alk <- stats_2023_310_alk$H[1]
# 
# stats_2023_312_alk <- refine_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res")
# H_2023_312_alk <- stats_2023_312_alk$H[1]
# 
# stats_2023_star_alk <- refine_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Alkaloids_Res")
# H_2023_star_alk <- stats_2023_star_alk$H[1]
# 
# stats_2023_310_ct <- refine_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res", 15, .6)
# H_2023_310_ct <- stats_2023_310_ct$H[1]
# 
# stats_2023_312_ct <- refine_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2023_312_ct <- stats_2023_312_ct$H[1]
# 
# stats_2023_star_ct <- refine_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res", 15, .60)
# H_2023_star_ct <- stats_2023_star_ct$H[1]
# 
# # Just the 2024 Heritabilities
# stats_2024_310_alk <- refine_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res", 15, .60)
# H_2024_310_alk <- stats_2024_310_alk$H[1]
# 
# stats_2024_312_alk <- refine_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res", 20, .60)
# H_2024_312_alk <- stats_2024_312_alk$H[1]
# 
# stats_2024_star_alk <- refine_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Alkaloids_Res", 20, .60)
# H_2024_star_alk <- stats_2024_star_alk$H[1]
# 
# stats_2024_310_ct <- refine_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res", 5, .60)
# H_2024_310_ct <- stats_2024_310_ct$H[1]
# 
# stats_2024_312_ct <- refine_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2024_312_ct <- stats_2024_312_ct$H[1]
# 
# stats_2024_star_ct <- refine_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res")
# H_2024_star_ct <- stats_2024_star_ct$H[1]


########## Plotting Heritability platoe ################

plot_avg_310_alk <- graph_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm_314x310)/2, .04)
plot_avg_312_alk <- graph_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm_314x312)/2, .04)
plot_avg_star_alk <- graph_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Alkaloids_Res", nrow(Residual_data_avg_outliars_rm)/2)
plot_avg_310_ct <- graph_heritability(Residual_data_avg_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm_314x310)/2)
plot_avg_312_ct <- graph_heritability(Residual_data_avg_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm_314x312)/2)
plot_avg_star_ct <- graph_heritability(Residual_data_avg_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_data_avg_outliars_rm)/2)

plot_2023_310_alk <- graph_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm_314x310)/2)
plot_2023_312_alk <- graph_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm_314x312)/2)
plot_2023_star_alk <- graph_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_23_outliars_rm)/2)
plot_2023_310_ct <-graph_heritability(Residual_Data_23_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm_314x310)/2)
plot_2023_312_ct <- graph_heritability(Residual_Data_23_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm_314x312)/2)
plot_2023_star_ct <- graph_heritability(Residual_Data_23_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_23_outliars_rm)/2)

plot_2024_310_alk <- graph_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm_314x310)/2) ## This one is all 0's
plot_2024_312_alk <- graph_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm_314x312)/2)
plot_2024_star_alk <- graph_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Alkaloids_Res",nrow(Residual_Data_24_outliars_rm)/2)
plot_2024_310_ct <- graph_heritability(Residual_Data_24_outliars_rm_314x310, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_24_outliars_rm_314x310)/2)
plot_2024_312_ct <- graph_heritability(Residual_Data_24_outliars_rm_314x312, geno_matrix, trait = "Delta_CT_adj_Res", nrow(Residual_Data_24_outliars_rm_314x312)/2)
plot_2024_star_ct <- graph_heritability(Residual_Data_24_outliars_rm, geno_matrix, trait = "Delta_CT_adj_Res",nrow(Residual_Data_24_outliars_rm)/2)




# Plot for avg data
# Create overall titles
title_avg <- textGrob("Average Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_avg_310_alk$graph, plot_avg_312_alk$graph, plot_avg_star_alk$graph,
              plot_avg_310_ct$graph,plot_avg_312_ct$graph,plot_avg_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)

# Plot for 2023 data
# Create overall titles
title_avg <- textGrob("2023 Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_2023_310_alk$graph, plot_2023_312_alk$graph, plot_2023_star_alk$graph,
              plot_2023_310_ct$graph,plot_2023_312_ct$graph,plot_2023_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)


# Plot for 2024 data
# Create overall titles
title_avg <- textGrob("2024 Data", gp = gpar(fontsize = 20, fontface = "bold"))
x_axis <- textGrob("310 \t\t\t\t\t 312 \t\t\t\t\t Star", gp = gpar(fontsize = 16), vjust = -.5)
y_axis <- textGrob("CT Ratio \t\t Alkaloids", gp = gpar(fontsize = 16), rot = 90)
# Arrange plots with grid.arrange
avg_H_plot <- grid.arrange(
  arrangeGrob(plot_2024_310_alk$graph, plot_2024_312_alk$graph, plot_2024_star_alk$graph,
              plot_2024_310_ct$graph,plot_2024_312_ct$graph,plot_2024_star_ct$graph, 
              ncol = 3, nrow = 2),  # Add your plots here
  top = title_avg,                             # Add the top title
  bottom = x_axis,                             # Add the x-axis title
  left = y_axis                                # Add the y-axis title
)


####################### Saving the data ########################################
output_dir <- "../../Data/Heritability_Outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to save plots and data tables
save_results <- function(result, name, output_dir) {
  # Save the plot (assuming the first element is the plot)
  plot_file <- file.path(output_dir, paste0(name, ".png"))
  ggsave(plot_file, result[[1]])
  
  # Save the data tables (assuming the second and third elements are data tables)
  if (length(result) > 1) {
    for (i in 2:length(result)) {
      data_file <- file.path(output_dir, paste0(name, "_data_", i - 1, ".csv"))
      write.csv(result[[i]], data_file, row.names = FALSE)
    }
  }
}

# Save all results
save_results(plot_avg_310_alk, "plot_avg_310_alk_gen_removed", output_dir)
save_results(plot_avg_312_alk, "plot_avg_312_alk_gen_removed", output_dir)
save_results(plot_avg_star_alk, "plot_avg_star_alk_gen_removed", output_dir)
save_results(plot_avg_310_ct, "plot_avg_310_ct_gen_removed", output_dir)
save_results(plot_avg_312_ct, "plot_avg_312_ct_gen_removed", output_dir)
save_results(plot_avg_star_ct, "plot_avg_star_ct_gen_removed", output_dir)

save_results(plot_2023_310_alk, "plot_2023_310_alk_gen_removed", output_dir)
save_results(plot_2023_312_alk, "plot_2023_312_alk_gen_removed", output_dir)
save_results(plot_2023_star_alk, "plot_2023_star_alk_gen_removed", output_dir)
save_results(plot_2023_310_ct, "plot_2023_310_ct_gen_removed", output_dir)
save_results(plot_2023_312_ct, "plot_2023_312_ct_gen_removed", output_dir)
save_results(plot_2023_star_ct, "plot_2023_star_ct_gen_removed", output_dir)

save_results(plot_2024_310_alk, "plot_2024_310_alk_gen_removed", output_dir)
save_results(plot_2024_312_alk, "plot_2024_312_alk_gen_removed", output_dir)
save_results(plot_2024_star_alk, "plot_2024_star_alk_gen_removed", output_dir)
save_results(plot_2024_310_ct, "plot_2024_310_ct_gen_removed", output_dir)
save_results(plot_2024_312_ct, "plot_2024_312_ct_gen_removed", output_dir)
save_results(plot_2024_star_ct, "plot_2024_star_ct_gen_removed", output_dir)



################# Genetic relatedness matrix ###################################

################## Heatmap and PCA Analysis Function ###########################
base_kin_analysis <- function(pheno_data, geno_matrix, cross_list, title, removed_data_alk, removed_data_ct) {
    
  # Ensure both datasets have the same IDs
  common_IDs <- intersect(pheno_data$ID, rownames(geno_matrix))
  # Subset data to include only common IDs
  pheno_data <- pheno_data[pheno_data$ID %in% common_IDs, ]
  geno_matrix <- geno_matrix[common_IDs, ]
  # Check if row names match
  if (!all(rownames(geno_matrix) == pheno_data$ID)) {
    stop("Row names of geno_data do not match the ID in pheno_data.")
  }
  # Convert geno_matrix to data frame, and then convert genotypes
  geno_data <- as.data.frame(geno_matrix)
  geno_data <- geno_data %>% mutate_all(convert_genotypes)
  # Convert geno_data back to matrix for kinship matrix calculation
  geno_data <- as.matrix(geno_data)
  # Create kinship matrix using Gmatrix function
  kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")
  pca_result <- prcomp(kinship_matrix)
  diag(kinship_matrix) <- NA
  
  kin_heatmap <- pheatmap(kinship_matrix, 
               cluster_rows = TRUE, 
               cluster_cols = TRUE, 
               color = colorRampPalette(c("blue", "white", "red"))(50),
               main = paste0("Kinship Matrix Heatmap ", title))
  
  
  # Perform PCA on the kinship matrix (or SNP matrix if needed)
  pca_scores <- pca_result$x  # Principal component scores
  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pca_df <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2])
  pca_df <-merge(pca_df, cross_list, by.x = "row.names", by.y = "V1" )
  pca_df <- pca_df %>% rename(Cross = V2)
  
  # Create the scatter plot of PC1 vs PC2
  kin_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
            geom_point(aes(color = Cross), size = 3, alpha = 0.3) +  # Plot the points
            geom_text(
              data = subset(pca_df, grepl("parent", Cross, ignore.case = TRUE)),  # Filter for "parent"
              aes(label = Row.names),
              vjust = -1,  # Adjust text position
              size = 3,
              color = "black"
            ) +
            labs(title = paste0("PCA:", title), x = "PC1", y = "PC2") +
            theme_minimal()
  
  # PCA of data being removed
  data_path <- "../../Data/Heritability_Outputs/"
  number_heri_alk <- read.csv(paste0(data_path, removed_data_alk, "_data_1.csv"), header = TRUE)
  indviduals_removed_alk <- read.csv(paste0(data_path, removed_data_alk, "_data_2.csv"), header = TRUE)
  number_heri_ct <- read.csv(paste0(data_path, removed_data_ct, "_data_1.csv"), header = TRUE)
  indviduals_removed_ct <- read.csv(paste0(data_path, removed_data_ct, "_data_2.csv"), header = TRUE)
  
  NHI_data_alk <- cbind(number_heri_alk, indviduals_removed_alk)
  NHI_data_ct <- cbind(number_heri_ct, indviduals_removed_ct)
  
  
  pca_df_2 <- merge(pca_df,NHI_data_alk, by.x = "Row.names", by.y = "x", all = TRUE )
  pca_df_2$color <- ifelse(is.na(pca_df_2$N), "black", "purple")
  
  pca_df_3 <- merge(pca_df,NHI_data_ct, by.x = "Row.names", by.y = "x", all = TRUE )
  pca_df_3$color <- ifelse(is.na(pca_df_3$N), "black", "purple")
  
  #Create plot that shows removed individuals
  removed_data_pca <- ggplot(pca_df_2, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = color), size = 3) + # Color based on N
    scale_color_manual(values = c("black", "purple"),
                       labels = c("black" = "Stayed", "purple" = "Removed")) + # Set color scheme
    geom_text(aes(label = ifelse(!is.na(N), N, "")), vjust = -1, size = 3) + # Add N value above points
    labs(title = paste0("PCA Plot Colored by Removed Individuals ", "Alkaloids" ), x = "PC1", y = "PC2") +
    theme_minimal() 
  
  removed_data_pca2 <- ggplot(pca_df_3, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = color), size = 3) + # Color based on N
    scale_color_manual(values = c("black", "purple"),
                       labels = c("black" = "Stayed", "purple" = "Removed")) + # Set color scheme
    geom_text(aes(label = ifelse(!is.na(N), N, "")), vjust = -1, size = 3) + # Add N value above points
    labs(title = paste0("PCA Plot Colored by Removed Individuals ", "CT" ), x = "PC1", y = "PC2") +
    theme_minimal() 
  
  # PCA graphic showing individuals removed from both alkaloids and CT
  # Combine the IDs from both datasets
  combined_ids <- data.frame(ID = unique(c(NHI_data_alk$x, NHI_data_ct$x)))
  
  # Check where each ID appears
  combined_ids$Source <- ifelse(
    combined_ids$ID %in% NHI_data_alk$x & combined_ids$ID %in% NHI_data_ct$x, "Both",
    ifelse(combined_ids$ID %in% NHI_data_alk$x, "NHI_data_alk", "NHI_data_ct")
  )
  
  pca_df_4 <- merge(pca_df,combined_ids, by.x = "Row.names", by.y = "ID", all = TRUE )
  pca_df_4$Source[is.na(pca_df_4$Source)] <- "Not_Removed"
  
  
  
  #Create plot that shows removed individuals
  removed_data_by_pheno <- ggplot(pca_df_4, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Source), size = 4, alpha = 0.75) + # Color based on N
    scale_color_manual(values = c("purple", "red", "blue", "black"),
                       labels = c("Removed in Both", "Removed in Alkaloids", "Removed in CT", "Not Removed")) + 
    labs(title = paste0("PCA Plot Colored by Removed Individuals by Phenotype" ), x = "PC1", y = "PC2") +
    theme_bw()
  
  
  # return plot list
  return(list(plot1 = kin_heatmap, plot2 = kin_pca, plot3 = removed_data_pca, plot4 = removed_data_pca2, plot5 = removed_data_by_pheno))

}
##################### Function End #############################################

nrow(Residual_Data_23_outliars_rm)
plots_star_23 <- base_kin_analysis(Residual_Data_23_outliars_rm, geno_matrix, cross_list, title = "Star Cross 2023", removed_data_alk = "plot_2023_star_alk", removed_data_ct = "plot_2023_star_ct")
plots_314x310_23 <- base_kin_analysis(Residual_Data_23_outliars_rm_314x310, geno_matrix, cross_list, title = "314x310 2023", removed_data_alk = "plot_2023_310_alk", removed_data_ct = "plot_2023_310_ct")
plots_314x312_23 <- base_kin_analysis(Residual_Data_23_outliars_rm_314x312, geno_matrix, cross_list, title = "314x312 2023",removed_data_alk = "plot_2023_312_alk", removed_data_ct = "plot_2023_312_ct")

nrow(Residual_Data_24_outliars_rm)
plots_star_24 <- base_kin_analysis(Residual_Data_24_outliars_rm, geno_matrix, cross_list, title = "Star Cross 2024", removed_data_alk = "plot_2024_star_alk", removed_data_ct = "plot_2024_star_ct")
plots_314x310_24 <- base_kin_analysis(Residual_Data_24_outliars_rm_314x310, geno_matrix, cross_list, title = "314x310 2024", removed_data_alk = "plot_2024_310_alk", removed_data_ct = "plot_2024_310_ct")
plots_314x312_24 <- base_kin_analysis(Residual_Data_24_outliars_rm_314x312, geno_matrix, cross_list, title = "314x312 2024",removed_data_alk = "plot_2024_312_alk", removed_data_ct = "plot_2024_312_ct")

nrow(Residual_data_avg_outliars_rm)
plots_star_avg <- base_kin_analysis(Residual_data_avg_outliars_rm, geno_matrix, cross_list, title = "Star Cross Years Avraged", removed_data_alk = "plot_avg_star_alk", removed_data_ct = "plot_avg_star_ct")
plots_314x310_avg <- base_kin_analysis(Residual_data_avg_outliars_rm_314x310, geno_matrix, cross_list, title = "314x310 Years Avraged", removed_data_alk = "plot_avg_310_alk", removed_data_ct = "plot_avg_310_ct")
plots_314x312_avg <- base_kin_analysis(Residual_data_avg_outliars_rm_314x312, geno_matrix, cross_list, title = "314x312 Years Avraged", removed_data_alk = "plot_avg_312_alk", removed_data_ct = "plot_avg_312_ct")

plot1 <- plots_star_avg$plot1
plot2 <- plots_star_avg$plot2
plot3 <- plots_star_avg$plot3
plot4 <- plots_star_avg$plot4
plot5 <- plots_star_avg$plot5


grid.arrange(plot2, plot3, plot4, plot5, ncol = 2, nrow = 2) 




# ../../Data/Heritability_Outputs/plot_2023_310_alk_data_1.csv

##################### PCA Graph Iterator #######################################
# This function will color the points that get removed in the PCA and number em




##################### Function End #############################################













































#################### fixing functions         ##################################

pheno_data <- Residual_data_avg_outliars_rm
geno_matrix2 <- geno_matrix
cross_list2 <- cross_list
title <- "Star Cross Years Avraged"
removed_data_alk <- "plot_avg_star_alk" 
removed_data_ct <- "plot_avg_star_ct"

# Ensure both datasets have the same IDs
common_IDs <- intersect(pheno_data$ID, rownames(geno_matrix2))
# Subset data to include only common IDs
pheno_data <- pheno_data[pheno_data$ID %in% common_IDs, ]
geno_matrix2 <- geno_matrix2[common_IDs, ]
# Check if row names match
if (!all(rownames(geno_matrix2) == pheno_data$ID)) {
  stop("Row names of geno_data do not match the ID in pheno_data.")
}
# Convert geno_matrix2 to data frame, and then convert genotypes
geno_data <- as.data.frame(geno_matrix2)
geno_data <- geno_data %>% mutate_all(convert_genotypes)
# Convert geno_data back to matrix for kinship matrix calculation
geno_data <- as.matrix(geno_data)
# Create kinship matrix using Gmatrix function
kinship_matrix <- Gmatrix(SNPmatrix = geno_data, method = "VanRaden")
pca_result <- prcomp(kinship_matrix)
diag(kinship_matrix) <- NA

kin_heatmap <- pheatmap(kinship_matrix, 
                        cluster_rows = TRUE, 
                        cluster_cols = TRUE, 
                        color = colorRampPalette(c("blue", "white", "red"))(50),
                        main = paste0("Kinship Matrix Heatmap ", title))


# Perform PCA on the kinship matrix (or SNP matrix if needed)
pca_scores <- pca_result$x  # Principal component scores
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pca_df <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2])
pca_df <-merge(pca_df, cross_list2, by.x = "row.names", by.y = "V1" )
pca_df <- pca_df %>% rename(Cross = V2)

# Create the scatter plot of PC1 vs PC2
kin_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Cross), size = 3, alpha = 0.3) +  # Plot the points
  geom_text(
    data = subset(pca_df, grepl("parent", Cross, ignore.case = TRUE)),  # Filter for "parent"
    aes(label = Row.names),
    vjust = -1,  # Adjust text position
    size = 3,
    color = "black"
  ) +
  labs(title = paste0("PCA:", title), x = "PC1", y = "PC2") +
  theme_minimal()

# PCA of data being removed
data_path <- "../../Data/Heritability_Outputs/"
number_heri_alk <- read.csv(paste0(data_path, removed_data_alk, "_data_1.csv"), header = TRUE)
indviduals_removed_alk <- read.csv(paste0(data_path, removed_data_alk, "_data_2.csv"), header = TRUE)
number_heri_ct <- read.csv(paste0(data_path, removed_data_ct, "_data_1.csv"), header = TRUE)
indviduals_removed_ct <- read.csv(paste0(data_path, removed_data_ct, "_data_2.csv"), header = TRUE)

NHI_data_alk <- cbind(number_heri_alk, indviduals_removed_alk)
NHI_data_ct <- cbind(number_heri_ct, indviduals_removed_ct)


pca_df_2 <- merge(pca_df,NHI_data_alk, by.x = "Row.names", by.y = "x", all = TRUE )
pca_df_2$color <- ifelse(is.na(pca_df_2$N), "black", "purple")

#Create plot that shows removed individuals
removed_data_pca <- ggplot(pca_df_2, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = color), size = 3) + # Color based on N
  scale_color_manual(values = c("black", "purple"),
                     labels = c("black" = "Stayed", "purple" = "Removed")) + # Set color scheme
  geom_text(aes(label = ifelse(!is.na(N), N, "")), vjust = -1, size = 3) + # Add N value above points
  labs(title = paste0("PCA Plot Colored by Removed Individuals ", "Alkaloids" ), x = "PC1", y = "PC2") +
  theme_minimal() 

# PCA graphic showing individuals removed from both alkaloids and CT
# Combine the IDs from both datasets
combined_ids <- data.frame(ID = unique(c(NHI_data_alk$x, NHI_data_ct$x)))

# Check where each ID appears
combined_ids$Source <- ifelse(
  combined_ids$ID %in% NHI_data_alk$x & combined_ids$ID %in% NHI_data_ct$x, "Both",
  ifelse(combined_ids$ID %in% NHI_data_alk$x, "NHI_data_alk", "NHI_data_ct")
)

pca_df_4 <- merge(pca_df,combined_ids, by.x = "Row.names", by.y = "ID", all = TRUE )
pca_df_4$Source[is.na(pca_df_4$Source)] <- "Not_Removed"



#Create plot that shows removed individuals
removed_data_by_pheno <- ggplot(pca_df_4, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Source), size = 4, alpha = 0.75) + # Color based on N
  scale_color_manual(values = c("purple", "red", "blue", "black"),
                     labels = c("Removed in Both", "Removed in Alkaloids", "Removed in CT", "Not Removed")) + 
  labs(title = paste0("PCA Plot Colored by Removed Individuals by Phenotype" ), x = "PC1", y = "PC2") +
  theme_bw() 




# What are the 3 to 4 biological insights that we can confirm based off of this data
# 1 fungal load not correlated to alkaloid amount
# 2 grass genetics seemingly have no influence on fungal laod
# Can try an anova by year or plant ID to see the environmental impact on the data (residuals are assumed to be a normal distribution)
