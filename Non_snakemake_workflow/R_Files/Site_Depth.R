#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a file path was provided
if (length(args) == 0) {
  stop("Usage: calculate_avg_depth.R <path_to_depth_file>")
}

# Read the file path
file_path <- args[1]

# Check if the file exists
if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

# Load the SNP depth file into a data frame
snp_depth <- read.table(file_path, header = FALSE)

# Select the depth columns (ignoring the first two columns)
depth_data <- snp_depth[, -c(1:2)]

# Calculate the average depth per site (row)
row_averages <- rowMeans(depth_data)

# Calculate the overall average depth across all sites
overall_average <- mean(row_averages)

# Print the result
cat("The overall average site depth is:", overall_average, "\n")



