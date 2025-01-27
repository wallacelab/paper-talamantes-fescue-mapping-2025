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
snp_depth <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)

# Select the depth columns (ignoring the first two columns)
depth_data <- snp_depth[, -c(1:2)]

# Ensure all columns are numeric
depth_data <- as.data.frame(lapply(depth_data, function(col) as.numeric(as.character(col))))

# Check for non-numeric values and report
if (anyNA(depth_data)) {
  cat("Warning: Non-numeric values detected. Rows with NA values will be removed.\n")
}

# Remove rows with NA values
depth_data <- na.omit(depth_data)

# Check if the dataset is empty after cleaning
if (nrow(depth_data) == 0) {
  stop("Error: No valid numeric data found after cleaning. Check your input file.")
}

# Calculate the average depth per site (row)
row_averages <- rowMeans(depth_data)

# Calculate the overall average depth across all sites
overall_average <- mean(row_averages)

# Print the result
cat("The overall average site depth is:", overall_average, "\n")