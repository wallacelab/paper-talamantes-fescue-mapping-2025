# Objective: This script evaluates how normal my allel frequencies are 
# Below are some basic info on the 320 x 315 cross
# The parents will dictate the percent chance of a success based on mendalian genetics
# The progeny are the number of succeses
# total number of progeny is 188
library(tidyverse)
library(reshape2) 

# Load in data
allel_freq_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/allel_frequencies.csv"
allel_freq <- read.csv(allel_freq_loc, header = TRUE, strip.white=TRUE)

#Fixing the number labels
allel_freq$X <- allel_freq$X +1

# Creating a column for the CDF
allel_freq$Parents2 <- allel_freq$Parents -.001
allel_freq$binom <- dbinom(allel_freq$Progeny, size=188, prob=allel_freq$Parents) 



