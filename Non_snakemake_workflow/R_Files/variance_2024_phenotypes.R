# Objective: This will calculate heritablitlity for 2024 data Biomass and alkaloid levels
library(tidyverse)


# Loading in data 
MLM_stats_loc = "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Tassel_Outputs/2024_only/MLM_stats_MCR50.txt"
MLM_stats <- read.table(MLM_stats_loc, header = TRUE, sep = "\t")

variance <- MLM_stats[ , c(1, 16, 17)]

v_alkaloids <- variance[375222,]
v_biomass <- variance[1,]

# Vg / (Vp + Vg)
v_alkaloids$Genetic.Var / (v_alkaloids$Residual.Var + v_alkaloids$Genetic.Var)
v_biomass$Genetic.Var / (v_biomass$Residual.Var + v_biomass$Genetic.Var)
