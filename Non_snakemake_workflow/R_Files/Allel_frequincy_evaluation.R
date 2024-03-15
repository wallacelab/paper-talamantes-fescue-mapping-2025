# Objective: This script evaluates how normal my allel frequencies are 
# Below are some basic info on the 320 x 315 cross
# The parents will dictate the percent chance of a success based on mendalian genetics
# The progeny are the number of succeses
# total number of progeny is 188
library(tidyverse)
library(reshape2) 
library(tidyr)


# Load in data
allel_freq_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/allel_frequencies.csv"
allel_freq <- read.csv(allel_freq_loc, header = TRUE, strip.white=TRUE)

beagle_vs_vcf_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/beagle_vs_vcf.csv"
beagle_vs_vcf <- read.csv(beagle_vs_vcf_loc, header = TRUE, strip.white=TRUE)

genotype_depth_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/genotype_depth.csv"
genotype_depth <- read.csv(genotype_depth_loc, header = FALSE, strip.white=TRUE)

prog_counts_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/Homo_Hetero_Numbers.csv"
prog_count <- read.csv(prog_counts_loc, header = TRUE, strip.white=TRUE)

# Save locations
pseudo_testcrosses_save <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/psedotest_cross_list.txt"
  
  
#Fixing the number labels
allel_freq$X <- allel_freq$X +1


# Separate the ones and zeroes on the parent side from the other data
allel_freq_1_0 <- allel_freq[allel_freq$Parents %in% c(0, 1), ]
allel_freq_not_100 <- allel_freq[allel_freq$Parents %in% c(.25, .75, .5), ]
print(nrow(allel_freq_1_0) + nrow(allel_freq_not_100))
nrow(allel_freq)

# Creating a column for the CDF
allel_freq_not_100$binom <- pbinom(allel_freq_not_100$Progeny, size=188, prob=allel_freq_not_100$Parents) 
pbinom(47, size=188, .25) 
pbinom(94, size=188, .5) 
pbinom(141, size=188, .75) 

# Making a progeny frequncy column
allel_freq$prog_freq <- allel_freq$Progeny/188

allel_freq_not_100 <- allel_freq_not_100 %>%
  mutate(binom2sided = ifelse(binom > 0.5, (1- binom)*2, binom*2))


################################################################################
### Working with progeny counts to find test crosses
################################################################################
head(prog_count)
prog_count$Het1_percent <- prog_count$Het1 / (94-prog_count$Missing)
prog_count$Het2_percent <- prog_count$Het2 / (94-prog_count$Missing)
prog_count$Homo1_percent <- prog_count$Homo1 / (94-prog_count$Missing)
prog_count$Homo2_percent <- prog_count$Homo2 / (94-prog_count$Missing)

prog_percents <- subset(prog_count, select = c ('ID', 'Het1_percent','Het2_percent','Homo1_percent','Homo2_percent'))
head(prog_percents)
cross_type <- data.frame(ID = character(length = nrow(prog_percents)))
cross_type$ID <- prog_percents$ID


# This loops tests the percents of progeny outcomes in each site and if they match  
# with an expected percent from a cross it will add that site to the predicted cross
# This is to clear the lists before the loop is run
for (x in 1:1){
tolerance <- .08
no_hits <- 0
PureHomo0 <- list()
PureHomo1 <- list()
PureHetero <- list()
Half_het_homo0 <- list()
Half_het_homo1 <- list()
HalfHomo0_HalfHomo1 <- list()
No_Hit <- list()

  for (row in 1:nrow(prog_percents)){
    Het <- prog_percents$Het1_percent[row]
    Homo0 <- prog_percents$Homo1_percent[row]
    Homo1 <- prog_percents$Homo2_percent[row]
    Site <- prog_percents$ID[row]
    print(Site)
    # Half hetero with homo 0|0 parents #1
    if (Het >= (0.5 - tolerance) & Het <= (0.5 + tolerance) & Homo0 >= (0.5 - tolerance) & Homo0 <= (0.5 + tolerance) ){
      Half_het_homo0 <- c(Half_het_homo0, paste(Site, row))
      print("HIT Half_het_homo0!")
      
    } # Half hetero with homo 1|1 parents #2
      else if (Het >= (0.5 - tolerance) & Het <= (0.5 + tolerance) & Homo1 >= (0.5 - tolerance) & Homo1 <= (0.5 + tolerance) ){
      Half_het_homo1 <- c(Half_het_homo1, paste(Site, row))
      print("HIT Half_het_homo1!")
      
    }# All Homo 0|0 parents #3
      else if (Homo0 >= (1 - tolerance) & Homo0 <= (1 + tolerance)){
      PureHomo0 <- c(PureHomo0, paste(Site, row))
      print("HIT PureHomo0!")
      
    }# All Homo 1|1 parents #4
      else if (Homo1 >= (1 - tolerance) & Homo1 <= (1 + tolerance)){
      PureHomo1 <- c(PureHomo1, paste(Site, row))
      print("HIT PureHomo1!")
      
      } # All Hetero parents #5
      else if (Het >= (0.5 - tolerance) & Het <= (0.5 + tolerance) & Homo0 >= (0.25 - tolerance) & Homo0 <= (0.25 + tolerance) & Homo1 >= (0.25 - tolerance) & Homo1 <= (0.25 + tolerance)){
      PureHetero <- c(PureHetero, paste(Site, row))
      print("HIT PureHetero!")
      
      } # Half Homo1 Half Homo2 parents
      else if (Homo0 == (0.5 - tolerance) & Homo1 == (0.5 - tolerance)){
      HalfHomo0_HalfHomo1 <- c(HalfHomo0_HalfHomo1, paste(Site, row))
      print("HIT HalfHomo0_HalfHomo1!")
      }
    
      else{
      no_hits = no_hits + 1
      cat("no hit", no_hits, "\n")
      No_Hit <- c(No_Hit, paste(Site, row))
      }
  }
}
# Grabbing the pseudo test cross and exporting it
Half_het_homo0_df <- as.data.frame(t(Half_het_homo0))
Half_het_homo0_df <- as.vector(t(Half_het_homo0_df))

Half_het_homo1_df <- as.data.frame(t(Half_het_homo1))
Half_het_homo1_df <- as.vector(t(Half_het_homo1_df))

pseudo_testcrosses <- as.data.frame(rbind(Half_het_homo0_df, Half_het_homo1_df))
pseudo_testcrosses$V1 <- as.character(pseudo_testcrosses$V1)
pseudo_testcrosses <- separate(pseudo_testcrosses, V1, into = c("V1", "Number"), sep = "\\ ")
pseudo_testcrosses <- separate(pseudo_testcrosses, V1, into = c("CHROM", "POS"), sep = "l")
pseudo_testcrosses$CHROM <- paste0(pseudo_testcrosses$CHROM, "l")
pseudo_testcrosses <- subset(pseudo_testcrosses, select = -c(Number))

write.table(pseudo_testcrosses, pseudo_testcrosses_save, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


### Secondary method (uses chi square test and p adjust)
# p.adjust

# Creating a data frame for the results 
chi_squared_results <- data.frame(ID = character(length = nrow(prog_percents)))
chi_squared_results$ID <- prog_percents$ID
chi_squared_results$Half_het_homo1 <- NA
chi_squared_results$Half_het_homo2 <- NA

# Extracting data for test 
chi_squared_data <- subset(prog_percents, select = c(Het1_percent,Homo1_percent,Homo2_percent))

# creating a function to do a chi squared test on each row.
# It looks for a specific distribution
perform_chi_squared_test0 <- function(row) {
  expected_distribution <- c(0.475, 0.475, 0.05)
  chi_squared_result <- chisq.test(row, p = expected_distribution)
  return(chi_squared_result$p.value)
}
perform_chi_squared_test1 <- function(row) {
  expected_distribution <- c(0.475, 0.05, 0.475)
  chi_squared_result <- chisq.test(row, p = expected_distribution)
  return(chi_squared_result$p.value)
}

# # This loop ensures that I get rid of any 1's in the data that would make the chi sq not work.
# # This did not work.
# for (row in range(1:nrow(chi_squared_data))){
#   het <- chi_squared_data$Het1_percent[row]
#   homo1 <- chi_squared_data$Homo1_percent[row]
#   homo2 <- chi_squared_data$Homo2_percent[row]
#   if (het == 1) {
#     chi_squared_data$Het1_percent[row] <- chi_squared_data$Het1_percent[row] -.0002
#     chi_squared_data$Homo1_percent[row] <-chi_squared_data$Homo1_percent[row] + .0001
#     chi_squared_data$Homo2_percent[row] <- chi_squared_data$Homo2_percent[row] + .0001
#   }
#   else if (homo1 == 1) {
#     chi_squared_data$Het1_percent[row] <- chi_squared_data$Het1_percent[row] + .0001
#     chi_squared_data$Homo1_percent[row] <-chi_squared_data$Homo1_percent[row] - .0002
#     chi_squared_data$Homo2_percent[row] <- chi_squared_data$Homo2_percent[row] + .0001
#   }
#   else if (homo2 == 1) {
#     chi_squared_data$Het1_percent[row] <- chi_squared_data$Het1_percent[row] + .0001
#     chi_squared_data$Homo1_percent[row] <-chi_squared_data$Homo1_percent[row] + .0001
#     chi_squared_data$Homo2_percent[row] <- chi_squared_data$Homo2_percent[row] - .0002
#   }
# }


p_values_test0 <- apply(chi_squared_data, 1, perform_chi_squared_test0)
p_values_test1 <- apply(chi_squared_data, 1, perform_chi_squared_test1)

chi_squared_results$Half_het_homo1 <- round(p_values_test0, digits = 4)
chi_squared_results$Half_het_homo2 <- round(p_values_test1, digits = 4)

expected_distribution <- c(0.475, 0.475, 0.05)
simudata <- data.frame(.50,.50,.0)
simudata <- data.frame(.50,.0,.50)
simudata <- data.frame(0,1,.0)
simudata <- data.frame(0,0,1)
simudata <- data.frame(.25,.50,.25)
simudata <- data.frame(.30,.10,.60)

chi_squared_result <- chisq.test(simudata, p = expected_distribution)
chi_squared_result$p.value
################################################################################
###  Making graphs
################################################################################
ggplot(allel_freq_not_100, aes(x=binom)) +
  geom_histogram(bins = 50) +
  labs(x = "Cumulative Distibution", y = "Count")

# This plot shows the results of a 2 sided bionomial distibution test
ggplot(allel_freq_not_100, aes(x=binom2sided)) +
  geom_histogram(bins = 50) +
  labs(x = "Cumulative Distibution", y = "Count")

# Plotting the frequencies on the same plot
Frequncies_only <- subset(allel_freq, select = c("Site", "prog_freq", "Parents" ))
Frequncies_only <- gather(Frequncies_only, key = "origin", value = "value", prog_freq, Parents)
ggplot(Frequncies_only, aes(x = value, fill = origin, color = origin)) +
  geom_histogram(alpha=.5, position = "identity", bins = 50)
allel_freq_small <- head(allel_freq, 100)

ggplot(allel_freq_small, aes(x = Site)) +
  geom_bar(aes(y = prog_freq), stat = "identity", position = "dodge", fill = "blue", alpha = 0.7) +
  geom_bar(aes(y = Parents), stat = "identity", position = "dodge", fill = "red", alpha = 0.7)

#Plotting beagle vs vcf files
beagle_vs_vcf_tots <- colSums(beagle_vs_vcf[, c("Same", "Filled", "Homo2Het", "Het2Homo","Het2Het")])
beagle_vs_vcf_tots <- melt(beagle_vs_vcf_tots, id = c("Same","Filled", "Homo2Het", "Het2Homo", "Het2Het"))
beagle_vs_vcf_tots$Category <- rownames(beagle_vs_vcf_tots)
beagle_vs_vcf_tots$value <- as.numeric(beagle_vs_vcf_tots$value)

ggplot(data=beagle_vs_vcf_tots, aes(x=Category, y=value, fill=Category)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Comparison Result") + ylab("Count") +
  ggtitle("VCF to Beagle Allele Changes")


# Plotting the genotype depth 
genotype_depth <- data.frame(rowname = rownames(genotype_depth), genotype_depth)

ggplot(data = genotype_depth, aes(x = rowname, y = V2)) +
  geom_point(shape = 14) +
  ggtitle("Distribution of SNP Location Depths") +
  xlab("SNP") + ylab("Count")

# Plotting predicted cross data
all_cross_types <- c("PureHomo0", "PureHomo1", "PureHetero", "Half_het_homo0", "Half_het_homo1", "HalfHomo0_HalfHomo1", "No_Hit")
cross_counts <- tibble(
  cross_type = c(1:7),
  Count = c(length(PureHomo0),length(PureHomo1),length(PureHetero),length(Half_het_homo0),length(Half_het_homo1),length(HalfHomo0_HalfHomo1),length(No_Hit))
)
cross_counts$cross_type <- all_cross_types

ggplot(data = cross_counts, aes(x = cross_type, y = Count, fill = cross_type)) +
  geom_bar(stat = "identity", color = "black") +
  ggtitle("Classification of Sites with a tolerance of ", tolerance) +
  xlab("Predicted Cross Type") + ylab("Count")


