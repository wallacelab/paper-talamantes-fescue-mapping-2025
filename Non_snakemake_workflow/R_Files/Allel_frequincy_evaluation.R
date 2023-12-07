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

beagle_vs_vcf_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/beagle_vs_vcf.csv"
beagle_vs_vcf <- read.csv(beagle_vs_vcf_loc, header = TRUE, strip.white=TRUE)

genotype_depth_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data/genotype_depth.csv"
genotype_depth <- read.csv(genotype_depth_loc, header = FALSE, strip.white=TRUE)
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


