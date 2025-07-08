# Purpose: This script will analyze the entire set of 2023 data
# This works with the data that is already residuals.

library(tidyverse)
library(gridExtra)
library(lme4)
library(qqman)
library(multcompView)   
library(gridExtra)

################################################################################
# Load Data
################################################################################
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"


phenotype_Data_2023 <- read.table(paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_only.txt") , header = TRUE,sep = '\t')
Meta_Data_2023 <- read.table(paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_data_and_metadata.txt") , header = TRUE,sep = '\t')
PCA_FlexSeq <- read.table(paste0(data_folder,"/Tassel_Outputs/FlexSeq_Genos/PCA_flexseq_genotypes.txt") , header = TRUE,sep = '\t')



# Merge the Data
phenotype_Data_2023 <- merge(phenotype_Data_2023, Meta_Data_2023, by = "ID")
# Remove duplicate columns that have ".y" in the name, keeping ".x"
phenotype_Data_2023 <- phenotype_Data_2023[, !grepl("\\.y$", colnames(phenotype_Data_2023))]
colnames(phenotype_Data_2023) <- sub("\\.x$", "", colnames(phenotype_Data_2023))



################################################################################
# Scatter Plot
################################################################################head(phenotype_Data_2023)
nrow(phenotype_Data_2023)

# linear model
model <- lm(Alkaloids_Res ~ Delta_CT_adj_Res, data = phenotype_Data_2023)
rsq <- summary(model)$r.squared #This is low so I think I should use spearman.


# Function to calculate r squared and make scatter plot
scatterplot_phenos <- function(dataset, Alkaloidcol, DeltaCTcol, Title) {
  # Perform Spearman correlation
  spearmodel <- cor.test(dataset[[Alkaloidcol]], 
                         dataset[[DeltaCTcol]], 
                         method = "spearman", 
                         use = "complete.obs")  # Ignore NAs
  rho <- spearmodel$estimate
  rsq <- rho^2
  p <- spearmodel$p.value
  
  # Create the scatter plot
  plot1 <- ggplot(dataset, aes(x = .data[[DeltaCTcol]], y = .data[[Alkaloidcol]])) +
    geom_point() +  # Scatter plot points
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
    annotate("text", 
             x = min(dataset[[DeltaCTcol]], na.rm = TRUE), 
             y = max(dataset[[Alkaloidcol]], na.rm = TRUE), 
             label = paste("Spearman R-squared = ", round(rsq, 3)), 
             hjust = 0, vjust = -1, size = 5, color = "red") +
    annotate("text", 
             x = min(dataset[[DeltaCTcol]], na.rm = TRUE), 
             y = max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
               0.05 * (max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
                         min(dataset[[Alkaloidcol]], na.rm = TRUE)), 
             label = paste("P-value = ", format(p, digits = 3, scientific = TRUE)), 
             hjust = 0, vjust = -1, size = 5, color = "red") +
    labs(title = Title, 
         x = "Efficiency Adjusted CT Ratio", 
         y = "Residual Alkaloids") + 
    theme_bw() +
    theme(text = element_text(size = 20))
  
  return(plot1)
}


scatterplot_phenos(phenotype_Data_2023,"Alkaloids_Res","Delta_CT_adj_Res", "Alkaloid vs Relative Biomass of The Large Population")


################################################################################
# Making Boxplots
################################################################################

# Adding Cross and Paternal Data
Fathers <- read.table(paste0(data_folder,"/Phenotype_Data/Meta_Data/Mother_Father_Data.csv") , header = TRUE,sep = ',')
Fathers <- Fathers %>%
  rename(Maternal_Parent = Mother, Paternal_Parent = Father)

phenotype_Data_2023 <- merge(phenotype_Data_2023, Fathers, by = c("ID", "Maternal_Parent"))
phenotype_Data_2023$Cross <- paste(phenotype_Data_2023$Maternal_Parent, phenotype_Data_2023$Paternal_Parent, sep = "x")


# Makes a graph from the Crosses from ALkaloid Data
anova <- aov(Alkaloids_Res ~ Cross, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Cross) %>%
  summarise(mean=mean(Alkaloids_Res), quant = quantile(Alkaloids_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Cross)
Tk$cld <- cld$Letters
print(Tk)
ggplot(phenotype_Data_2023, aes(x=Cross, y=Alkaloids_Res, fill = Cross, group=Cross)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("ALkaloid Residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Cross, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

# Makes a graph from the Maternal Parents and alkaloids
phenotype_Data_2023$Maternal_Parent <- as.factor(phenotype_Data_2023$Maternal_Parent)
anova <- aov(Alkaloids_Res ~ Maternal_Parent, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Maternal_Parent) %>%
  summarise(mean=mean(Alkaloids_Res), quant = quantile(Alkaloids_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent)
Tk$cld <- cld$Letters
print(Tk)
ggplot(phenotype_Data_2023, aes(x=Maternal_Parent, y=Alkaloids_Res, fill = Maternal_Parent, group=Maternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("ALkaloid Residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Maternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

# Makes a graph from the Maternal Parents and ADJ delta CT Ratios
phenotype_Data_2023$Maternal_Parent <- as.factor(phenotype_Data_2023$Maternal_Parent)
anova <- aov(Delta_CT_adj_Res ~ Maternal_Parent, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Maternal_Parent) %>%
  summarise(mean=mean(Delta_CT_adj_Res), quant = quantile(Delta_CT_adj_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent)
Tk$cld <- cld$Letters
print(Tk)
ggplot(phenotype_Data_2023, aes(x=Maternal_Parent, y=Delta_CT_adj_Res, fill = Maternal_Parent, group=Maternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Relative Biomass") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Maternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)




################################################################################
# Making PCA plot 
################################################################################

#Subsetting maternal parent
Mothers <- subset(Meta_Data_2023, select = c(ID, Maternal_Parent))
PCA_FlexSeq <- PCA_FlexSeq %>%
  rename(ID = Taxa)
PCA_FlexSeq <- merge(PCA_FlexSeq, Mothers, by = "ID")
PCA_FlexSeq$Maternal_Parent <- as.factor(PCA_FlexSeq$Maternal_Parent)

# The plot
ggplot(PCA_FlexSeq, aes(x = PC1, y = PC2, color = Maternal_Parent)) +
  geom_point(alpha = 0.5) +
  theme_minimal()


################################################################################
# Hapmat plot
################################################################################

########### Example 
MLM_DRT_Filters_residuals <- read.table(paste0(data_folder,"/Tassel_Outputs/FlexSeq_Genos/MLM_2023_phenotypes_FlexSeq_Genos.txt") , header = TRUE,sep = '\t')
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[MLM_DRT_Filters_residuals$Trait == "Alkaloids_Res", ]
MLM_DRT_Filters_residuals$Chr <- as.numeric(factor(MLM_DRT_Filters_residuals$Chr))
MLM_DRT_Filters_residuals <- MLM_DRT_Filters_residuals[-c(1), ]
# BOnferonii line
alpha <- 0.05
num_tests <- nrow(MLM_DRT_Filters_residuals)  # Number of SNPs or tests
bonferroni_threshold <- alpha / num_tests
# Making FDR threshold
MLM_DRT_Filters_residuals$padj <- p.adjust(MLM_DRT_Filters_residuals$p, method = "BH")
FDR_threshold <- max(MLM_DRT_Filters_residuals$p[MLM_DRT_Filters_residuals$padj <= alpha], na.rm = TRUE)
if (!is.numeric(FDR_threshold) || is.na(FDR_threshold) || FDR_threshold == -Inf) {
  FDR_threshold <- 1e-5
}
print(FDR_threshold)
Alkaloids <- manhattan(MLM_DRT_Filters_residuals, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
                       ylim = c(0, 6), main = "MLM on Alkaloid Residuals in FlexSeq Genotypes",
                       genomewideline = -log10(bonferroni_threshold),
                       suggestiveline = -log10(FDR_threshold),
                       col = c("blue", "red", "darkgrey", "purple"))


################################################################################
# Box plots but for paternal parents
################################################################################


# Makes a graph from the Crosses from ALkaloid Data
anova <- aov(Alkaloids_Res ~ Cross, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Cross) %>%
  summarise(mean=mean(Alkaloids_Res), quant = quantile(Alkaloids_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Cross)
Tk$cld <- cld$Letters
print(Tk)
ggplot(phenotype_Data_2023, aes(x=Cross, y=Alkaloids_Res, fill = Cross, group=Cross)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Paternal Lines") +
  ylab("ALkaloid Residuals") +
  scale_fill_discrete(name = "Paternal Lines") +
  geom_text(data = Tk, aes(label = cld, x = Cross, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

# Makes a graph from the Maternal Parents and alkaloids
phenotype_Data_2023$Paternal_Parent <- as.factor(phenotype_Data_2023$Paternal_Parent)
anova <- aov(Alkaloids_Res ~ Paternal_Parent, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Paternal_Parent) %>%
  summarise(mean=mean(Alkaloids_Res), quant = quantile(Alkaloids_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Paternal_Parent)
Tk$cld <- cld$Letters
print(Tk)
alk <- ggplot(phenotype_Data_2023, aes(x=Paternal_Parent, y=Alkaloids_Res, fill = Paternal_Parent, group=Paternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Paternal Lines") +
  ylab("Alkaloid Residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  ggtitle("A") + 
  theme(legend.position = "none") + 
  geom_text(data = Tk, aes(label = cld, x = Paternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

# Makes a graph from the Maternal Parents and ADJ delta CT Ratios
phenotype_Data_2023$Paternal_Parent <- as.factor(phenotype_Data_2023$Paternal_Parent)
anova <- aov(Delta_CT_adj_Res ~ Paternal_Parent, data = phenotype_Data_2023)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(phenotype_Data_2023, Paternal_Parent) %>%
  summarise(mean=mean(Delta_CT_adj_Res), quant = quantile(Delta_CT_adj_Res, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Paternal_Parent)
Tk$cld <- cld$Letters
print(Tk)
bio <- ggplot(phenotype_Data_2023, aes(x=Paternal_Parent, y=Delta_CT_adj_Res, fill = Paternal_Parent, group=Paternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Paternal Lines") +
  ylab("Relative Biomass") +
  scale_fill_discrete(name = "Paternal Lines") +
  ggtitle("B") + 
  geom_text(data = Tk, aes(label = cld, x = Paternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)


grid.arrange(alk, bio, ncol = 2, widths = c(1.6, 1.9))
