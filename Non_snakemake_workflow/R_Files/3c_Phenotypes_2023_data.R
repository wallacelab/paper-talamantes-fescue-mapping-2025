# Purpose: This script will analyze the entire set of 2023 data
# This works with the data that is already residuals.

library(tidyverse)
library(gridExtra)
library(lme4)
################################################################################
# Load Data
################################################################################
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"


phenotype_Data_2023 <- read.table(paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_only.txt") , header = TRUE,sep = '\t')
Meta_Data_2023 <- read.table(paste0(data_folder,"/Phenotype_Data/2023_Data/3a_all_2023_residual_data_and_metadata.txt") , header = TRUE,sep = '\t')

# Merge the Data
phenotype_Data_2023 <- merge(phenotype_Data_2023, Meta_Data_2023, by = "ID")
# Remove duplicate columns that have ".y" in the name, keeping ".x"
phenotype_Data_2023 <- phenotype_Data_2023[, !grepl("\\.y$", colnames(phenotype_Data_2023))]
colnames(phenotype_Data_2023) <- sub("\\.x$", "", colnames(phenotype_Data_2023))



############### Scatter plots #######################
head(phenotype_Data_2023)
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
