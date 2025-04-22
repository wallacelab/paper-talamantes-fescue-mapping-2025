# Purpose of this is to create boxplots of parental data from initial tests


library(tidyverse)



# Load in data
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

parental_alk <- read.csv(paste0(data_folder,"/Phenotype_Data/Preliminary_Phenotype_Data/Cleaned/Parent_Alkaloids_R.csv"), header = TRUE)
parental_CT <- read.csv(paste0(data_folder,"/Phenotype_Data/Preliminary_Phenotype_Data/Cleaned/Clone_Parent_Biomass_Data.csv"), header = TRUE)
parental_phenos <- read.csv(paste0(data_folder,"/Phenotype_Data/Preliminary_Phenotype_Data/Cleaned/Parental_Phenotypes_Clean.csv"), header = TRUE)

# Clean up data
parental_alk$Parent <- as.character(parental_alk$Parent)
parental_CT$ID <- as.character(parental_CT$ID)
parental_CT <- parental_CT %>%
  rename(Parent = ID)

################################################################################
#Box plot data 
################################################################################


# Alkaloid Box Plots
anova <- aov(ng.g ~ Parent, data = parental_alk)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(parental_alk, Parent) %>%
  summarise(mean=mean(ng.g), quant = quantile(ng.g, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Parent)
Tk$cld <- cld$Letters
print(Tk)
alk <- ggplot(parental_alk, aes(x=Parent, y=ng.g, fill = Parent, group=Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  theme(
    text = element_text(size = 20)) +
  xlab("Parent") +
  ylab("Ergot Alkaloid Amount") +
  scale_fill_discrete(name = "Parent") +
  geom_text(data = Tk, aes(label = cld, x = Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)



# CT ratio boxplots
anova <- aov(CP_Ratio ~ Parent, data = parental_CT)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(parental_CT, Parent) %>%
  summarise(mean=mean(CP_Ratio), quant = quantile(CP_Ratio, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Parent)
Tk$cld <- cld$Letters
print(Tk)
CT <- ggplot(parental_CT, aes(x=Parent, y=CP_Ratio, fill = Parent, group=Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  theme(
    text = element_text(size = 20)) +
  xlab("Parent") +
  ylab("Delta CT") +
  ylim(.7,1.4) +
  scale_fill_discrete(name = "Parent") +
  geom_text(data = Tk, aes(label = cld, x = Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)



#Combine the plots
combined_plot <- grid.arrange(
  alk + theme(legend.position = "none"), 
  CT + theme(legend.position = "none"), 
  ncol = 2,
  top = textGrob("Differences in Maternal Parent Clones", 
                 gp = gpar(fontsize = 25, fontface = "bold")))



################################################################################
# Scatter Plots and Correlation
################################################################################
# Removing Outliars
parental_phenos <- parental_phenos[parental_phenos$Treatment != "314C1R1", ]

model <- lm(Ng.g ~ CP_Ratio, data = parental_phenos)
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
             hjust = 0, vjust = 0, size = 5, color = "red") +
    annotate("text", 
             x = min(dataset[[DeltaCTcol]], na.rm = TRUE), 
             y = max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
               0.05 * (max(dataset[[Alkaloidcol]], na.rm = TRUE) - 
                         min(dataset[[Alkaloidcol]], na.rm = TRUE)), 
             label = paste("P-value = ", format(p, digits = 3, scientific = TRUE)), 
             hjust = 0, vjust = .5, size = 5, color = "red") +
    labs(title = Title, 
         x = "Efficiency Adjusted CT Ratio", 
         y = "Residual Alkaloids") + 
    theme_bw() +
    theme(text = element_text(size = 20))
  
  return(plot1)
}


scatterplot_phenos(parental_phenos,"Ng.g","CP_Ratio", "Alkaloid vs Relative Biomass \n of Progenitors")




