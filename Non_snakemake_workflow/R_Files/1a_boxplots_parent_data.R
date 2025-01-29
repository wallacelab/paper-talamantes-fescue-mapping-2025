# Purpose of this is to create boxplots of parental data from initial tests


library(tidyverse)



# Load in data
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

parental_alk <- read.csv(paste0(data_folder,"/Phenotype_Data/Preliminary_Phenotype_Data/Cleaned/Parent_Alkaloids_R.csv"), header = TRUE)
parental_CT <- read.csv(paste0(data_folder,"/Phenotype_Data/Preliminary_Phenotype_Data/Cleaned/Clone_Parent_Biomass_Data.csv"), header = TRUE)

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







