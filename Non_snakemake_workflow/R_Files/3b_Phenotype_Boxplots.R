# Purpose: This will make boxplots of the phenotypes 


library(tidyverse)
library(multcompView)
library(gvlma)
library(ggpubr)
library(car)
library(reshape2)
library(grid)
library(gridExtra)
library(ggpubr)

################################################################################
# Load Data
################################################################################
data_folder = "/home/darrian/Documents/Mapping_and_QTL/Data"

pg_Residual_data_avg <- read.table(paste0(data_folder,"/Phenotype_Data/pg_Residual_data_avg"), header = TRUE)
parent_list <- read.table(paste0(data_folder, "/Lists/Parent_Progeny_Lists/314_Star_Cross_Parents.txt"), header = FALSE)
parents <- read.csv(paste0(data_folder, "/Lists/usable_predicted_parents_double.csv"), header = TRUE)


################################################################################
# Fixing the data sets

colnames(parent_list) <- c("ID", "Cross")

head(parent_list)

pg_Residual_data_avg <- merge(pg_Residual_data_avg, parent_list, by = "ID")
head(pg_Residual_data_avg)

# Setting Mother and Father
data_with_parents <- data.frame()
data_with_parents <- parents %>%
  mutate(
    Mother = as.numeric(sub("-.*", "", Progeny)), # Extract maternal parent from Progeny
    Father = ifelse(Mother == Parent1, Parent2, Parent1) # Assign the paternal parent
  ) %>%
  select(ID = Progeny, Mother, Father)
data_with_parents <- data_with_parents %>%
  mutate(ID = gsub("-", "_", ID))
# View the result
head(data_with_parents)
head(pg_Residual_data_avg)
pg_Residual_data_avg <- merge(pg_Residual_data_avg,data_with_parents, by = "ID")
pg_Residual_data_avg$Mother <- as.character(pg_Residual_data_avg$Mother)
pg_Residual_data_avg$Father <- as.character(pg_Residual_data_avg$Father)



################################################################################
# Making Boxplots
################################################################################

pg_Residual_data_avg <- pg_Residual_data_avg %>%
  filter(!ID %in% c(310, 312, 314))

anova <- aov(Alkaloids_Res_avg ~ Cross, data = pg_Residual_data_avg)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(pg_Residual_data_avg, Cross) %>%
  summarise(mean=mean(Alkaloids_Res_avg), quant = quantile(Alkaloids_Res_avg, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Cross)
Tk$cld <- cld$Letters
print(Tk)
ggplot(pg_Residual_data_avg, aes(x=Cross, y=Alkaloids_Res_avg, fill = Cross, group=Cross)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("ALkaloid Residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Cross, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)





anova <- aov(DeltaCT_adj_Res_avg ~ Cross, data = pg_Residual_data_avg)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(pg_Residual_data_avg, Cross) %>%
  summarise(mean=mean(DeltaCT_adj_Res_avg), quant = quantile(DeltaCT_adj_Res_avg, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Cross)
Tk$cld <- cld$Letters
print(Tk)
ggplot(pg_Residual_data_avg, aes(x=Cross, y=DeltaCT_adj_Res_avg, fill = Cross, group=Cross)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Delta CT residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Cross, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

################################################################################
# Looking at mothers of individuals

anova <- aov(Alkaloids_Res_avg ~ Mother, data = pg_Residual_data_avg)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(pg_Residual_data_avg, Mother) %>%
  summarise(mean=mean(Alkaloids_Res_avg), quant = quantile(Alkaloids_Res_avg, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Mother)
Tk$cld <- cld$Letters
print(Tk)
alkaloid_plot <- ggplot(pg_Residual_data_avg, aes(x=Mother, y=Alkaloids_Res_avg, fill = Mother, group=Mother)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Mother Plant") +
  ylab("Alkaloid Residuals") +
  scale_fill_discrete(name = "Mother Plant") +
  theme(text = element_text(size = 20)) +
  geom_text(data = Tk, aes(label = cld, x = Mother, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)


anova <- aov(DeltaCT_adj_Res_avg ~ Mother, data = pg_Residual_data_avg)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(pg_Residual_data_avg, Mother) %>%
  summarise(mean=mean(DeltaCT_adj_Res_avg), quant = quantile(DeltaCT_adj_Res_avg, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Mother)
Tk$cld <- cld$Letters
print(Tk)
CT_plot <- ggplot(pg_Residual_data_avg, aes(x=Mother, y=DeltaCT_adj_Res_avg, fill = Mother, group=Mother)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Mother Plant") +
  ylab("CT Ratio Residuals") +
  scale_fill_discrete(name = "Mother Plant") +
  theme(text = element_text(size = 20)) +
  geom_text(data = Tk, aes(label = cld, x = Mother, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)


combined_plot <- grid.arrange(
  alkaloid_plot + theme(legend.position = "none"), 
  CT_plot + theme(legend.position = "none"), 
  ncol = 2,
  top = textGrob("Differences in Maternal Parents", 
                 gp = gpar(fontsize = 25, fontface = "bold"))
) # Combine plots


# Show the plot
combined_plot_with_title


combined_plot
head(pg_Residual_data_avg)


