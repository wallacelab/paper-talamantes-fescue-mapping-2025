# Objective: This code will explore the phenotype data, specifically the alklaoid data
# Import libraries
library(tidyverse)
library(multcompView)

#Loading in data
alklaoid_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Alkaloid_Data.csv"
CT_Values_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Phenotype_Data_Delta_CT.txt"
alklaoid <- read.csv(alklaoid_loc, header = TRUE, strip.white=TRUE)
CT_Values <- read.table(CT_Values_loc, header = TRUE, strip.white=TRUE)

# Giving data parents
alklaoid[c('Maternal_Parent', 'DeleteMe')] <- str_split_fixed(alklaoid$ID, '-', 2)
alklaoid <- subset(alklaoid, select = -c(DeleteMe))

#Merging the data
colnames(CT_Values)[colnames(CT_Values) == "Treatment"] ="ID"
phenotype_Data <- merge(CT_Values, alklaoid, by.x = c("ID"), by.y = c("ID"))
colnames(phenotype_Data)[colnames(phenotype_Data) == "Plate"] ="Alkaloid_Plate"
phenotype_Data <- subset(phenotype_Data, select = -c(Maternal_Parent.y))
phenotype_Data$Maternal_Parent.x <- as.character(phenotype_Data$Maternal_Parent.x)




# Making letters for groups
# Making letters for box plot alkaloids
anovaA <- aov(ng.g ~ Maternal_Parent.x, data = phenotype_Data)
summary(anovaA)
tukeyA <- TukeyHSD(anovaA)
print(tukeyA)
TkA <- group_by(phenotype_Data, Maternal_Parent.x) %>%
  summarise(mean=mean(ng.g), quant = quantile(ng.g, probs = .75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anovaA, tukeyA) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent.x)
TkA$cld <- cld$Letters
print(TkA)

# Making letters for box plot biomass
anovaD <- aov(Delta_CT ~ Maternal_Parent.x, data = phenotype_Data)
summary(anovaD)
tukeyD <- TukeyHSD(anovaD)
print(tukeyD)
TkD <- group_by(phenotype_Data, Maternal_Parent.x) %>%
  summarise(mean=mean(Delta_CT), quant = quantile(Delta_CT, probs = .75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anovaD, tukeyD) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent.x)
TkD$cld <- cld$Letters
print(TkD)

# Boxplots Alkaloids
ggplot(phenotype_Data, aes(x=Maternal_Parent.x, y=ng.g, fill = Maternal_Parent.x, group=Maternal_Parent.x)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Alkaloid") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = TkA, aes(label = cld, x = Maternal_Parent.x, y = quant), 
            vjust = 2.5 , size = 5)

# Boxplots biomass
ggplot(phenotype_Data, aes(x=Maternal_Parent.x, y=Delta_CT, fill = Maternal_Parent.x, group=Maternal_Parent.x)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Alkaloid") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = TkD, aes(label = cld, x = Maternal_Parent.x, y = quant),
            vjust = 1 , size = 5)



