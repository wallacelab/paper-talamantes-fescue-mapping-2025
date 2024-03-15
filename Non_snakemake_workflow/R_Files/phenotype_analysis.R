# Objective: This code will explore the phenotype data, specifically the alklaoid data
# Import libraries
library(tidyverse)
library(multcompView)
library(gvlma)
library(mbQTL)
library(ggpubr)
library(car)

###############################################################################
# CT Data filtering and adding
###############################################################################
# Adding in the data sets

# This first set of data is all the data corresponding to the first set of standards
all_2x2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_2x2_Red_Blue.csv"
all_g3p4_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/all_g3p4_correspond_to_Red_Blue_Standards.csv"
# This set of data corresponds to the second set of standards and a few redoes. Will be added towrds end
all_2x2_2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/1041-end-redoes_2x2.csv"
all_g3p4_2_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/1041-end-redoes_g3p4.csv"  
# This set of data is for parents and cross 305x320, will be added twords end
all_315x320_2x2_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Raw_Data/CT_Values/315x320_2x2.csv"
all_315x320_g3p4_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Raw_Data/CT_Values/315x320_G3P4.csv"

all_2x2_1 <- read.csv(all_2x2_loc, header = TRUE, strip.white=TRUE)
all_g3p4_1 <- read.csv(all_g3p4_loc, header = TRUE, strip.white=TRUE)
all_2x2_2 <- read.csv(all_2x2_2_loc, header = TRUE, strip.white=TRUE)
all_g3p4_2 <- read.csv(all_g3p4_2_loc, header = TRUE, strip.white=TRUE)
all_315x320_2x2 <- read.csv(all_315x320_2x2_loc, header = TRUE, strip.white=TRUE)
all_315x320_g3p4 <- read.csv(all_315x320_g3p4_loc, header = TRUE, strip.white=TRUE)

all_2x2_1$Standard <- "First_combo"
all_2x2_2$Standard <- "Second_combo"
all_g3p4_1$Standard <- "First_combo"
all_g3p4_2$Standard <- "Second_combo"
all_315x320_2x2$Standard <- "Pre_combo"
all_315x320_g3p4$Standard <- "Pre_combo"

all_2x2_samples <- rbind(all_2x2_1, all_2x2_2, all_315x320_2x2)
all_g3p4_samples <- rbind(all_g3p4_1, all_g3p4_2, all_315x320_g3p4)

########### 
# Functions
###########
# This function takes the data set I combined using cat, all_2x2 is an example
# It will filter out NA's, 0's, standards
# It spits out a table of means for all samples that passed filtering. 
Data_Filtering <- function(Data){
  Data$Cp <- as.numeric(Data$Cp)
  Data[Data == ''] <- NA
  # filtering data 2x2
  for (x in nrow(Data):1){
    if (is.na(Data[x,3])){
      Data <- Data[-x,]
    }
    else if (is.na(Data[x,2])){
      Data <- Data[-x,]
    }
    else if (Data[x,2]=="0"){
      Data <- Data[-x,]
    }
    else if (Data[x,1]=="Pos"){
      Data <- Data[-x,]
    }
    else if (str_detect(Data[x,3], "std") | Data[x,3] == 0 | Data[x,3] == "water"){
      Data <- Data[-x,]
    }
  }
  Data_Mean <-
    Data %>%
    group_by(Data_Set,Treatment) %>%
    summarise(MeanCP = mean(Cp))
  return(Data_Mean)
}

# This pair of functions allows you to put in a CP value and it will output an adjusted CP value that takes efficiency into consideration
# Step 1 Calculate line for a single set of standards
# Step 2 Calculate efficiency for standards
# Step 2.5 Calculate the outliars and remove them
# Step 3 Use efficiency to find adjusted CP value
CpAdjuster <- 
  function(model, Data){
    Data$adjCP <- 0
    slope <- model$coefficients[2]
    E = -1+10^(-1/slope)
    
    for (i in 1:nrow(Data)){
      ogCP <- Data[i,3]
      ratio =  (2^ogCP) / (2*E)^ogCP # Finding ratio of expected vs actual
      adjuster = log(ratio, base = 2) # Finding amount needed to adjust CP videa
      adjCP = round(ogCP - adjuster, 3) #finalizes CP adjustment
      Data$adjCP[i] <- adjCP
    }
    Data$Efficiency <- E
    return(Data)
  }



# This next function works with 2 data sets. Mean standards and mean data. 
# It will use the standards to calculate efficiency and use that to spit out the data with adjusted CP values
# This function uses another function called "CpAdjuster"
# THis function also removes outliars who are out of standard range or 3x std away from mean
CpAdjusterP2 <- function(Mean_Standards, Mean_Data, Data_Sets){
  Data_Sets <- as.matrix(Data_Sets)
  All_Data = list()
  for (i in 1:nrow(Data_Sets)){
    
    Set <- Data_Sets[i]
    Current_Standard <- subset(Mean_Standards, Data_Set == Set)
    Current_Data <- subset(Mean_Data, Data_Set == Set)
    
    #Removing outliars
    meanCP <- mean(Current_Data$MeanCP)
    SdCP <- sd(Current_Data$MeanCP)
    cutoff1 <- meanCP + 3*SdCP
    cutoff2 <- meanCP - 3*SdCP
    cutoffmin <- min(Current_Standard$MeanCP)
    cutoffmax <- max(Current_Standard$MeanCP)
    Before <- nrow(Current_Data)
    Current_Data <- subset(Current_Data, MeanCP < cutoff1 & MeanCP > cutoff2 & MeanCP > cutoffmin & MeanCP < cutoffmax)
    after <- nrow(Current_Data)
    
    print(paste0(Set, " Dropped ", Before-after))
    
    modelSet <- lm(Current_Standard$MeanCP ~ Current_Standard$LogCopyNumber)
    y_intercept <- modelSet$coefficients[1]
    slope <- modelSet$coefficients[2]
    summary(modelSet)
    Current_Data$LogCopyNumber <- round(((Current_Data$MeanCP)-y_intercept)/slope, digits = 3)
    Current_Data$CopyNumber <- 10^Current_Data$LogCopyNumber # Reverse the log
    
    All_DataC <- CpAdjuster(modelSet, Current_Data) # Need this to remove outliars too
    All_Data[[i]] <-  All_DataC
  }
  All_Data <- bind_rows(All_Data)
  return(All_Data)
}

# seperating the standards and getting the means of them, calculating ng of DNA
# function takes full data set and seperates standards to find their means
# The function also calculates the copy number and log copy number
FindStandardMeans <- function(stdData, length){
  stdData$Cp <- as.numeric(stdData$Cp)
  stdData$Concentration <- as.numeric(stdData$Concentration)
  for (x in nrow(stdData):1){
    if (stdData[x,2]=="" | stdData[x,2]=="0" | nchar(stdData[x,3])!=4){
      stdData <- stdData[-x,]}
  }
  Standard_Means <-
    stdData %>%
    group_by(Data_Set,Treatment) %>%
    summarise(MeanCP = mean(Cp), NgDNA = mean(Concentration)*5) # we add 5 ul of sample
  Standard_Means$CopyNumber <- (Standard_Means$NgDNA * 6.02214076*10^23/ (length * 650 * 10^9))
  # NgDNA*1 mole / length * 1 mole of base pairs in ng
  Standard_Means$LogCopyNumber <- log10(Standard_Means$CopyNumber)
  return(Standard_Means)
}

## Removing Redone samples
redoes_loc <- "/home/drt06/Documents/QPCR_Data_Wrangler/Program/int_files/Data_for_Project/Sample_Redo_List.csv"
redoes <- read.csv(redoes_loc, header = FALSE)
redoner<- function(redoes,sampledata){
  for (x in nrow(redoes):1){
    for (y in nrow(sampledata):1){
      if(redoes[x,] == sampledata[y,3]){
        if(sampledata[y,7] != "1281-1325-R"){
          sampledata <- sampledata[-y,]
        }
      }  
    }
  }
  return(sampledata)
}

############
# End of functions
############
all_2x2_samples <- redoner(redoes,all_2x2_samples)
all_g3p4_samples <- redoner(redoes,all_g3p4_samples)

# Removing specific redone samples
all_2x2_samples <- all_2x2_samples[!(all_2x2_samples$Treatment == "303-1-41" & all_2x2_samples$Data_Set == "81-160" ),]

all_2x2_samples$Concentration <- as.numeric(all_2x2_samples$Concentration)
all_g3p4_samples$Concentration <- as.numeric(all_g3p4_samples$Concentration)

# Adding in Standard combination data
standard_group <- subset(all_2x2_samples, select = c("Data_Set","Standard"))
standard_group <- standard_group[!duplicated(standard_group), ]
standard_group <- standard_group[!(is.na(standard_group$Data_Set) | standard_group$Data_Set==""), ]

# Filtering and getting means of all CP values
# takes out obviously bad samples, 0's, missing, ect.
Data_g3p4_Means_all <- Data_Filtering(all_g3p4_samples) # method filters data and leaves only samples
Data_2x2_Means_all <- Data_Filtering(all_2x2_samples)

# Method seperates standards from data,  to get their mean and adds copy number and log copy number
std_rest_g3p4_means <- FindStandardMeans(all_g3p4_samples,230)
# Needed to double filter 2x2 for some reason
all_2x2_samples2 <- all_2x2_samples
for (x in nrow(all_2x2_samples2):1){
  if (all_2x2_samples2[x,2]=="" | all_2x2_samples2[x,2]=="0" | nchar(all_2x2_samples2[x,3])!=4){
    all_2x2_samples2 <- all_2x2_samples2[-x,]
  }
}
std_rest_2x2_means <- FindStandardMeans(all_2x2_samples2,118)

# method adjusts the CP value based on efficiency
# Also filters out samples that are 3x away from std, too low, too high, ect. 
Data_Sets <- unique(std_rest_g3p4_means$Data_Set) #Creating list of data sets for next method

Data_g3p4_AdjCP_all <- CpAdjusterP2(std_rest_g3p4_means, Data_g3p4_Means_all,Data_Sets) 
Data_2x2_AdjCP_all <- CpAdjusterP2(std_rest_2x2_means, Data_2x2_Means_all,Data_Sets) 


################ Combining data and getting the delta CT values  ###############

#Combining data set and naming columns
all_Data <- merge(Data_2x2_AdjCP_all,Data_g3p4_AdjCP_all, by.x = c("Treatment", "Data_Set"), by.y = c("Treatment", "Data_Set"))
Columns <- colnames(all_Data)
Columns <- gsub("\\.y$", ".Fescue", Columns)
Columns <- gsub("\\.x$", ".Epichloe", Columns)
colnames(all_Data) <- Columns 
all_Data$adjCP.Fescue <- as.numeric(all_Data$adjCP.Fescue)
all_Data$adjCP.Epichloe <- as.numeric(all_Data$adjCP.Epichloe)

# Dela CT is Fescue - EPichloe (adjusted)
all_Data$Delta_CT <- all_Data$adjCP.Fescue - all_Data$adjCP.Epichloe

# Delta Ratio is Fescue / Epichloe (adjusted)
all_Data$Fes_to_Epi_Ratio <- all_Data$adjCP.Fescue / all_Data$adjCP.Epichloe

# Introducing Standard Combination Data
all_Data <- merge(all_Data,standard_group, by.x = c("Data_Set"), by.y =c("Data_Set"))
# Introducing Maternal Parent Data
all_Data[c('Maternal_Parent', 'DeleteMe')] <- str_split_fixed(all_Data$Treatment, '-', 2)
all_Data <- all_Data[ , !names(all_Data) %in% 
                        c("DeleteMe" )]

# Using OG Data 
colnames(all_Data)[colnames(all_Data) == "MeanCP.Fescue"] ="OG_CP_Fescue"
colnames(all_Data)[colnames(all_Data) == "MeanCP.Epichloe"] ="OG_CP_Epichloe"
all_Data$Delta_CT_OG <- all_Data$OG_CP_Fescue - all_Data$OG_CP_Epichloe

# Making phenotype data to export
Data_to_export <- subset(all_Data, select = c("Treatment", "Delta_CT", "Data_Set", "Standard", "Maternal_Parent", "Delta_CT_OG"))
write.table(Data_to_export, file = '/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Phenotype_Data_Delta_CT.txt', sep = '\t', row.names=FALSE)

# Checking the amount of progeny each parent has
parents_one_row_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/all_used_parents_one_row.csv"
parents_one_row <- read.csv(parents_one_row_loc, header = FALSE, strip.white=TRUE)
parents_one_row %>% count(V1, sort = TRUE)

######
# Quick graph checks to make sure data is correct
######

#HIstograms with different coloring based on different groupings
plot1 <- ggplot(all_Data, aes(x=Delta_CT, color=Data_Set)) + 
  geom_histogram(binwidth=.2, fill = "white")

ggplot(all_Data, aes(x=Delta_CT, color=Maternal_Parent)) + 
  geom_histogram(binwidth=.2, fill = "white")

#Using the OG CP values (not efficiency adjusted)
plot2 <- ggplot(all_Data, aes(x=Delta_CT_OG, color=Data_Set)) + 
  geom_histogram(binwidth=.2, fill = "white")

ggarrange(plot2, plot1, ncol = 1, nrow = 2)

################################################################################
# Adding in the Alkaloid data
################################################################################

#Loading in data
alklaoid_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Alkaloid_Data_Combined.csv"
alklaoid <- read.csv(alklaoid_loc, header = TRUE, strip.white=TRUE)

# Giving data parents
alklaoid[c('Maternal_Parent', 'DeleteMe')] <- str_split_fixed(alklaoid$ID, '-', 2)
alklaoid <- subset(alklaoid, select = -c(DeleteMe))

# Removing outliars based on plate they were done on 
# Outliars removed are any samples 3x std away from the mean by each set. 
Data_Sets <- as.matrix(unique(alklaoid$Plate))
column_names <- colnames(alklaoid)
All_Data <- data.frame(matrix(ncol = length(column_names)))
colnames(All_Data) <- column_names

for (i in 1:nrow(Data_Sets)){
  
  Set <- Data_Sets[i]
  Current_Data <- subset(alklaoid, Plate == Set)
  
  #Removing outliars
  meanalki <- mean(Current_Data$ng.g)
  Sdalki <- sd(Current_Data$ng.g)
  cutoff1 <- meanalki + 3*Sdalki
  cutoff2 <- meanalki - 3*Sdalki
  Before <- nrow(Current_Data)
  Current_Data <- subset(Current_Data, ng.g < cutoff1 & ng.g > cutoff2)
  after <- nrow(Current_Data)
  
  print(paste0(Set, " Dropped ", Before-after))
  All_Data <- bind_rows(All_Data, Current_Data)
  }
alklaoid_filtered <- All_Data
alklaoid_filtered <- na.omit(alklaoid_filtered)
alklaoid_filtered_parents <- alklaoid_filtered
alklaoid_filtered <- alklaoid_filtered[-c(1307:1323), ]

alklaoid_filtered$Plate <- as.character(alklaoid_filtered$Plate)  

# Quick alkaloid data graph
ggplot(alklaoid_filtered, aes(x=ng.g, fill=Plate)) + 
  geom_histogram(binwidth=1000, color = "black")

ggplot(alklaoid_filtered, aes(x=ng.g, fill=Maternal_Parent)) + 
  geom_histogram(binwidth=1000, color = "black")  +
  ggtitle("Alkaloid Data Outliars Removed")


# Making letters for groups
# Making letters for box plot alkaloids
anovaA <- aov(ng.g ~ Maternal_Parent, data = alklaoid_filtered)
summary(anovaA)
tukeyA <- TukeyHSD(anovaA)
print(tukeyA)
TkA <- group_by(alklaoid_filtered, Maternal_Parent) %>%
  summarise(mean=mean(ng.g), quant = quantile(ng.g, probs = .75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anovaA, tukeyA) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent)
TkA$cld <- cld$Letters
print(TkA)

# Boxplots Alkaloids
ggplot(alklaoid_filtered, aes(x=Maternal_Parent, y=ng.g, fill = Maternal_Parent, group=Maternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Alkaloid") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = TkA, aes(label = cld, x = Maternal_Parent, y = quant), 
            vjust = 2.5 , size = 5)


################################################################################
# Merging all data together with MetaData
################################################################################
# Merging the data alkaloid and CT data together
CT_Values <- subset(all_Data, select = c(Treatment, Data_Set, Delta_CT, Delta_CT_OG, Fes_to_Epi_Ratio, Standard, Maternal_Parent))
colnames(CT_Values)[colnames(CT_Values) == "Treatment"] ="ID"
phenotype_Data <- merge(CT_Values, alklaoid, by.x = c("ID", "Maternal_Parent"), by.y = c("ID", "Maternal_Parent" ) ,all = TRUE)

colnames(phenotype_Data)[colnames(phenotype_Data) == "Plate"] ="Alkaloid_Plate"
phenotype_Data$Maternal_Parent <- as.character(phenotype_Data$Maternal_Parent)

# adding the metadata and merging it with the phenotype data
extraction_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Meta_Data/Extraction_info.csv"
harvest_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/Meta_Data/Harvest_Info.csv"
extraction <- read.csv(extraction_loc, header = TRUE, strip.white=TRUE)
harvest <- read.csv(harvest_loc, header = TRUE, strip.white=TRUE)

metadata <- merge(extraction, harvest, by.x = c("ID"), by.y = c("ID"), all = TRUE)
phenotype_Data <- merge(phenotype_Data, metadata, by = c("ID"))
phenotype_Data$Alkaloid_Plate <- as.character(phenotype_Data$Alkaloid_Plate)

# Saving the data to make a phenotype file
phenotype_Data$ng.g <- round(phenotype_Data$ng.g, 3)
write.table(phenotype_Data, file = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/phenotype_data.txt", row.names = FALSE)

################################################################################
# Gathering Residual Data
################################################################################
# getting rid of out liar data 
#plots will show us which data pints to turn into NAs
plot(phenotype_Data$Delta_CT_OG, col = as.factor(phenotype_Data$Extraction_Date))
plot(phenotype_Data$Delta_CT, col = as.factor(phenotype_Data$Extraction_Date))
plot(phenotype_Data$ng.g, col = as.factor(phenotype_Data$Extraction_Date))

# Removing other specific duplicates and outliars
phenotype_Data$Delta_CT[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "306-3-8"] <- NA
phenotype_Data$Delta_CT_OG[phenotype_Data$ID == "320-5-26"] <- NA
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "314" & phenotype_Data$Data_Set != "315x320"), ]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Data_Set == "315x320"),]
phenotype_Data <- phenotype_Data[!(phenotype_Data$ID == "315-1-8" & phenotype_Data$Extraction_Date == "03/16/23"),]

# This is to get the residuals of maternal parent while taking everything out that affects CT values
filtered_CT_model <- aov(Delta_CT ~ Harvest_Date + Standard + Extraction_Date + Extractor+ Data_Set, data = phenotype_Data)
summary(filtered_CT_model) # The filtered data 

raw_CT_model <- aov(Delta_CT_OG ~ Harvest_Date + Standard + Extraction_Date + Extractor+ Data_Set, data = phenotype_Data)
summary(raw_CT_model)




alkaloid_model <- lm(ng.g ~ Harvest_Date + Extractor + Extraction_Date + Alkaloid_Plate, data = phenotype_Data)
summary(alkaloid_model)

# extreacting all the residual data and combining it together 
Delta_CT_Data <- phenotype_Data[!is.na(phenotype_Data$Delta_CT), ]
Delta_CT_OG_Data <- phenotype_Data[!is.na(phenotype_Data$Delta_CT_OG), ]
Alkaloid_Data <- phenotype_Data[!is.na(phenotype_Data$ng.g), ]

filtered_CT_model <- lm(Delta_CT ~ Harvest_Date + Standard + Extraction_Date + Extractor + Data_Set, data = Delta_CT_Data)
raw_CT_model <- lm(Delta_CT_OG ~ Harvest_Date + Standard + Extraction_Date + Extractor + Data_Set, data = Delta_CT_OG_Data)
alkaloid_model <- lm(ng.g ~ Harvest_Date  + Alkaloid_Plate, data = Alkaloid_Data)

data1 <- data.frame(Delta_CT_Data$ID, filtered_CT_model$residuals)
data2 <- data.frame(Delta_CT_OG_Data$ID, raw_CT_model$residuals)
data3 <- data.frame(Alkaloid_Data$ID, alkaloid_model$residuals)

names(data1)[names(data1) == "Delta_CT_Data.ID"] <- "ID"
names(data2)[names(data2) == "Delta_CT_OG_Data.ID"] <- "ID"
names(data3)[names(data3) == "Alkaloid_Data.ID"] <- "ID"

residual_data <- merge(data1, data2, by = "ID", all = TRUE)
residual_data <- merge(residual_data, data3, by = "ID", all = TRUE)

# This data set is your final phenotype data
######################
residual_data
###################### 
plot(residual_data$filtered_CT_model.residuals)
plot(residual_data$raw_CT_model.residuals)
plot(residual_data$alkaloid_model.residuals)
























################################################################################
# no filtering analysis
################################################################################
#load in CT data 
ct_data_raw <- merge(Data_g3p4_Means_all,Data_2x2_Means_all, by.x = c("Treatment", "Data_Set"), by.y = c("Treatment", "Data_Set"))
names(ct_data_raw)[names(ct_data_raw) == "MeanCP.x"] <- "CP_Fescue"
names(ct_data_raw)[names(ct_data_raw) == "MeanCP.y"] <- "CP_Epichloe"
names(ct_data_raw)[names(ct_data_raw) == "Treatment"] <- "ID"
ct_data_raw[c('Maternal_Parent', 'DeleteMe')] <- str_split_fixed(ct_data_raw$ID, '-', 2)
ct_data_raw <- subset(ct_data_raw, select = -c(DeleteMe))
ct_data_raw$DeltaCT_Raw <- ct_data_raw$CP_Fescue - ct_data_raw$CP_Epichloe

# loading in alkaloid data
alklaoid_raw <- read.csv(alklaoid_loc, header = TRUE, strip.white=TRUE)
alklaoid_raw[c('Maternal_Parent', 'DeleteMe')] <- str_split_fixed(alklaoid_raw$ID, '-', 2)
alklaoid_raw <- subset(alklaoid_raw, select = -c(DeleteMe))
alklaoid_raw$Plate <- as.character(alklaoid_raw$Plate)
names(alklaoid_raw)[names(alklaoid_raw) == "Plate"] <- "Alkaloid_Plate"

# Merge all data together 
phenotype_Data_raw <- merge(ct_data_raw, alklaoid_raw, by.x = c("ID", "Maternal_Parent"), by.y = c("ID", "Maternal_Parent") ,all = TRUE)
phenotype_Data_raw <- merge(phenotype_Data_raw, metadata, by.x = c("ID"), by.y = c("ID") ,all = TRUE)
phenotype_Data_raw <- merge(phenotype_Data_raw,standard_group, by.x = c("Data_Set"), by.y =c("Data_Set"), all =  TRUE)


# Plots comparing raw to unraw (no outliar removal and no efficiency adjustements)
rawscatter <- ggplot(phenotype_Data_raw, aes(x = DeltaCT_Raw, y = ng.g, color = Maternal_Parent)) +
  geom_point()
cleanscatter <- ggplot(phenotype_Data, aes(x = Delta_CT, y = ng.g, color = Maternal_Parent)) +
  geom_point()

ggarrange(rawscatter,
          cleanscatter,
          nrow = 2,
          ncol = 1)

# Extracting residuals and putting them all together
Delta_CT_Data_raw <- phenotype_Data_raw[!is.na(phenotype_Data_raw$DeltaCT_Raw), ]
Alkaloid_Data_raw <- phenotype_Data_raw[!is.na(phenotype_Data_raw$ng.g), ]

CT_model_raw <- lm(DeltaCT_Raw ~ Harvest_Date + Standard + Extraction_Date + Extractor + Data_Set, data = Delta_CT_Data_raw)
alkaloid_model_raw <- lm(ng.g ~ Harvest_Date  + Alkaloid_Plate, data = Alkaloid_Data_raw)

data4 <- data.frame(Delta_CT_Data_raw$ID, CT_model_raw$residuals)
data5 <- data.frame(Alkaloid_Data_raw$ID, alkaloid_model_raw$residuals)

nrow(Alkaloid_Data_raw)
length(alkaloid_model_raw$residuals)

names(data4)[names(data4) == "Delta_CT_Data_raw.ID"] <- "ID"
names(data5)[names(data5) == "Alkaloid_Data_raw.ID"] <- "ID"

residual_data_raw <- merge(data4, data5, by = "ID", all = TRUE)

################################################################################
# making graphs to compare residual raw vs residual filterd
################################################################################
metadata <- subset(phenotype_Data, select = c("ID", "Maternal_Parent", "Alkaloid_Plate", "Extraction_Date", "Extractor", "Data_Set", "Standard", "Harvest_Date") )
residual_data_M <- merge(residual_data, metadata, by="ID")
residual_data_raw_M <- merge(residual_data_raw, metadata, by="ID")

# Scatter plot of delta CT residuals
CT_Data <- ggplot(residual_data_M, aes(x=ID, y=filtered_CT_model.residuals)) +
  geom_point() +
  ylim(-5,5) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Efficiency filtered CT values",
       x = "ID",
       y = "Delta CT")

CT_Data_raw <- ggplot(residual_data_raw_M, aes(x=ID, y=CT_model_raw.residuals)) +
  geom_point() +
  ylim(-5,5) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Raw Data",
       x = "ID",
       y = "Delta CT")
ggarrange(CT_Data, CT_Data_raw, ncol = 1, nrow = 2)

# Histogram plot of delta CT residuals
plot1 <- ggplot(residual_data_M, aes(x=filtered_CT_model.residuals)) +
  geom_histogram(bins = 60, fill = "blue", color = "black") +
  ylim(0,250) +
  xlim(-8,8) +
  labs(title = "Efficiency filtered CT values",
       x = "Delta CT Residuals",
       y = "Frequency")
plot2 <- ggplot(residual_data_M, aes(x=raw_CT_model.residuals)) +
  geom_histogram(bins = 60, fill = "green", color = "black") +
  ylim(0,250) +
  xlim(-8,8) +
  labs(title = "filtered raw CT values",
       x = "Delta CT Residuals",
       y = "Frequency")
plot3 <- ggplot(residual_data_raw_M, aes(x=CT_model_raw.residuals)) +
  geom_histogram(bins = 60, fill = "red", color = "black") +
  ylim(0,250) +
  xlim(-8,8) +
  labs(title = "Unfiltered raw CT values",
       x = "Delta CT Residuals",
       y = "Frequency")
ggarrange(plot1, plot2, plot3, ncol = 1, nrow = 3)

# Histogram plot of alkaloid residuals
plot2 <- ggplot(residual_data, aes(x=alkaloid_model.residuals)) +
  geom_histogram(bins = 60, fill = "blue", color = "black") +
  labs(title = "Filtered alkaloid values",
       x = "ng per gram",
       y = "Frequency")
plot3 <- ggplot(residual_data_raw, aes(x=alkaloid_model_raw.residuals)) +
  geom_histogram(bins = 60, fill = "red", color = "black") +
  labs(title = "Unfiltered raw alkaloid values",
       x = "ng per gram",
       y = "Frequency")
ggarrange(plot2, plot3, ncol = 1, nrow = 2)

# Scatter plot to see if we have groupings of stuff
plot5 <- ggplot(residual_data_M, aes(x=Data_Set, y=filtered_CT_model.residuals, color=Data_Set)) +
  geom_point() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Efficiency Filtered Delta CT Values",
       x = "ID",
       y = "Delta CT")
plot6 <- ggplot(residual_data_M, aes(x=Data_Set, y=raw_CT_model.residuals, color=Data_Set)) +
  geom_point() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Filtered Raw Delta CT Values",
       x = "ID",
       y = "Delta CT")
plot7 <- ggplot(residual_data_raw_M, aes(x=Data_Set, y=CT_model_raw.residuals, color=Data_Set)) +
  geom_point() +
  labs(title = "Unfiltered Delta CT Values",
       x = "ID",
       y = "Delta CT")
ggarrange(plot5, plot6, plot7, ncol = 1, nrow = 3)

# scatterplot plot of alkaloid residuals
plot7 <- ggplot(residual_data_M, aes(x=Alkaloid_Plate, y=alkaloid_model.residuals, color=Alkaloid_Plate)) +
  geom_point() +
  labs(title = "Filtered Alkaloid Values",
       x = "ID",
       y = "ng/g")
plot8 <- ggplot(residual_data_raw_M, aes(x=Alkaloid_Plate, y=alkaloid_model_raw.residuals, color=Alkaloid_Plate)) +
  geom_point() +
  labs(title = "Unfiltered Alkaloid Values",
       x = "ID",
       y = "ng/g")
ggarrange(plot7, plot8, ncol = 1, nrow = 2)

################################################################################
# Making graphs more usefull for presentations, with residuals
################################################################################
write.table(residual_data_M, file = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Phenotype_Data/All_Data_Filtered/residual_data.txt", row.names = FALSE)

# is the residual data linearly related
plot(residual_data$filtered_CT_model.residuals, residual_data$alkaloid_model.residuals) # fuck no

# Make the above graph just better
model <- lm(filtered_CT_model.residuals ~ alkaloid_model.residuals, data = residual_data_M)
# Extract the R-squared value
rsquared <- summary(model)$r.squared
# Create a scatter plot
ggplot(residual_data_M, aes(x = alkaloid_model.residuals, y = filtered_CT_model.residuals, color = Maternal_Parent)) +
  geom_point() +
  annotate("text", x = -12000, y = 2.4, label = paste("R-squared =", round(rsquared, 3))) +
  theme_bw() +
  labs(title = "Scatter Plot with R-squared Value",
       x = "Alkaloid Residuals",
       y = "Delta CT Residuals")
## THis next plot looks at the non residual data
# Make the above graph just better
model <- lm(Delta_CT ~ ng.g, data = phenotype_Data)
# Extract the R-squared value
rsquared <- summary(model)$r.squared
# Create a scatter plot
ggplot(phenotype_Data, aes(x = ng.g, y = Delta_CT)) +
  geom_point() +
  annotate("text", x = 35000, y = 5, label = paste("R-squared =", round(rsquared, 3))) +
  labs(title = "Scatter Plot with R-squared Value",
       x = "Alkaloid Residuals",
       y = "Delta CT Residuals")

# Showing us getting rid of terrible batch effects.
#Using the OG CP values (not efficiency adjusted)
plot1 <- ggplot(all_Data, aes(x=Delta_CT, color=Data_Set)) + 
  geom_histogram(binwidth=.2, fill = "white")
plot2 <- ggplot(residual_data_M, aes(x=filtered_CT_model.residuals, color=Data_Set)) + 
  geom_histogram(binwidth=.2, fill = "white") +
  xlab("Delta CT Residuals")
ggarrange(plot1, plot2, ncol = 1, nrow = 2)

# Showing Box plots by maternal parent
residual_data_M$Maternal_Parent <- as.factor(residual_data_M$Maternal_Parent)
removeems <- c("311","318","319")
residual_data_M_2 <- subset(residual_data_M, !(Maternal_Parent %in% removeems))
# Getting the letters for the plots
anova <- aov(alkaloid_model.residuals ~ Maternal_Parent, data = residual_data_M_2)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(residual_data_M_2, Maternal_Parent) %>%
  summarise(mean=mean(alkaloid_model.residuals), quant = quantile(alkaloid_model.residuals, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent)
Tk$cld <- cld$Letters
print(Tk)


### Making Boxplot ###

ggplot(residual_data_M_2, aes(x=Maternal_Parent, y=alkaloid_model.residuals, fill = Maternal_Parent, group=Maternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Alkaloid Residuals") +
  scale_fill_discrete(name = "Plant Lines") + 
  geom_text(data = Tk, aes(label = cld, x = Maternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)

# Getting the letters for the plots
anova <- aov(filtered_CT_model.residuals ~ Maternal_Parent, data = residual_data_M_2)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
Tk <- group_by(residual_data_M_2, Maternal_Parent) %>%
  summarise(mean=mean(filtered_CT_model.residuals), quant = quantile(filtered_CT_model.residuals, probs = .75, na.rm = TRUE)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(anova, tukey) #forgot why we do this
print(cld)
cld <- as.data.frame.list(cld$Maternal_Parent)
Tk$cld <- cld$Letters
print(Tk)
ggplot(residual_data_M_2, aes(x=Maternal_Parent, y=filtered_CT_model.residuals, fill = Maternal_Parent, group=Maternal_Parent)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
  theme_bw() +
  xlab("Plant Lines") +
  ylab("Delta CT residuals") +
  scale_fill_discrete(name = "Plant Lines") +
  geom_text(data = Tk, aes(label = cld, x = Maternal_Parent, y = quant), 
            vjust = -1.3, hjust = 1.1, size = 5)



gvlma(alkaloid_model)
gvlma(alkaloid_model_raw)














