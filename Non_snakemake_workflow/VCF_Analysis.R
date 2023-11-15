library(tidyverse)
install.packages("qqman")
library(edgebundleR)
library(igraph)
library(ggraph)
library(qqman)

# This makes manhattan plots

# Load in sample
depth = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/Variants_all.depth.txt"
depth <- read.table(depth, header = FALSE)
head(depth)


# Counting up the depth of samples
class(depth_counts$V3)
depth_counts$V3 <- as.numeric(depth_counts$V3)

depth_counts <- depth %>% count(V3, sort = TRUE)
depth_counts <- subset(depth_counts,(V3!="INDEL"))
depth_counts <- subset(depth_counts,(n>5))

depth_counts_zoomed <- subset(depth_counts,(V3 < 500))


# Graphs for depth
ggplot(data=depth_counts, aes(x=V3)) + 
  geom_histogram(colour="blue", fill="red", bins = 50) + 
  xlab("Depth") + ylab("Count") +
  ggtitle("Depth Counts")  

ggplot(data=depth_counts_zoomed, aes(x=V3)) + 
  geom_histogram(colour="blue", fill="red", bins = 40) + 
  xlab("Depth") + ylab("Count") +
  ggtitle("Depth Counts Zoomed")  

############################################################
## Cross Graphic maker ##
############################################################

# loading in list of predicted parents
pred_parents_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Lists/usable_predicted_parents_double.csv"
pred_parents <- read.table(pred_parents_loc, header = TRUE, sep = ",")
head(pred_parents)

# counting up unique parents
pred_parents$parent_combo <- 0
for (i in 1:nrow(pred_parents)){
  print(i)
  p1 = pred_parents[i,2]
  p2 = pred_parents[i,3]
  if (p1 > p2){
    pred_parents[i,4] = paste0(p2,"x",p1)
  }
  if (p2 > p1){
    pred_parents[i,4] = paste0(p1,"x",p2)
  }
}

parent_counts <- pred_parents %>% count(parent_combo, sort = TRUE)
parent_counts[c('Parent 1', 'Parent 2')] <- str_split_fixed(parent_counts$parent_combo, 'x', 2)

nodes <- t(cbind(t(pred_parents[,2]), t(pred_parents[,3])))
nodes <- data.frame(nodes) %>% group_by(nodes) %>%
  summarise(count=n())

# Graph
parent_counts2 <- subset(parent_counts, n >= 10)
ggplot(data=parent_counts2, aes(x=reorder(parent_combo, -n), y=n)) +
  geom_bar(colour="grey", fill = "blue", stat = "identity") +
  xlab("Parental Combination") + ylab("Cross Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1),
  panel.background = element_rect(fill = 'white', color = 'grey')) 

g <-graph.data.frame(pred_parents[,2:3], directed=F, vertices=nodes)
edgebundle( g )  

plot(g, vertex.size=10, edge.width=edge.betweenness(g))

#########################
# hap map analysis
##########################

#importing data and getting rid of file names
hap_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Hap_Map/all_maf_5_1960_min.hmp.txt"
hap <- read.table(hap_loc, header = FALSE, sep = '\t', skip = 1)
head(hap)

# Counting up chromosome appearances
variant_locations <- subset(hap, select = c(V1))
variant_locations[c('chrom', 'pos')] <- str_split_fixed(variant_locations$V1, '_', 2)
Chrom_counts <- variant_locations %>% count(chrom, sort = TRUE)


# graphs
Chrom_counts2 <- subset(Chrom_counts, n >= 10)
ggplot(data=Chrom_counts2, aes(x=reorder(chrom, -n), y=n)) +
  geom_bar(colour="grey", fill = "blue", stat = "identity") +
  xlab("Chromosome") + ylab("SNP Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1),
  panel.background = element_rect(fill = 'white', color = 'grey'))

############################
# Depth by site and individual
############################


data(mtcars)
dummy_data <- head(mtcars)
dummy_data_num_only <- subset(dummy_data, select = -c(mpg))
individuals <- ncol(dummy_data_num_only)
dummy_data_num_only$sum <- rowSums(dummy_data_num_only, )
dummy_data_num_only$count.0 <- apply(dummy_data_num_only, 1, function(x) length(which(x=="0")))
dummy_data_num_only$diviser <- individuals - dummy_data_num_only$count.0
dummy_data_num_only$Avg.Depth <- (dummy_data_num_only$sum / dummy_data_num_only$diviser)



ggplot(data=dummy_data_num_only, aes(x=Avg.Depth)) + 
  geom_histogram(colour="blue", fill="red", binwidth = 1) + 
  xlab("Depth") + ylab("Count") +
  ggtitle("Depth Counts")  

############################
# Manhattan plot in R
############################

mlm_stats_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Hap_Map/Linear_Models/MLM_16Pcs_stats.txt"
mlm_stats <- read.table(mlm_stats_loc, header = TRUE, sep = '\t')

mlm_stats$Chr<-gsub("PTG","",as.character(mlm_stats$Chr))
mlm_stats$Chr<-gsub("L","",as.character(mlm_stats$Chr))
mlm_stats$Chr<- as.numeric(mlm_stats$Chr)
mlm_stats <- mlm_stats[-c(1), ]

manhattan(mlm_stats, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          main = "P values of Delta CT", cex = 1, col = c("blue", "red","darkgrey","purple"))
# suggestiveline = F, genomewideline = F (We can take away lines with this) (1e-5, 5e-8)

# Put a line on graph
# geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed")

### The 16 pc GLM with 50 permutations
GLM_16_PCs_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Hap_Map/Linear_Models/GLM_16PC_50_Perm_Stats.txt"
GLM_16_PCs <- read.csv(GLM_16_PCs_loc, header = TRUE, strip.white=TRUE, sep = '\t')
GLM_16_PCs$Chr<-gsub("PTG","",as.character(GLM_16_PCs$Chr))
GLM_16_PCs$Chr<-gsub("L","",as.character(GLM_16_PCs$Chr))
GLM_16_PCs$Chr<- as.numeric(GLM_16_PCs$Chr)
GLM_16_PCs <- GLM_16_PCs[-c(1), ]
manhattan(GLM_16_PCs, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", 
          main = "P values of Delta CT", cex = 1, col = c("blue", "red","darkgrey","purple"))

### The 16 pc GLM with 1 permutations residule data
GLM_16_PCs_loc <- "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Hap_Map/Linear_Models/GLM_16PCs_1_perm_residules.txt"
GLM_16_PCs <- read.csv(GLM_16_PCs_loc, header = TRUE, strip.white=TRUE, sep = '\t')
# GLM_16_PCs <- GLM_16_PCs[GLM_16_PCs$Trait == "filtered_CT_model.residuals", ]
GLM_16_PCs <- GLM_16_PCs[GLM_16_PCs$Trait == "alkaloid_model.residuals", ]
GLM_16_PCs$Chr<-gsub("PTG","",as.character(GLM_16_PCs$Chr))
GLM_16_PCs$Chr<-gsub("L","",as.character(GLM_16_PCs$Chr))
GLM_16_PCs$Chr<- as.numeric(GLM_16_PCs$Chr)
GLM_16_PCs <- GLM_16_PCs[-c(1), ]
manhattan(GLM_16_PCs, chr = "Chr", bp = "Pos", p = "p", snp = "Marker", ylim = c(0, 4.2), 
          main = "P values of Alkaloid Data", cex = 1, col = c("blue", "red","darkgrey","purple")) 







