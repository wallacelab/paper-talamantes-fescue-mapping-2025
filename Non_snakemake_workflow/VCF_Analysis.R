library(tidyverse)
install.packages("ggraph")
library(edgebundleR)
library(igraph)
library(ggraph)


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
hap_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Hap_Map/all_maf_5_150_min.hmp.txt"
hap <- read.table(hap_loc, header = FALSE, sep = '\t', skip = 1)
head(hap)

# Counting up chromosome appearances
variant_locations <- subset(hap, select = c(V1))
variant_locations[c('chrom', 'pos')] <- str_split_fixed(variant_locations$V1, '_', 2)
Chrom_counts <- variant_locations %>% count(chrom, sort = TRUE)


# graphs
ggplot(data=Chrom_counts, aes(x=reorder(chrom, -n), y=n)) +
  geom_bar(colour="grey", fill = "blue", stat = "identity") +
  xlab("Parental Combination") + ylab("Cross Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))



