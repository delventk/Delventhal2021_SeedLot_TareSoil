#Code is adapted from Shade and Stopnisek 2019, which uses occupancy abundance curves to define the core microbiome https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro

#BACTERIA

library(magrittr)
library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(microbiome)

rm(list = ls())


phy <- readRDS("./phy_object/phy_bac.rds")
otu <- otu_table(phy) %>% data.frame() %>% t()
map <- sample_data(phy) %>% data.frame()

check <- otu_table(phy) %>% data.frame

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- rowMeans(otu)    # mean relative abundance
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame


otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank=otu_occ) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start=otu_ranked$otu[1]                   
start_matrix <- as.matrix(otu[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(sum(otu[,x])))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
# calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. Can be set to the entire length of OTUs in the dataset, however it might take some time if more than 5000 OTUs are included.
for(i in 2:500){                              
  otu_add=otu_ranked$otu[i]                       
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(sum(otu[,x])))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}
# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(sum(otu[,x])))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs

Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)


elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])))
lastCall5 <-last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.05)])))
lastCall3 <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))

#Creating plot of Bray-Curtis similarity
ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept= lastCall, lty=3, col='blue', cex=.5) +
  geom_vline(xintercept= lastCall3, lty=3, col='darkgreen', cex=.5) +
  geom_vline(xintercept= lastCall5, lty=3, col= 'black', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+ 12, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall + 8, y=.425, label=paste("Last 2% increase (",lastCall,")",sep=''), color="blue")+
  annotate(geom="text", x=lastCall3 + 8, y=.35, label=paste("Last 3% increase (",lastCall3,")",sep=''), color="darkgreen")+
  annotate(geom="text", x=lastCall5 +8, y=.25, label=paste("Last 5% increase (",lastCall5,")",sep=''), color="black")

#ggsave("Core_Bac_Thresholds.jpeg", path = "./output/Figures/", width = 7.5, height = 4.5, units = "in")

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:7]] <- 'elbow'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[8:last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))]] <- 'core3'


occ_abun$core <- 'no'
occ_abun$core[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))]] <- 'core'


KeepOTUs <- occ_abun$otu[occ_abun$core == "core"]
otucore <- subset(otu, rownames(otu) %in% KeepOTUs)
otu_table(phy) <- otu_table(otucore, taxa_are_rows = TRUE)
phy

saveRDS(phy, "./phy_object/phy_bac_core.rds")

#saveRDS(occ_abun, "./output/Rds/occ_abun_bacteria.rds")
#saveRDS(BC_ranked, "./output/Rds/BC_ranked_bacteria.rds")
