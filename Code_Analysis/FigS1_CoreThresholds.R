#Code of creating Supplementary Figure S1. Code is adapted from Shade and Stopnisek 2019, and the saved objects used here are created in more detail in the CoreBacteria.R and CoreFungi.R files of this project.

library(patchwork)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(reshape2)

##################
# FUNGI
##################

#read in objects
occ_abun_bac <- readRDS("./output/Rds/occ_abun_bacteria.rds")
BC_ranked_bac <- readRDS("./output/Rds/BC_ranked_bacteria.rds")

#A) elbow method
elbowbac <- which.max(BC_ranked_bac$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall3_bac <- last(as.numeric(as.character(BC_ranked_bac$rank[(BC_ranked_bac$IncreaseBC>=1.03)])))

#Creating plot of Bray-Curtis similarity
bacteria <- 
  ggplot(BC_ranked_bac[1:100,], aes(x=factor(rank[1:100], levels=rank[1:100]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbowbac, lty=3, col='red', cex=.5) +
  geom_vline(xintercept= lastCall3_bac, lty=3, col='blue', cex=.5) +
  labs(x='Ranked Bacterial ASVs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbowbac+ 6, y=.125, label=paste("Elbow method"," (",elbowbac,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall3_bac +7, y=.075, label=paste("Last 3% increase (",lastCall3_bac,")",sep=''), color="blue")


##################
# BACTERIA
##################

#read in objects
occ_abun_fun <- readRDS("./output/Rds/occ_abun_fungi.rds")
BC_ranked_fun <- readRDS("./output/Rds/BC_ranked_fungi.rds")


elbowfun <- which.max(BC_ranked_fun$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall3_fun <- last(as.numeric(as.character(BC_ranked_fun$rank[(BC_ranked_fun$IncreaseBC>=1.03)])))


fungi <- ggplot(BC_ranked_fun[1:100,], aes(x=factor(rank[1:100], levels=rank[1:100]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbowfun, lty=3, col='red', cex=.5) +
  geom_vline(xintercept= lastCall3_fun, lty=3, col='blue', cex=.5) +
  labs(x='Ranked Fungal ASVs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbowfun+ 6, y=.4, label=paste("Elbow method"," (",elbowfun,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall3_fun + 7, y=.15, label=paste("Last 3% increase (",lastCall3_fun,")",sep=''), color="blue")


fungi/bacteria + plot_annotation(tag_levels = 'A')


#ggsave("Core_Thresholds.jpeg", path = "./output/Figures/", width = 7, height = 5, units = "in")
