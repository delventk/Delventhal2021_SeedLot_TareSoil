#Figure 4 includes ordinations of tare soil microbial communities trough non-metric multidimensional scaling, grouped by seed grower source and state of origin. 
#Also included is code for the permutation analysis of variance (PERMANOVA) on the subset data

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(microbiome)
library(patchwork)
library(ggordiplots)

#clear workspace
rm(list = ls())

#read in objects
phyF <- readRDS("./phy_object/phy_fun.rds")
phyB <- readRDS("./phy_object/phy_bac.rds")

metaF <- meta(phyF)
metaB <- meta(phyB)

# Colorbind pallette
cbbPalette <- c("grey33", "#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2", "#D55E00", "#CC79A7", "tomato4")
# another pallette for State factor
cbbPalette2 <- c("darksalmon", "royalblue3", "green4", "#000000","gold")

theme_set(theme_bw())

################################
###FUNGI
################################

###############
###PERMANOVA

#first, subset so that there are only cultivars (i.e., varieties) included for which there are more than 2 instances or replicates. then, the same process repeats for the seed grower factor.
keep <- sample_data(phyF)$Variety %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1]
phy1 <- phyF %>% subset_samples(., Variety %in% keep) %>% prune_taxa(taxa_sums(.) > 0,.)

keep2 <- sample_data(phy1)$Seed.Grower %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1]
phy3 <- phy1 %>% subset_samples(., Seed.Grower %in% keep2) %>% prune_taxa(taxa_sums(.) > 0,.)

#after subsetting, there are two cultivars that are no longer represented more than twice, so we need to remove them. Same with one seed grower
phy3 %<>% subset_samples(Variety != "Snowden" & Variety != "Lamoka")
phy3 %<>% subset_samples(Seed.Grower != "Grower 12")#also get rid of seed grower 12

meta3 <- meta(phy3) #store meta data in an object
dist3<- phy3 %>% phyloseq::distance("bray") #community dissimilarity, bray-curtis

perm.fungi <- adonis2(dist3 ~ Year + Variety + Seed.Grower, data = meta3, by = "margin")
perm.fungi

###############
###ORDINATIONS

# ord <- metaMDS(dist3, 
#                distance= "bray", 
#                trymax = 250, 
#                maxit = 999,
#                engine = "monoMDS")
#saveRDS(ord, "./output/Rds/Ord_fungi.rds")

ordF<- readRDS("./output/Rds/Ord_fungi.rds")

my.plotF1<- gg_ordiplot(ordF, meta3$Seed.Grower, hull = TRUE, spiders = FALSE, ellipse = FALSE, plot = FALSE, label = FALSE)

plotF1 <- my.plotF1$plot

GrowerF<- plotF1 +
  theme(legend.position = "none") +
  scale_color_manual(values= cbbPalette) 


my.plotF2<- gg_ordiplot(ordF, meta3$State, hull = TRUE, spiders = FALSE, ellipse = FALSE, plot = FALSE, label = FALSE)

plotF2 <- my.plotF2$plot

StateF<- 
  plotF2 + 
  scale_color_manual(values= cbbPalette2) +
  theme(legend.position = "none")

GrowerF + StateF

################################
###BACTERIA
################################

###############
###PERMANOVA

#same filtering process above but repeated for bacteria
keepB <- sample_data(phyB)$Variety %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1]
phyB1 <- phyB %>% subset_samples(., Variety %in% keepB) %>% prune_taxa(taxa_sums(.) > 0,.)

keepB2 <- sample_data(phyB1)$Seed.Grower %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1]
phyB3 <- phyB1 %>% subset_samples(., Seed.Grower %in% keepB2) %>% prune_taxa(taxa_sums(.) > 0,.)

phyB3 %<>% subset_samples(Variety != "Snowden" & Variety != "Lamoka") #get rid of snowden and lamoka 
phyB3 %<>% subset_samples(Seed.Grower != "Grower 12")#also get rid of seed grower 12

metaB3 <- meta(phyB3)

distB3<- phyB3 %>% phyloseq::distance("bray") #distances 


perm.bacteria <- adonis2(distB3 ~ Year + Variety + Seed.Grower, data = metaB3, by = "margin")
perm.bacteria

###############
###ORDINATIONS

# ordB <- metaMDS(distB3, 
#                 distance= "bray",
#                 trymax = 250, 
#                 maxit = 999,
#                 engine = "monoMDS")
# #saveRDS(ordB, "./output/Rds/ord_bacteria.rds")

ordB<- readRDS("./output/Rds/ord_bacteria.rds")

my.plotB1<- gg_ordiplot(ordB, metaB3$Seed.Grower, hull = TRUE, spiders = FALSE, ellipse = FALSE, plot = FALSE, label = FALSE)

plotB1 <- my.plotB1$plot

GrowerB<- 
  plotB1 +
  labs(color = "Seed Grower")+
  scale_color_manual(values= cbbPalette) 


my.plotB2<- gg_ordiplot(ordB, metaB3$State, hull = TRUE, spiders = FALSE, ellipse = FALSE, plot = FALSE, label = FALSE)

plotB2 <- my.plotB2$plot

StateB<- 
  plotB2 + 
  labs(color = "State") +
  scale_color_manual(values= cbbPalette2) 


#####################
## create the plot

plot <- (GrowerF + GrowerB) / (StateF + StateB)
plot + plot_annotation(tag_levels = 'A')

#ggsave("NMS_SeedState.jpeg", path = "./output/Figures/", width = 7, height = 6, units = "in")
