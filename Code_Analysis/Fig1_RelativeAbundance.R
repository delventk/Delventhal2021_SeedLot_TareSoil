#Figure 1 - Bar plots depicting the relative abundances of fungal and bacterial Orders grouped by seed grower origin and state grown.

##load libraries
library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(microbiome)
#library(goeveg)
library(patchwork)

############################################################
###FUNGI 

phyITS<- readRDS("./phy_object/phy_fun.rds") #read in the phyloseq object

##Filtering samples 
phyITS %<>% subset_samples(., Seed.Grower != "") #remove samples that do not have seed grower source 
KEEPseedg <- sample_data(phyITS)$Seed.Grower %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1] #include only seed growers that were present at least twice
phyITS %<>% subset_samples(Seed.Grower %in% KEEPseedg)
metaITS<- meta(phyITS) #create meta file for reference

##merge samples by seed grower
phy_mergeITS <- merge_samples(phyITS, "Seed.Grower") 
tax <- tax_table(phy_mergeITS) %>% data.frame

#covert to proportional (relative) abundance
proportional <- function(phy.obj){transform_sample_counts(phy.obj,function(x) {x / sum(x) })}
phy_mergeITS %<>% proportional()

#aggregate to only the top 15 orders
phytopITS <- aggregate_top_taxa(phy_mergeITS, 15, "Order")
tax_table(phytopITS)[,"Order"] %<>% gsub("o__","", .) #gsub looks for something and replaces it. 

metatop <- meta(phytopITS) #make meta file
metatop$State <- c("OR", "AB", "ID", "AB", "OR", "WA", "NE", "OR", "ID") #replace state names for each seed grower
sample_data(phytopITS) <- sample_data(metatop)

#colorpallete:
cpal<- c("#CC6666", "#F0E442", "#66CC99", "chocolate","purple", "lightgoldenrod3", "#9999CC", "steelblue1", "#CC79A7", "wheat","goldenrod", "grey66","steelblue4", "olivedrab3", "#009E73", "orange")

fungi<- plot_bar(phytopITS, fill = "Order") + 
  theme_classic()+
  geom_bar(aes(color= Order, fill= Order), stat= "identity", position= "stack") +
  facet_grid(~ State, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_manual(values= cpal) +
  scale_color_manual(values= cpal) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        axis.text.x = element_text(angle = 90, hjust = 0 , vjust=0.5)) +
  xlab("Seed Grower") + ylab("Relative Abundance") +
  ggtitle("Fungi")


############################################################
###Bacteria

phy16s<- readRDS("./phy_object/phy_bac.rds")

#Filtering samples 
phy16s %<>% subset_samples(., Seed.Grower != "")
KEEP <- sample_data(phy16s)$Seed.Grower %>% tapply(.,.,length) %>% .[is_greater_than(.,2)] %>% names %>% .[-1]
phy16s %<>% subset_samples(Seed.Grower %in% KEEP)

meta16s- meta(phy16s)

#merge
phy_merge16s <- merge_samples(phy16s, "Seed.Grower") 
tax <- tax_table(phy_merge16s) %>% data.frame

#proportional abundance
phy_merge16s %<>% proportional()

#aggregate
phytop16s <- aggregate_top_taxa(phy_merge16s, 15, "Order")
tax_table(phytop16s)[,"Order"] %<>% gsub("D_3__","", .)

metatop2 <- meta(phytop16s)
metatop2$State <- c("OR", "AB", "ID", "AB", "OR", "WA", "NE", "OR", "ID")
sample_data(phytop16s) <- sample_data(metatop2)

#colorpallete:
cpal2<- c("#CC6666", "#F0E442", "#66CC99", "chocolate","purple", "lightgoldenrod3", "#9999CC", "steelblue1", "#CC79A7", "wheat","goldenrod", "grey66","steelblue4", "olivedrab3", "orange", "#009E73")

bacteria<- plot_bar(phytop16s, fill = "Order") + 
  theme_classic()+
  geom_bar(aes(color= Order, fill= Order), stat= "identity", position= "stack") +
  facet_grid(~ State, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_manual(values= cpal2) +
  scale_color_manual(values= cpal2) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        axis.text.x = element_text(angle = 90, hjust = 0 , vjust=0.5)) +
  xlab("Seed Grower") + ylab("Relative Abundance") +
  ggtitle("Bacteria")



fungi + bacteria

ggsave("Fig1_RelativeAbundances.jpeg", path= "./output/Figures/Manuscript/", width = 9, height = 5.5, units = "in")
