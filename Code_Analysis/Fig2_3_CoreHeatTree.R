#This code creates Fig,2 and 3, which are heat trees produced by metacoder, depicting the fungal and bacterial core microbiomes. 
#Code preceding this include:
#CoreBacteria.R
#CoreFungi.R

library(tidyverse)
library(stringr)
library(magrittr)
library(phyloseq)
library(ggthemes)
library(ggplot2)
library(vegan)
library(metacoder)
library(taxa)
library(microbiome)

rm(list = ls()) #clear workspace

##################
# FUngi
##################

#two percent threshold 
core.ITS<- readRDS("./phy_object/phy_fun_core.rds")

#change names of columns
tax_table(core.ITS)[,"Genus"] %<>% gsub("g__","",.)
tax_table(core.ITS)[,"Family"] %<>% gsub("f__","",.)
tax_table(core.ITS)[,"Order"] %<>% gsub("o__","",.)
tax_table(core.ITS)[,"Class"] %<>% gsub("c__","",.)
tax_table(core.ITS)[,"Phylum"] %<>% gsub("p__","",.)
tax_table(core.ITS)[,"Kingdom"] %<>% gsub("k__","",.)

#get rid of species column 
tax <- tax_table(core.ITS) %>% data.frame
tax$Species <- NULL
tax %<>% as.matrix
tax_table(core.ITS) <- tax_table(tax)

core.ITS.coder <- parse_phyloseq(core.ITS)
core.ITS.coder$data$tax_abund <- calc_taxon_abund(core.ITS.coder, "otu_table", cols= core.ITS.coder$data$sample_data$sample_id, groups = rep("total", nrow(core.ITS.coder$data$sample_data))) 


heat_tree(core.ITS.coder,
          node_size = n_obs,
          node_color = total,
          edge_size_range = c(0.008, 0.008),
          node_label = taxon_names,
          #node_color_range = c("grey77", "steelblue2", "steelblue3", "steelblue4"),
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Read count",
          #make_node_legend = F,
          #initial_layout = "re",
          layout = "re" )
#layout types, re reingold-tilford (the default),  fr, and da

#ggsave("heattree_core_fungi_genus_re.jpeg", path= "./output/Figures/", width = 7, height = 6, units = "in")


##################
# BACTERIA
##################

core.16S<- readRDS("./phy_object/phy_bac_core.rds")

tax_table(core.16S)<- tax_table(core.16S)[ ,1:6] # this gets rid of the species column 

#gsub looks for something and replaces it. 
tax_table(core.16S)[,"Genus"] %<>% gsub("D_5__","",.)
tax_table(core.16S)[,"Family"] %<>% gsub("D_4__","",.)
tax_table(core.16S)[,"Order"] %<>% gsub("D_3__","",.)
tax_table(core.16S)[,"Class"] %<>% gsub("D_2__","",.)
tax_table(core.16S)[,"Phylum"] %<>% gsub("D_1__","",.)
tax_table(core.16S)[,"Kingdom"] %<>% gsub("D_0__","",.)

taxcheck <- tax_table(core.16S) %>% data.frame


core.16S.coder <- parse_phyloseq(core.16S)

core.16S.coder$data$tax_abund <- calc_taxon_abund(core.16S.coder, "otu_table", cols= core.16S.coder$data$sample_data$sample_id, groups = rep("total", nrow(core.16S.coder$data$sample_data))) 

print(core.16S.coder$data$tax_abund)



heat_tree(core.16S.coder,
                    node_size = n_obs,
                    node_color = total,
                    edge_size_range = c(0.008, 0.008),
                    node_label = taxon_names,
                    #node_color_range = c("grey77", "steelblue2", "steelblue3", "steelblue4"),
                    node_size_axis_label = "ASV count",
                    node_color_axis_label = "Read count",
                    #make_node_legend = F,
                    layout = "re")


#ggsave("heattree_core_bac.jpeg", path= "./output/Figures/", width = 7, height = 6, units = "in")

