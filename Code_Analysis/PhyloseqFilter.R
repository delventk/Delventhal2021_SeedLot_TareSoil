#This code includes all of the data processing and filtering steps to create the phyloseq objects that were used for data analysis

#load libraries
library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(microbiome)
library(goeveg)

#clear R workspace environment
rm(list = ls())


##################
# FUNGI
##################

phyfun <- readRDS("./phy_object/phy_fun_prefilter.rds")

#made changes to meta file, replaced here
newmetafun <- read.csv("./phy_object/Meta_ITS.csv", row.names = 1)
sample_data(phyfun) <- sample_data(newmetafun)
check<- meta(phyfun)

phyfun %<>% subset_samples(., Seed.Lot != "Blank") #removes controls and blanks

###this is the phy object that has been hellinger transformed-- square root of relative abundances
phyfun %<>% transform_sample_counts(function(x) sqrt(x/sum(x)))
saveRDS(phyfun, "./phy_object/phy_fun_nofilt.rds")

#filtering step to reduce noise. filters to retains samples that are present in at least 1% of the dataset
phyfun %<>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
saveRDS(phyfun, "./phy_object/phy_fun.rds")

##################
# BACTERIA
##################

#load in raw phyloseq object 
phybac <- readRDS("./phy_object/phy_bac_prefilter.rds")

#made changes to meta file, needs to be replaced
newmetabac <- read.csv("./phy_object/Meta_16S.csv", row.names = 1)
sample_data(phybac) <- sample_data(newmetabac)

phybac %<>% subset_samples(., Seed.Lot != "Blank")

phybac %<>% transform_sample_counts(function(x) sqrt(x/sum(x)))
saveRDS(phybac, "./phy_object/phy_bac_nofilt.rds")

phybac %<>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
phybac
saveRDS(phybac, "./phy_object/phy_bac.rds")



##################
# Some useful code to check phy objects and check information associated with the meta files 
##################
rm(list = ls())

phyfun <- readRDS("./phy_object/phy_fun.rds")
phybac <- readRDS("./phy_object/phy_bac.rds")

metafun<- meta(phyfun)
metabac<- meta(phybac)

length(unique(metafun$State))
length(unique(metabac$State))

metafun$State %>% tapply(.,.,length)
metabac$State %>% tapply(.,.,length)

length(unique(metafun$Variety))
metafun$Variety %>% tapply(.,.,length)

length(unique(metafun$Seed.Grower))
metafun$Seed.Grower %>% tapply(.,.,length)
