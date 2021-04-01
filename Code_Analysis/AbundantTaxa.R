#Code looks at percentage of read coverage by Phylum for bacteria and fungi, and also identifies the top most abundant OTUs

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)



############################################################
###FUNGI 


#start with raw phyloseq object. do minimal filtering to retain original read counts
phyITS <- readRDS("./phy_object/phy_fun_prefilter.rds")
phyITS %<>% subset_samples(., Seed.Lot != "Blank")
#made changes to meta file, needs to be replaced
newmetafun <- read.csv("C:/Users/delventk/Box/Research/Projects/Seedlot_Observational/2019 Combined/Meta/Meta_TareSoil_ITS.csv", row.names = 1)
sample_data(phyITS) <- sample_data(newmetafun)


phylum <- phyITS %>% tax_glom("Phylum") #This collapses everything to the phylum level first.
tax<- tax_table(phylum) %>% data.frame #create tax table
phylum_otu <- otu_table(phylum)%>% data.frame #create OTU table

phylum_otu %<>% t() %>% data.frame #flip matrix and convert back to data frame
phylum_otu$OTUSums <- rowSums(phylum_otu) #get sums of read counts for each OTU

otunames <- row.names(phylum_otu)#extract names
Sums <- phylum_otu$OTUSums #extract sums

df <- data.frame(otunames, Sums) #create data frame
df$total <- 8243407 #store total read count variable
df$percent <- (df$Sums/df$total)*100 #make new column expressing percent
df$tax <- tax$Phylum

#save the dataframe
#saveRDS(df, "./output/Rds/ITS_Phylum.rds")

###############
#top most abundant OTUs
rm(list = ls())

phyITS<- readRDS("./phy_object/phy_fun_nofilt.rds")
meta <- meta(phyITS)

topNOTUs.ITS <- names(sort(taxa_sums(phyITS),TRUE)[1:10])
phy10.ITS <- prune_taxa(topNOTUs.ITS, phyITS)
tax10.ITS <- tax_table(phy10.ITS) %>% data.frame

#write.csv(tax10.ITS, "./output/top10_ITS.csv")

############################################################
###BACTERIA

phybac <- readRDS("./phy_object/phy_TSbac_prefilter.rds")
phybac %<>% subset_samples(., Seed.Lot != "Blank")
newmetabac <- read.csv("C:/Users/delventk/Box/Research/Projects/Seedlot_Observational/2019 Combined/Meta/Meta_TareSoil_16S.csv", row.names = 1) #made changes to meta file, needs to be replaced
sample_data(phybac) <- sample_data(newmetabac)

phylum16 <- phybac %>% tax_glom("Phylum") #This collapses everything to the Phylum level
tax16<- tax_table(phylum16) %>% data.frame
phyotu16 <- otu_table(phylum16)%>% data.frame
phyotu16 %<>% t() %>% data.frame
phyotu16$OTUSums <- rowSums(phyotu16)

#extract for data frame
otunames16 <- row.names(phyotu16)
Sums16 <- phyotu16$OTUSums
#build data frame
df16 <-data.frame(otunames16, Sums16)
df16$total <- 7594024
df16$percent <- (df16$Sums16/df16$total)*100
df16$tax <- tax16$Phylum

#saveRDS(df16, "./output/Rds/16S_Phylum.rds")

###############
#top most abundant OTUs
rm(list = ls())

phy16s<- readRDS("./phy_object/phy_bac_nofilt.rds")

topNOTUs.16s <- names(sort(taxa_sums(phy16s),TRUE)[1:10])
phy10.16s <- prune_taxa(topNOTUs.16s, phy16s)
tax10.16s <- tax_table(phy10.16s) %>% data.frame

#write.csv(tax10.16s, "./output/top10_16s.csv")