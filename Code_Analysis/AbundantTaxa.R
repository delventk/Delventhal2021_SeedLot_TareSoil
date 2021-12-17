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
#meta file needs to be replaced
newmetafun <- read.csv("./Meta/Meta_ITS.csv", row.names = 1)
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

#Read in object
phyITS<- readRDS("./phy_object/phy_fun.rds")
tax <- tax_table(phyITS) %>% data.frame # use tax table to assign taxonomy to top OTUs

topNOTUs.ITS <- names(sort(taxa_sums(phyITS),TRUE)) %>% data.frame

############################################################
###BACTERIA

phybac <- readRDS("./phy_object/phy_TSbac_prefilter.rds")
phybac %<>% subset_samples(., Seed.Lot != "Blank")
newmetabac <- read.csv("./Meta/Meta_16S.csv", row.names = 1) #made changes to meta file, needs to be replaced
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
###BACTERIA
rm(list = ls())

phybac<- readRDS("./phy_object/phy_bac.rds")
otu <- otu_table(phybac) %>% data.frame
tax <- tax_table(phybac) %>% data.frame # use tax table to assign taxonomy to top OTUs
topNOTUs.bac <- names(sort(taxa_sums(phybac),TRUE)) %>% data.frame

#OR curve of relative abundance
library(goeveg)
curve<- racurve(otu)
relabun <-curve$rel.abund %>% data.frame
