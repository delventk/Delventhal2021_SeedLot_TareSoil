#Pathogen profile

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(microbiome)

############################################ 
#fungi

phyITS<- readRDS("./phy_object/phy_fun_prefilter.rds")
phyITS %<>% subset_samples(., Seed.Lot != "Blank")

meta <- meta(phyITS)

#gsub looks for something and replaces it. Here we remove the first few characters of the column names for taxonomy table
tax_table(phyITS)[,"Species"] %<>% gsub("s__","",.)
tax_table(phyITS)[,"Genus"] %<>% gsub("g__","",.)
tax_table(phyITS)[,"Family"] %<>% gsub("f__","",.)
tax_table(phyITS)[,"Order"] %<>% gsub("o__","",.)
tax_table(phyITS)[,"Class"] %<>% gsub("c__","",.)
tax_table(phyITS)[,"Phylum"] %<>% gsub("p__","",.)
tax_table(phyITS)[,"Kingdom"] %<>% gsub("k__","",.)

#change "unidenfitied" to NA for uniformity  and better merging
tax_table(phyITS)[,] %<>% gsub("unidentified",NA,.)

taxonomy <- tax_table(phyITS) %>% data.frame #create tax table
phyglom <- tax_glom(phyITS, taxrank = "Genus", NArm = FALSE, bad_empty = FALSE) #group OTUs to Genus. this step is included because there are multiple OTUs associated with each pathogen
tax<- tax_table(phyglom) %>% data.frame 

#Check tax, identify the OTU number associated with each pathogen (KeepOTUs) and then subset
KeepOTUs <- c("OTU.7", "OTU.12", "OTU.24","OTU.47", "OTU.51","OTU.73",) #OTUs associated with putative pathogens
OTU <- otu_table(phyglom) %>%  t() %>% data.frame 
otupath <- subset(OTU, rownames(OTU) %in% KeepOTUs) #subset table to just putative pathogens
otupath$Sums <- rowSums(otupath) #get read counts across all samples for each OTU

#create parts to build a data.frame
otunames <- row.names(otupath)
Species<- c( "Thanatephorus","Fusarium", "Alternaria", "Colletotrichum", "Verticillium", "Helminthosporium")
reads <- otupath$Sums

df1 <- data.frame(otunames, Species, reads)


###to express it as a percent,we need to get the total number of reads in this dataset, which is outlined in the code below.
# otu_all<- OTU
# otu_all$SumOTUs <- rowSums(otu_all)
# otu_all %<>% t() %>% data.frame
# otu_all$SumAgain <- rowSums(otu_all)

#total reads are - 8243407
df1$percentreads <- (df1$reads/8243407)*100
df1$Kingdom <- "Fungi"

#to get the number of seed lots that the OTU occured in - we need to convert otupath to binary
otubinary <- otupath #create object to change
otubinary$Sums <- NULL
otubinary[otubinary > 0] <- 1
#then get a sum of the rows 
otubinary$Sums <- rowSums(otubinary)
SeedLots <- otubinary$Sums

#combine with data.frame
df1$SeedLots <- SeedLots
df1$occupancy <- (SeedLots/130)*100

#saveRDS(df1, "./output/Rds/PathogenTable_Fungi_genus.rds")

#repeat process with...
############################################ 
#bacteria 

phy16S<- readRDS("./phy_object/phy_TSbac_prefilter.rds")
phy16S %<>% subset_samples(., Seed.Lot != "Blank")

meta16 <- meta(phy16S)
taxcheck<- tax_table(phy16S) %>% data.frame
tax_table(phy16S)<- tax_table(phy16S)[ ,1:6] # this gets rid of the species column 

#gsub looks for something and replaces it. 
tax_table(phy16S)[,"Genus"] %<>% gsub("D_5__","",.)
tax_table(phy16S)[,"Family"] %<>% gsub("D_4__","",.)
tax_table(phy16S)[,"Order"] %<>% gsub("D_3__","",.)
tax_table(phy16S)[,"Class"] %<>% gsub("D_2__","",.)
tax_table(phy16S)[,"Phylum"] %<>% gsub("D_1__","",.)
tax_table(phy16S)[,"Kingdom"] %<>% gsub("D_0__","",.)
tax_table(phy16S)[,] %<>% gsub("unidentified",NA,.)

taxonomy16 <- tax_table(phy16S) %>% data.frame
phyglom16 <- tax_glom(phy16S, taxrank = "Genus", NArm = FALSE, bad_empty = FALSE)
tax16<- tax_table(phyglom16) %>% data.frame

#then identify OTUs of the plant pathogens and subset 
KeepOTU16s <- c("OTU.6","OTU.109","OTU.12964")
OTU16 <- otu_table(phyglom16) %>%  t() %>% data.frame
otupath16 <- subset(OTU16, rownames(OTU16) %in% KeepOTU16s)
otupath16$Sums <- rowSums(otupath16)

otu <- row.names(otupath16)
Genus<- c("Streptomyces","Pectobacterium", "Dickeya")

reads16 <- otupath16$Sums
df2 <- data.frame(otu, Genus, reads16)

# #to express it as a percent, get the total read count
# otu_ball<- OTU16
# otu_ball$SumOTUs <- rowSums(otu_ball)
# otu_ball %<>% t() %>% data.frame
# otu_ball$SumAgain <- rowSums(otu_ball)

#total reads are 7514851
df2$percentreads <- (df2$reads16/7514851)*100
df2$Kingdom <- "Bacteria"

#to get the number of seed lots each OTU occured in- we need to convert otupath to binary
otubinary16 <- otupath16
otubinary16$Sums <- NULL
otubinary16[otubinary16 > 0] <- 1
#then get a sum of the rows 
otubinary16$Sums <- rowSums(otubinary16)
SeedLots16 <- otubinary16$Sums

#add to dataframe
df2$SeedLots <- SeedLots16
df2$occupancy <- (SeedLots16/130)*100

#saveRDS(df2, "./output/Rds/PathogenTable_bacteria_genus.rds")

#merge dataframes
df<- merge(df1,df2, all.x = TRUE, all.y = TRUE)


#saveRDS(df, "./output/Rds/PathogenTable_full.rds")


#################################################################
# bargraph

library(ggplot2)
# Basic barplot
ggplot(data=df, aes(x=Taxon, y=Seed.Lots, color = Kingdom)) +
  geom_bar(stat="identity") + 
  theme_bw() +
  coord_flip()

#write.csv(df,"./output/CSV_Tables/PathogenTable.csv")
