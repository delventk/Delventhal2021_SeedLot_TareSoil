#this code creates Fig3, which includes a map of samples locations, a small table listing the amount of samples within each seed grower, nested within states, and plots depicting jaccard community distances between samples against geographic distances between samples (km). 
#before running this code, you will need to run "SourceCode_Phylogeo.R"

library(geosphere)
library(rgeos)
library(microbiome)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

################## Create map

theme_set(theme_bw())

#Create map of sample locations 
phybac <- readRDS("./phy_object/phy_bac.rds")
phymap <- subset_samples(phybac, Latitude != "NA") #get rid of samples that don't have geographic coordinates
sample <- meta(phymap)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
sites<- data.frame(sample$Longitude, sample$Latitude)
(sites <- st_as_sf(sites, coords = c("sample.Longitude", "sample.Latitude"), 
                   crs = 4326, agr = "constant"))

#create map and store object
map<-  ggplot(data = world) +
  geom_sf() +
  geom_sf(data = sites, size = 3, shape = 21, fill = "black") +
  coord_sf(xlim = c(-135, -60), ylim = c(25, 60), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")


################## Test for distances
library(magrittr)
library(phylogeo)
library(geosphere)
library(vegan)
library(ggthemes)
#read in the unfiltered object

#load fungal "raw" phyloseq object (not hellinger transformed). Rarefy in order to measure community distance with Jaccard
phyfun <- readRDS("./phy_object/phy_fun_prefilter.rds")
phyfun %<>% subset_samples(., Seed.Lot != "Blank")
phyfun
phyfun %<>% subset_samples(Latitude != "NA")
phyfun %<>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
phyfun
phyfun %<>% rarefy_even_depth(rngseed = 711)

#load bacterial phyloseq, filter, rarefy
phybac <- readRDS("./phy_object/phy_bac_prefilter.rds")
phybac %<>% subset_samples(., Seed.Lot != "Blank")
phybac %<>% subset_samples(Latitude != "NA")
phybac %<>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
phybac %<>% rarefy_even_depth(rngseed = 711)

#creating the two distances to compare in a mantel test for fungi and then bacteria
distGeo_fun <- cbind(sample_data(phyfun)$Longitude, sample_data(phyfun)$Latitude) %>% distm %>% as.dist
distCom_fun <- phyfun %>% phyloseq::distance("jaccard", binary=T)

distGeo_bac <- cbind(sample_data(phybac)$Longitude, sample_data(phybac)$Latitude) %>% distm %>% as.dist
distCom_bac <- phybac %>% phyloseq::distance("jaccard", binary=T)

#jaccard distance - measure of how dissimilar two sets are, compared to actual (km) distance. jaccard will do presence absense, that is why is it important to rarefy to even depth versus running a proportional function--converitn to relative abundances (which you could do for Bray since it takes into account abundance)

#mantel test between distances
mantel(distGeo_fun,distCom_fun)
mantel(distGeo_bac,distCom_bac)

################## Create plot of distances

#you need to run code on "SourceCode_Phylogeo.R" first

FUNplot<- plot_distance(phyfun) +
  theme_few()+
  xlab("Geographic Distance (km)") +
  ylab("Jaccard Distance") +
  geom_point(size = 1.5) 

BACplot<- plot_distance(phybac) +
  theme_few()+
  xlab("Geographic Distance (km)") +
  ylab("Jaccard Distance") +
  geom_point(size = 1.5) 

################## Create plot of distances
meta <- sample_data(phybac) %>% data.frame
table <- meta$State %>% tapply(.,.,length) %>% data.frame
table$Number.of.Samples <- table$.
table$. <- NULL
table$Seed.Growers <- c(5,9, 7, 3, 1, 1, 5, 3, 1)
table$Samples <- table$Number.of.Samples

################## Combine
library(patchwork)

top<- map + gridExtra::tableGrob(table[,c('Seed.Growers', 'Samples')])
distance<- FUNplot + BACplot
#patchwork<- (map + plot_spacer())/ distance
patchwork <- top / distance
patchwork + plot_annotation(tag_levels = 'A')

#ggsave("Fig5_Distances.jpeg", path= "./output/Figures/", width = 8, height = 6, units = "in" )
