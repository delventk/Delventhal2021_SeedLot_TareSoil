# Run this code before creating plots of community distance vs geographic (km) distances between are soil samples in the R code Fig5_DistanceCorrelation.R

#This source code is from the phylogeo package https://rdrr.io/github/zachcp/phylogeo/src/R/plot-methods.R
# the last function for plot_distance has been modified



setClass("phylogeo",
         representation(
           latitude = "character",
           longitude = "character"),
         contains = "phyloseq",
         prototype = prototype(otu_table = NULL,
                               tax_table = NULL,
                               sam_data = NULL,
                               phy_tree = NULL,
                               refseq = NULL,
                               latitude = NULL,
                               longitude = NULL)
)

#' initialize a phylogeo object from a phyloseq object
#'
#'
#' \code{phylogeo()} is a constructor method used internally to create a phylogeo
#' object and to ensure data integrity.
#'
#' @import phyloseq
#' @import dplyr
phylogeo <- function(physeq) {
  
  # check for sample_data slot
  if (!"sam_data" %in% phyloseq::getslots.phyloseq(physeq)) {
    stop("Mapping requires that phyloseq objects have Sample_Data with
              Latitude and Longitude")
  }
  
  # get lat/long columns. (based on https://github.com/rstudio/leaflet/blob/4ef0023c9fefa00a64e382ccd77d34d1413c47dc/R/normalize.R)
  sampledata <- data.frame(sample_data(physeq))
  sdfnames <- names(sampledata)
  lats = sdfnames[grep("^(lat|latitude)$", sdfnames, ignore.case = TRUE)]
  lngs = sdfnames[grep("^(lon|lng|long|longitude)$", sdfnames, ignore.case = TRUE)]
  if (!(length(lats) == 1 && length(lngs) == 1)) stop("Couldn't infer longitude/latitude columns")
  
  # Update sample_data and add it back to the phyloseq object
  sampledata <- sampledata %>%
    coerce_latlon_columns(lats) %>%
    coerce_latlon_columns(lngs)
  sample_data(physeq) <- sampledata
  
  # get samples without NA values in Lat/Lng columns
  samples_to_keep <- row.names(sampledata[ !is.na(sampledata[[lats]]) & !is.na(sampledata[[lngs]]), ])
  
  # if samples have been dropped prune OTUs belonging to them
  if (length(row.names) > 0 ) {
    #update sample data
    prune_samples(samples_to_keep, physeq)
    
    #update otu/sequence data
    if (.hasSlot(physeq, "otu_table")) {
      prune_taxa(taxa_sums(physeq) > 0, physeq)
    }
  }
  
  # return the class
  new('phylogeo',
      latitude = lats,
      longitude = lngs,
      otu_table = physeq@otu_table,
      tax_table = physeq@tax_table,
      sam_data = physeq@sam_data,
      phy_tree = physeq@phy_tree,
      refseq = physeq@refseq)
}

#' Address Lat//Long Columns:
#'  "None" -> NA
#'  factor -> vector
#'  vector -> numeric
#'
#' @keywords internal
coerce_latlon_columns <- function(df, col){
  colvals <- df[[col]]
  if (is.factor(colvals)) colvals <- as.vector(colvals) #convert to vector
  colvals[ colvals == "None"] <- NA  #some data has "None" so be sure to replace with NA
  df[[col]] <- as.numeric(as.character(colvals))
  df
}

##########################################################################
dist_to_edge_table = function(Dist, MaxDistance=NULL){
  dmat <- as.matrix(Dist)
  # Set duplicate entries and self-links to Inf
  dmat[upper.tri(dmat, diag = TRUE)] <- Inf
  distdf = reshape2::melt(dmat, as.is = TRUE)
  # Eliminate Inf Values (melt's third column is "value")
  distdf <- distdf[is.finite(distdf$value), ]
  names(distdf) <- c("Var1","Var2", "distance")
  # Remove entries above the threshold, MaxDistance
  if (!is.null(MaxDistance)) {
    distdf <- distdf[distdf$distance < MaxDistance, ]
  }
  return(distdf)
}

###########################################################################################
#THIS PLOT_DISTANCE FUNCTION HAS BEEN SLIGHTLY MODIFIED, some labeling of ggplot2 object is removed here and modified later

plot_distance <- function(physeq, distancemethod="jaccard"){
  #convert to phylogeo
  phygeo <- phylogeo(physeq)
  
  #get bigcircle distances using spDists
  #spDists expects the first column to be longitude
  df <- data.frame(sample_data(phygeo))
  df2 <- df[ c(phygeo@longitude, phygeo@latitude) ]
  names(df2) <- c("lon", "lat")
  df2 <- as.matrix(df2)
  geodistances <- spDists(df2, longlat = TRUE)
  colnames(geodistances)  <- row.names(df2)
  row.names(geodistances) <- row.names(df2)
  geodistances <- dist_to_edge_table(geodistances)
  names(geodistances) <- c("Var1","Var2","geodist")
  
  #get ecologicaldistances
  ecodistance <- distance(physeq, method = distancemethod)
  ecodistance <- dist_to_edge_table(ecodistance)
  names(ecodistance) <- c("Var1","Var2","ecodist")
  
  #make mergeable names for the two distance functions and merge
  concatvals <- function(x,y){ return(paste(x,"_",y,sep = ""))}
  geodistances['pairs'] <- mapply(concatvals, geodistances$Var1, geodistances$Var2)
  ecodistance['pairs'] <- mapply(concatvals, ecodistance$Var1, ecodistance$Var2)
  df <- merge(geodistances, ecodistance, by = "pairs")
  
  #make the plot
  p <- ggplot(df, aes(y = ecodist, x = geodist)) + geom_point()
  
  return(p)
}