## README description of R code for analyses 
Small note: Throughout the code, the proper name ASV (amplicon sequence variant) used in the manuscript is actually referred to here as OTU (operational taxonomic unit). 

AbundantTaxa.R is used to create a list of the most abundant ASVs and extract information on their ubiquity and read count. (corresponds to Table 1 in manuscript)

Core_Fungi.R and Core_Bacteria.R are the code used to determine the bacterial and fungal core microbiome of tare soil across the included seed lots. Visualizations of this process are available from the code FigS1_CoreThresholds.R. The core microbiome determined here is then used as input for Fig2_3_HeatTree.R, in which heat trees are created (corresponding to Figure 2 and 3 in the manuscript)

SourceCode_Phylogeo.R is modified source code from the R-package phylogeo, and must be run before running code for Fig5_DistanceCorrelation.R

PathogenProfile.R extracts information on OTUs that were assigned taxonomies of common potato pathogens.

All other R code files start with Fig, and correspond to the analysis associated with each Figure of the manuscript
