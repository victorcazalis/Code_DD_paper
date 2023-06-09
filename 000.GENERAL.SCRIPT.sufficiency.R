

setwd("H:/Data sufficiency")
Main_dir<-"H:/"


library(plyr) ; library(dplyr) ; library(sf) ; library(rredlist) ; library(raster) ; library(exactextractr) ; library(rgbif) ; library(stringr) ; library(openxlsx) ; library(robis)
library(httr) ; library(jsonlite) ; library(fasterize)


# Create a function that works as the opposite of %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Create a function to report progress in loops
print_loop<-function(Index, Num){if((Index/Num)==round(Index/Num)){cat(Index)}}



#########################################
### RUN SOME DATA PREPARATION SCRIPTS ###
#########################################

# ### 1: Create a table with all species from the Red List
source("2.Scripts/00.Extract.species.R")

# Subset (all species I want in the analyses)
species<-readRDS("1.Tables/Species.merged.april22.rds")
species<-subset(species, (species$phylum=="CHORDATA" & species$class != "AVES") | species$order=="ODONATA")

# Assign distribution groups and RLA
source("2.Scripts/00.Extract.Groups.RLA.R")


### 2: Prepare SIS and GIS data
source("2.Scripts/00.Merge fish distributions.R") # Merges the fish distributions in 2 shapefiles
source("2.Scripts/00.Merge SIS data.R") # Merges data from SIS (habitat preferences, countries of occurrence...)
source("2.Scripts/00.Process distributions.R") # Processes distributions: filter polygons, rename columns, and edit distribution for Odonata
source("2.Scripts/00.Prepare.human.pop.rasters.R") # Prepares the raster of rural human population






######################################
### EXTRACT COVARIATES OF INTEREST ###
######################################


### 3: Calculate variables that concern all species
source("2.Scripts/01a.Extract.data.All.R")


### 4: Calculate variables with range size, per distribution folders
for(Group in c("FW_ODONATA", "AMPHIBIANS", "REPTILES", "MAMMALS", "FISH")){
  
  cat("      Starting the group: ", Group, "      ")
  
  species_raw<-read.csv("1.Tables/Species.characteristics.datasufficiency.Script1a.csv")
  species<-species_raw[species_raw$Dist_group==Group,]
  
  source("2.Scripts/01b.Extract.data.Range.R")
  if(TRUE %in% c(species$terrestrial_system, species$freshwater_system)){source("2.Scripts/01c.Extract.data.Terrestrial.R")} # If at least one species in terrestrial or freshwater
  if(TRUE %in% species$marine_system){source("2.Scripts/01d.Extract.data.Marine.R")} # If at least one marine species
  
  saveRDS(species, paste0("1.Tables/Species.characteristics.Script1d.", Group, ".rds"))

}







########################
### RUN AOH ANALYSES ###
########################

for(Group in c("MAMMALS", "REPTILES", "AMPHIBIANS", "FW_ODONATA")){

  source("2.Scripts/02a.AoH_distributions.R")
  source("2.Scripts/02b.AoH_CCI.R")
  source("2.Scripts/02c.AoH_forest.R")

}





################################
### RUN MAIN ANALYSES SCRIPT ###
################################
source("2.Scripts/03.Analysis.sufficiency.R")






####################################################
### COMPARE RESULTS WITH UPDATE RED LIST  2022.1 ###
####################################################
source("2.Scripts/04.Compare.reassessments.R")










