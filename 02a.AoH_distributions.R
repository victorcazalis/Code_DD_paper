

### Charge species
if(Group=="FW_ODONATA"){species<-readRDS("1.Tables/Species.characteristics.Script1d.FW_ODONATA.rds")}
if(Group=="REPTILES"){species<-readRDS("1.Tables/Species.characteristics.Script1d.REPTILES.rds")}
if(Group=="AMPHIBIANS"){species<-readRDS("1.Tables/Species.characteristics.Script1d.AMPHIBIANS.rds")}
if(Group=="FISH"){species<-readRDS("1.Tables/Species.characteristics.Script1d.FISH.rds")}
if(Group=="MAMMALS"){species<-readRDS("1.Tables/Species.characteristics.Script1d.MAMMALS.rds")}


### Charge distributions: for Odonata I take the processed map (with the new polygons) and for others I re_run process without combining into a single polygon (to keep seasons etc)
source_lines <- function(file, lines){source(textConnection(readLines(file)[lines]))}

if(Group=="FW_ODONATA"){
  distributions<-st_read(paste0("0.Data/Distributions/processed/Dist.processed.", Group, ".shp"))
  distributions$presence<-distributions$origin<-distributions$seasonal<-1
  ### Change name of taxon identifier to fit with aoh package
  names(distributions)[names(distributions)=="binomil"]<-"binomial"
  distributions$SISID<-species$taxonid[match(distributions$binomial, species$scientific_name)]
  distributions$category<-species$category[match(distributions$binomial, species$scientific_name)]
  distributions$order<-"Odonata"
}

if(Group != "FW_ODONATA"){
  source_lines("2.Scripts/00.Process distributions.R", 4:23)  
}

sf::sf_use_s2(FALSE) 


### Subset distributions
speciesDD<-subset(species, species$category=="DD" &
                           species$Dist_exist==T &
                           (species$terrestrial_system==TRUE | species$freshwater_system==TRUE))

distDD<-subset(distributions, distributions$binomial %in% speciesDD$scientific_name)
rm(distributions, species)


### Project in Mollweide
CRSMOLL<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
distDD<-st_transform(distDD, CRSMOLL)
#ggplot(distDD)+geom_sf(fill="darkred", col=NA, alpha=0.2)+ggtitle(Group)

### Create empty columns
speciesDD$Forest2021<-speciesDD$ForestLoss<-speciesDD$AOHlost<-speciesDD$AOH2020<-NA


