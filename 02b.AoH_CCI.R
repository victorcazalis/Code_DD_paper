

library(aoh)
library(terra)
library(rappdirs)
library(ggplot2)
library(raster)
library(ncdf4)
library(exactextractr)

cache_dir <- user_data_dir("aoh")
if (!file.exists(cache_dir)) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
}

output_dir <- "1.Tables/AOH"

is_iucn_rl_api_available()


### Charge CCI and altitude data
cci2010<-rast(paste0(Main_dir, "Postdoc/Platform/Rasters Mollweide/CCI2010_reprojMollweide.tif")) ; crs(cci2010)<-CRSMOLL
cci2020<-rast(paste0(Main_dir, "Postdoc/Platform/Rasters Mollweide/CCI2020_reprojMollweide.tif")) ; crs(cci2020)<-CRSMOLL
alt<-rast(paste0(Main_dir, "Postdoc/Platform/Rasters Mollweide/Elevation_MollweidePlatform.tif")) ; crs(alt)<-CRSMOLL
crosswalk_to_use<- read.csv("0.Data/Crosswalk_CCI_IUCN_Lumbierres.csv") ; crosswalk_to_use$code<-as.character(crosswalk_to_use$code)



### Add generation length information
GL<-read.csv("0.Data/Generation_length_sRedList.csv")
speciesDD$GL<-GL$GL_estimate[match(speciesDD$scientific_name, GL$internal_taxon_name)]


### Empty folders
unlink(output_dir, recursive=T)
dir.create(paste0(output_dir, "/Initial"), recursive=T) ; dir.create(paste0(output_dir, "/Current")) ; dir.create(paste0(output_dir, "/Temporary"))
terraOptions(tempdir=paste0(output_dir, "/Temporary"))




### Start species by species analyses (I have to loop to use the different initial CCI layers)

for(SP in 1:nrow(speciesDD)){
  tryCatch({
  
  SP_name<-speciesDD$scientific_name[SP]

  ### Clean range
  rangeSP_clean<-create_spp_info_data(distDD[distDD$binomial==SP_name,], 
                                    keep_iucn_rl_presence = 1:2, 
                                    keep_iucn_rl_seasonal = 1:2, 
                                    keep_iucn_rl_origin = c(1,2,5,6),
                                    crs=st_crs(CRSMOLL))
  
  ### Prepare altitude and CCIs
  # Altitude
  alt_crop<-crop(alt, extent(rangeSP_clean))
  
  # CCI 2020
  cci2_crop<-crop(cci2020, extent(rangeSP_clean))
  
  # CCI old
  GL_species<-speciesDD$GL[SP]
  if(is.na(GL_species)){GL_species<-1}
  
  Year1<-(2020-max(10, round(3*GL_species))) %>% max(., 1992)
  if(Year1==2010){cci1<-cci2010} else {
  cci1<-rast(sub("XXXX", Year1, paste0(Main_dir, "Postdoc/Platform/Rasters Mollweide/CCIXXXX_reprojMollweide.tif"))) ; crs(cci1)<-CRSMOLL} # I ensure the CRS is correctly assigned
  
  # Crop CCI1
  cci1_crop<-crop(cci1, extent(rangeSP_clean))
  
  
  ### Calculate AOH1
  AOH1=create_spp_aoh_data(
    rangeSP_clean, 
    elevation_data = alt_crop,
    habitat_data = cci1_crop,
    crosswalk_data = crosswalk_to_use,
    output_dir=paste0(output_dir, "/Initial"),
    engine="terra")$path %>% terra::rast(.)
  
  
  ### Calculate AOH2
  AOH2=create_spp_aoh_data(
      rangeSP_clean, 
      elevation_data = alt_crop,
      habitat_data = cci2_crop,
      crosswalk_data = crosswalk_to_use,
      output_dir=paste0(output_dir, "/Current"),
      engine="terra")$path %>% terra::rast(.)
    
  
  ### Calculate trends
  AOH1_area<-global(AOH1, "sum", na.rm=T) %>% as.numeric(.)
  AOH2_area<-global(AOH2, "sum", na.rm=T) %>% as.numeric(.)
  speciesDD$AOH2020[SP]<- AOH2_area 
  speciesDD$AOHlost[SP]<- (AOH1_area-AOH2_area)/(AOH1_area)
  
  } ,error=function(e){cat("\n Error (likely no habitat suitable in crosswalk) \n")})
  
  if((SP/10)==round(SP/10)){logger::log_info(SP)}
}


write.csv(speciesDD, paste0("1.Tables/AOH.", Group, "(TEMPO).csv"), row.names=F)
