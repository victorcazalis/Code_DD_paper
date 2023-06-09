

### Subset the distributions with only species occurring on land
distributionsMAR<-subset(distributions, distributions$binomial %in% species$scientific_name[species$marine_system==TRUE])



### Extract realms

# Rasterize realms
realm_mar<-read_sf("0.Data/Ocean_realms/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012.shp")
realm_mar<-realm_mar %>% dplyr::group_by(REALM) %>% dplyr::summarise(N= n())
st_write(realm_mar, "0.Data/Ocean_realms/Realms_oceans_COMBINED_by_Victor.shp")
realm_mar<-read_sf("0.Data/Ocean_realms/Realms_oceans_COMBINED_by_Victor.shp")
realm_mar<-st_transform(realm_mar, crs(distributions))

realm_mar$REALM_qtt<-revalue(realm_mar$REALM, c("Arctic"='1', "Atlantic Warm Water"='2', "Central Indo-Pacific"='3', "Eastern Indo-Pacific"='4', "Indo-Pacific Warm Water"='5', "Northern Cold Water"='6', "Southern Cold Water"='7', "Southern Ocean"='8', "Temperate Australasia"='9', "Temperate Northern Atlantic"='10', "Temperate Northern Pacific"='11', "Temperate South America"='12', "Temperate Southern Africa"='13', "Tropical Atlantic"='14', "Tropical Eastern Pacific"='15', "Western Indo-Pacific"="16")) %>% as.numeric(.)
realm_rastMAR<-fasterize(realm_mar, raster(extent(realm_mar), res=c(10000,10000), crs = crs(realm_mar)), field="REALM_qtt", fun="first")
writeRaster(realm_rastMAR, "0.Data/Ocean_realms/Marine_realms_rasterised.tif")
realm_rastMAR<-raster("0.Data/Ocean_realms/Marine_realms_rasterised.tif") %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), method="ngb")
distributionsMAR$Realms<-NA

# Extract
distributionsMAR[, c("realmsMAR_mode", "realmsMAR_var")]<-exact_extract(realm_rastMAR, distributionsMAR, c("mode", "variety"))
species$realmsMAR_mode<-distributionsMAR$realmsMAR_mode[match(species$scientific_name, distributionsMAR$binomial)]
species$realmsMAR_mode<-revalue(as.character(species$realmsMAR_mode), c("1"="Arctic", '2'="Atlantic Warm Water", '3'="Central Indo-Pacific", '4'="Eastern Indo-Pacific", '5'="Indo-Pacific Warm Water", '6'="Northern Cold Water", '7'="Southern Cold Water", '8'="Southern Ocean", '9'="Temperate Australasia", '10'="Temperate Northern Atlantic", '11'="Temperate Northern Pacific", '12'="Temperate South America", '13'="Temperate Southern Africa", '14'="Tropical Atlantic", '15'="Tropical Eastern Pacific", '16'="Western Indo-Pacific"))
species$realmsMAR_var<-distributionsMAR$realmsMAR_var[match(species$scientific_name, distributionsMAR$binomial)]

# For the few species with NA (eg a small island not mapped in REALMS) get the closest 
for(i in which(is.na(species$realmsMAR_mode) & species$Dist_exist=="TRUE" & species$marine_system==T)){
  species$realmsMAR_var[i]=1
  
  dist<-st_distance(distributionsMAR[as.character(distributionsMAR$binomial)==as.character(species$scientific_name[i]),], realm_mar)
  species$realmsMAR_mode[i]=which(dist==min(dist))
}

write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1d(TEMPO1", Group, ").csv"), row.names=FALSE)






### Distance to ports (marine accessibility)
if("acc_mar" %not in% ls()){ # Charge only if not already charged
  acc_mar<-raster("0.Data/Port distance/port_distance.asc") 
  acc_mar<-replace(acc_mar, acc_mar==0, NA)
  acc_mar<- acc_mar %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
  }

distributionsMAR[, c("acc_marMED", "acc_marQ10", "acc_marQ90")]<-exact_extract(acc_mar, distributionsMAR, "quantile", quantiles=c(0.5, 0.1, 0.9))

species$acc_marMED<-distributionsMAR$acc_marMED[match(species$scientific_name, distributionsMAR$binomial)]
species$acc_marQ10<-distributionsMAR$acc_marQ10[match(species$scientific_name, distributionsMAR$binomial)]
species$acc_marQ90<-distributionsMAR$acc_marQ90[match(species$scientific_name, distributionsMAR$binomial)]




### Fishing intensity
fishing_raw<-raster("0.Data/Fishing_intensity/FishingEffort_log2015_5km.tif") 
fishing_raw<-exp(fishing_raw) # I have to back transform to add zeros
fishing_raw<-replace(fishing_raw, is.na(fishing_raw), 0) # I replace NAs by 0
fishingPROJ<-projectRaster(fishing_raw, acc_mar) # I reproject to use the map / continent map from accessibility map
fishing<-fishingPROJ * replace(acc_mar, is.na(acc_mar)==F, 1) # Transform terrestrial areas in NA
fishing<-log(fishing+1) # Log transform

distributionsMAR[, c("fishingMED", "fishingQ90")]<-exact_extract(fishing, distributionsMAR, "quantile", quantiles=c(0.5, 0.9))
species$fishingMED<-distributionsMAR$fishingMED[match(species$scientific_name, distributionsMAR$binomial)]
species$fishingQ90<-distributionsMAR$fishingQ90[match(species$scientific_name, distributionsMAR$binomial)]





write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1d(TEMPO2", Group, ").csv"), row.names=FALSE)



