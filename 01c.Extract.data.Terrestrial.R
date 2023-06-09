

### Subset the distributions with only species occurring on land
distributionsTERR<-subset(distributions, distributions$binomial %in% species$scientific_name[species$terrestrial_system==TRUE | species$freshwater_system==TRUE])


cat(paste("\n", "\n", "Starting realms", Sys.time(), "\n"))



### Extract realms 

# Rasterize realms
realm_terr<-read_sf("0.Data/WWF_ecoregions_realms/wwf_terr_ecos.shp") %>% st_transform(., crs(distributions))
realm_terr<-realm_terr %>% dplyr::group_by(REALM) %>% dplyr::summarise(N= n())
realm_terr$REALM_qtt<-revalue(realm_terr$REALM, c("AA"="1", "AN"="2", "AT"="3", "IM"="4", "NA"="5", "NT"="6", "OC"="7", "PA"="8")) %>% as.numeric(.)
realm_rast<-fasterize(realm_terr, raster(extent(realm_terr), res=c(10000,10000), crs = crs(realm_terr)), field="REALM_qtt", fun="first")
writeRaster(realm_rast, "0.Data/WWF_ecoregions_realms/Terrestrial_realms_rasterised.tif")
realm_rast<-raster("0.Data/WWF_ecoregions_realms/Terrestrial_realms_rasterised.tif") %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"), method='ngb')

# Extract
distributionsTERR[, c("realmsTERR_mode", "realmsTERR_var")]<-exact_extract(realm_rast, distributionsTERR, c("mode", "variety"))
species$realmsTERR_mode<-distributionsTERR$realmsTERR_mode[match(species$scientific_name, distributionsTERR$binomial)]
species$realmsTERR_var<-distributionsTERR$realmsTERR_var[match(species$scientific_name, distributionsTERR$binomial)]

# For the few species with NA (eg a small island not mapped in REALMS) get the closest 
for(i in which(is.na(species$realmsTERR_mode) & species$Dist_exist=="TRUE" & (species$terrestrial_system==TRUE | species$freshwater_system==TRUE))){
  species$realmsTERR_var[i]=1
  
  dist<-st_distance(distributionsTERR[as.character(distributionsTERR$binomial)==as.character(species$scientific_name[i]),], realm_terr)
  species$realmsTERR_mode[i]=which(dist==min(dist))
}

# Revalue with text names
species$realmsTERR_mode<-revalue(as.character(species$realmsTERR_mode), c("1"="AA", "2"="AN", "3"="AT", "4"="IM", "5"="NAm", "6"="NT", "7"="OC", "8"="PA", "9"="PA")) # 9 is returned for 9 fish species that occur in Caspian Sea or Baikal lake, they are Palearctic species






### Terrestrial accessibility
cat(paste("\n", "\n", "Starting Accessibility", Sys.time(), "\n"))

rm(realm_rast, sampling_effort, grid_group, Fetch_group, intersection_all) ; gc()
if("acc_terr" %not in% ls()){ # Charge only if not already charged
acc_terr<-raster("0.Data/Accessibility.to.cities.Weiss2018/accessibility_to_cities_2015_v1.0.tif") 
acc_terr<-replace(acc_terr, acc_terr==-9999, NA)
}

distributionsTERR[, c("acc_terrMED", "acc_terrQ10", "acc_terrQ90")]<-exact_extract(acc_terr, distributionsTERR, "quantile", quantiles=c(0.5, 0.1, 0.9))

species$acc_terrMED<-distributionsTERR$acc_terrMED[match(species$scientific_name, distributionsTERR$binomial)]
species$acc_terrQ10<-distributionsTERR$acc_terrQ10[match(species$scientific_name, distributionsTERR$binomial)]
species$acc_terrQ90<-distributionsTERR$acc_terrQ90[match(species$scientific_name, distributionsTERR$binomial)]




### Road density
roads<-raster("0.Data/GRIP4_density_total/grip4_total_dens_m_km2.asc")
crs(roads)<-"+init=epsg:4326"

distributionsTERR[, c("roadsMED", "roadsQ90")]<-exact_extract(roads, distributionsTERR, "quantile", quantiles=c(0.5, 0.9))

species$roadsMED<-distributionsTERR$roadsMED[match(species$scientific_name, distributionsTERR$binomial)]
species$roadsQ90<-distributionsTERR$roadsQ90[match(species$scientific_name, distributionsTERR$binomial)]





### Human population density
cat(paste("\n", "\n", "Starting Human population", Sys.time(), "\n"))

hpop<-raster("0.Data/GHS_population/GHS_POP_E2015_GLOBE_R2019A_54009_1K_V1_0.tif") # %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
distributionsTERR$hpop<-exact_extract(hpop, distributionsTERR, "mean")
species$hpop<-distributionsTERR$hpop[match(species$scientific_name, distributionsTERR$binomial)]

hrur<-raster("0.Data/GHS_population/Calculated_human_rural_population.tif")# %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
distributionsTERR$hrur<-exact_extract(hrur, distributionsTERR, "mean")
species$hrur<-distributionsTERR$hrur[match(species$scientific_name, distributionsTERR$binomial)]


write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1c(TEMPO3", Group, ").csv"), row.names=FALSE)




### Altitude
cat(paste("\n", "\n", "Starting Altitude", Sys.time(), "\n"))

alt<-raster("0.Data/GLOBEelevation/GLOBE_projected_Moll.tif") # %>% projectRaster(., crs=crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
distributionsTERR[, c("alt_med", "alt_Q10", "alt_Q90")]<-exact_extract(alt, distributionsTERR, "quantile", quantiles=c(0.5, 0.1, 0.9))

species$alt_med<-distributionsTERR$alt_med[match(species$scientific_name, distributionsTERR$binomial)]
species$alt_Q10<-distributionsTERR$alt_Q10[match(species$scientific_name, distributionsTERR$binomial)]
species$alt_Q90<-distributionsTERR$alt_Q90[match(species$scientific_name, distributionsTERR$binomial)]

write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1c(TEMPO4", Group, ").csv"), row.names=FALSE)





### Save
write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1c(TEMPO)", Group, ".csv"), row.names=FALSE)
