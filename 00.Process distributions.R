
for(Group in c("FW_ODONATA", "AMPHIBIANS", "MAMMALS", "REPTILES", "FISH")){


### Charge distributions
  if(Group != "FISH"){ # Charge distributions (in case of fishes, I need to combine two distributions)
    distributions<-st_read(paste0("0.Data/Distributions/", Group, "/", Group, ".shp")) 
  } else{
    distributions<-st_read("0.Data/Distributions/FISH/FISHa.shp") %>% mutate_if(is.factor, as.character)
    dist2<-st_read("0.Data/Distributions/FISH/FISHb.shp") %>% mutate_if(is.factor, as.character)
    
    distributions[(nrow(distributions)+1):(nrow(distributions)+nrow(dist2)),]<-dist2
  }


### Subset distributions
distributions<-subset(distributions, 
                      distributions$binomial %in% species$scientific_name &
                        distributions$presence %not in% c(3,4,5) & # Remove extinct and possibly extinct
                        distributions$seasonal %in% c(1,2) &  # Remove non-breeding grounds (keep only resident and breeding)
                        distributions$origin %in% c(1,2,5,6)) # Keep native, reintroduced, uncertain, assisted colonisation (remove introduced,  vagrant)


### Combine polygons
sf::sf_use_s2(FALSE) # Option to avoid bugs
distributions<-distributions %>% dplyr::group_by(binomial) %>% dplyr::summarise(N= n())





### Complete the maps of Odonata
if(Group=="FW_ODONATA"){

distributions<- distributions %>% st_transform(., st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

## Charge points data from the RL
pts_raw<-read.csv("0.Data/Distributions/FW_ODONATA/FW_ODONATA_points.csv")
tab_species<-as.data.frame(table(pts_raw$binomial))

pts<-subset(pts_raw, pts_raw$binomial %in% distributions$binomial | pts_raw$binomial %in% tab_species$Var1[tab_species$Freq>100]) # Keep only species with polygons, or species with at least 100 point data
ptsSHP<-SpatialPoints(cbind(pts$longitude, pts$latitude), proj4string=CRS("+init=epsg:4326")) %>% st_as_sf(.) %>% st_transform(., st_crs(distributions))
ptsSHP$binomial<-pts$binomial


## Charge hydrobasins
hydro<-st_read("0.Data/HydroBASINS_level08_w_attributes_2021_05/HydroBASINS_level08_w_attributes_2021_05.shp") %>% st_transform(., st_crs(distributions))





## Species loop
distributions$Check_Odonata<-NA
distributions$binomial<-as.character(distributions$binomial)
distributions$Check_Odonata<-as.character(distributions$Check_Odonata)

for (sp in 1:nlevels(droplevels(as.factor(ptsSHP$binomial)))){
  SP=levels(as.factor(ptsSHP$binomial))[sp]
  
  # Subset points and polygons
  ptsSP<-subset(ptsSHP, ptsSHP$binomial==SP)
  polySP<-subset(distributions, distributions$binomial==SP) 
  
  # Isolate points outside range
  inter<-st_join(ptsSP, polySP, join=st_intersects)
  ptsOUT<-subset(inter, is.na(inter$binomial.y)==T)
  
  # If no points outside range: we keep the current distribution
  if(nrow(ptsOUT)==0){
    distributions$Check_Odonata[distributions$binomial==SP]<-"AllPointsInRange"
  }
  
  # If points outside range + polygon
  if(nrow(ptsOUT)>0 & nrow(polySP)>0){
    interHyd<-st_join(ptsOUT, hydro, join=st_intersects) %>% subset(., is.na(.$hybas_id)==F) # Identify hydrobasins with data outside of range
    hydroSP<-subset(hydro, hydro$hybas_id %in% interHyd$hybas_id) # Isolate these hydrobasins
    
    if(nrow(hydroSP)==0){distributions$Check_Odonata[distributions$binomial==SP]<-"PointsOutInWater"} else{
    GHydro_added<-ggplot()+
           geom_sf(data=hydroSP, fill="gray95", colour="black", size=0.5)+
           geom_sf(data=polySP, col=NA, fill="#80cdc1")+
           geom_sf(data=ptsSP, col="forestgreen")+
           geom_sf(data=ptsOUT, col="#a6611a")+
           ggtitle(SP)+
           theme_void()%+replace%   theme(plot.background=element_rect(fill="white"))
    cowplot::save_plot(paste0("3.Figures/Check_Odonata_distribution/Hydro_added_to_poly/", SP, ".png"), GHydro_added, base_width=8, base_height=8) # Plot the distribution + points + hydrobasins
    
    distSP<-st_union(st_make_valid(polySP), hydroSP[, "geometry"]) %>% dplyr::group_by(binomial) %>% dplyr::summarise(N= n()) # Create a distribution that is the combined union of hydrobasins and polygon
    distSP$Check_Odonata<-"Hydro_added_to_poly"
    distributions[distributions$binomial==SP,]<-distSP # Change the distribution
    }}
  
  
  # If only points (in the end I did not keep these because the distributions were too incomplete)
  
  if(nrow(ptsOUT)>0 & nrow(polySP)==0){
    interHyd<-st_join(ptsOUT, hydro, join=st_intersects) %>% subset(., is.na(.$hybas_id)==F) # Identify hydrobasins with data outside of range
    hydroSP<-subset(hydro, hydro$hybas_id %in% interHyd$hybas_id) # Isolate these hydrobasins
    hydroSP$binomial<-SP
    
    if(nrow(hydroSP)==0){distributions$Check_Odonata[distributions$binomial==SP]<-"AllPointsInWater"} else{
    GHydro_tocreate<-ggplot()+
      geom_sf(data=hydroSP, fill="gray95", colour="black", size=0.5)+
      geom_sf(data=ptsOUT, col="#a6611a")+
      ggtitle(SP)+
      theme_void()%+replace%   theme(plot.background=element_rect(fill="white"))
    cowplot::save_plot(paste0("3.Figures/Check_Odonata_distribution/Hydro_to_create/", SP, ".png"), GHydro_tocreate, base_width=8, base_height=8) # Plot the distribution + points + hydrobasins
    
    distSP<-hydroSP[, c("binomial", "geometry")] %>% dplyr::group_by(binomial) %>% dplyr::summarise(N= n()) # Create a distribution that is the combined union of hydrobasins and polygon
    distSP$Check_Odonata<-"Hydro_to_create"
    distributions[(nrow(distributions)+1),]<-distSP # Change the distribution
  }}

  
  
  # Final check (for species with no polygon or with points outside polygons)
  if(nrow(ptsOUT)>0 | nrow(polySP)==0){
  GDistri_final<-ggplot()+
    geom_sf(data=distributions[distributions$binomial==SP,], col=NA, fill="#80cdc1")+
    geom_sf(data=ptsSP, col="forestgreen")+
    ggtitle(SP)+
    theme_void()%+replace%   theme(plot.background=element_rect(fill="white"))
  cowplot::save_plot(paste0("3.Figures/Check_Odonata_distribution/AllPointsInRange/", SP, ".png"), GDistri_final, base_width=8, base_height=8)
  }
  
  rm(polySP, ptsSP, ptsOUT, interHyd, inter, GDistri, distSP, GHydro_added, hydroSP, GDistri_final)
  print_loop(sp, 10)
}

} # End if Odonata








### Save distributions (in two files for fishes)
  if(Group=="FISH"){N_lim<-6999
  st_write(distributions[1:N_lim,], paste0("0.Data/Distributions/processed/Dist.processed.", Group, "a.shp"))
  st_write(distributions[(N_lim+1):nrow(distributions),], paste0("0.Data/Distributions/processed/Dist.processed.", Group, "b.shp"))
} else{
st_write(distributions, paste0("0.Data/Distributions/processed/Dist.processed.", Group, ".shp"))
}
}



