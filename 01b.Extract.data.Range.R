

### Charge distributions
if(Group =="FISH"){
  distributions<-st_read(paste0("0.Data/Distributions/processed/Dist.processed.", Group, "a.shp")) ; distributions$binomial<-as.character(distributions$binomial)
  distb<-st_read(paste0("0.Data/Distributions/processed/Dist.processed.", Group, "b.shp")) ;  distb$binomial<-as.character(distb$binomial)
  distributions[(nrow(distributions)+1):(nrow(distributions)+nrow(distb)),]<-distb ; distb$binomial<-as.factor(distb$binomial) ; rm(distb)
} else{
  distributions<-st_read(paste0("0.Data/Distributions/processed/Dist.processed.", Group, ".shp"))
}

names(distributions)<-revalue(names(distributions), c("binomil"="binomial", "Chck_Od"="Check_Odonata"))

sf::sf_use_s2(FALSE)
distributions<-st_transform(distributions, st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))


# Increase a bit the range of seasnakes to include the close coast
if(Group=="REPTILES"){
  distributions$tobuffer<-distributions$binomial %in% species$scientific_name[species$marine_system==T & species$terrestrial_system==F & species$freshwater_system==F]
  distributions[distributions$tobuffer==T,]<-st_buffer(distributions[distributions$tobuffer==T,], 5000)
}

# Create a variable saying if we have a distribution (and distribution quality for Odonata, i.e. did I expanded the distribution or not)
species$Dist_exist<-species$scientific_name %in% distributions$binomial ; table(species$Dist_exist)

if(Group=="FW_ODONATA"){species$Dist_quality<-distributions$Check_Odonata[match(species$scientific_name, distributions$binomial)] ; table(species$Dist_quality)}






### Range size
cat(paste("\n", "\n", "Starting Range size", Sys.time(), "\n"))

# Calculate range size
species$range_km2<-NA
distributions$binomial<-as.character(distributions$binomial) ; species$scientific_name<-as.character(species$scientific_name)

for(i in 1:nrow(species)){
  if(species$Dist_exist[i]==T){
    species$range_km2[i]<-as.numeric(st_area(distributions[distributions$binomial==species$scientific_name[i],]))/1000000
    if((i/100)==round(i/100)){cat(i)}
  }}






### GBIF sampling effort
cat(paste("\n", "\n", "Starting GBIF Sampling effort", Sys.time(), "\n"))


## Create a map of sampling effort for the broad taxonomic group

# Choose settings
CRS_chosen<-"EPSG:4326" 
FORMAT_chosen<-"@4x.png"

# Create empty grid (exactly the one used by mvt_fetch)
XMIN= -25046885 - 400 # The first coordinate is the minimum point of Fetch. I move a bit to align the points of Fetch with top left corner of the pixels, that's how I think it works (looking at specific localities of data on GBIF)
YMIN= -19959237
C_size=78272
N_x=641
N_y=511
poly<-st_sfc(st_polygon(list(rbind(c(XMIN,YMIN), c(XMIN+N_x*C_size,YMIN), c(XMIN+N_x*C_size, YMIN+N_y*C_size), c(XMIN,YMIN)))))
empty_grid<-st_make_grid(poly, n=c(N_x, N_y), what="polygons")
st_crs(empty_grid)<-st_crs("EPSG:3857") ; rm(XMIN, YMIN, poly )# Has to be Mercator 3857 to create a regular empty grid as in mvt_fetch, but Mercator creates bugs in some Pacific species (e.g., Acanthaeschna victoria) so I use PlateCarre for mvt_fetch
# To check that the grid is exactly the same: st_write(empty_grid, "Grid.shp", append=F) ; st_write(mvt_fetch(taxonKey = 0, srs = CRS_chosen, format=FORMAT_chosen), "Mvt_fetch_grid.shp", append=F)

r<-raster(extent(st_as_sf(empty_grid)), res=c(C_size,C_size), crs = crs(empty_grid))    
r[] <- 1



# Create the map of sampling effort
if(Group != "FISH"){
  
  if(Group=="MAMMALS"){keyGroup<-name_suggest(q = "Mammalia", rank = "Class")$data$key}
  if(Group=="REPTILES"){keyGroup<-name_suggest(q = "Reptilia", rank = "Class")$data$key}
  if(Group=="AMPHIBIANS"){keyGroup<-name_suggest(q = "Amphibia", rank = "Class")$data$key}
  if(Group=="FW_ODONATA"){keyGroup<-name_suggest(q = "Odonata", rank = "Order")$data$key}
  
  Fetch_group<-mvt_fetch(taxonKey = keyGroup, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  
  grid_group<-empty_grid %>% st_as_sf(.)
  inter_group<-st_join(grid_group, Fetch_group, join=st_intersects)
  
  sampling_effort<-fasterize(inter_group, r, field="total", fun="sum")
  
}

if(Group == "FISH"){ # FISHES: There are 5 classes in the RL: ACTINOPTERYGII, CEPHALASPIDOMORPHI, CHONDRICHTHYES (in GBIF split between Elasmobranchii and Holocephali), MYXINI, SARCOPTERYGII
  keyFish1<-name_suggest(q = "Actinopterygii", rank = "Class")$data$key
  keyFish2<-name_suggest(q = "Cephalaspidomorphi", rank = "Class")$data$key
  keyFish3a<-name_suggest(q = "Elasmobranchii", rank = "Class")$data$key
  keyFish3b<-name_suggest(q = "Holocephali", rank = "Class")$data$key
  keyFish4<-name_suggest(q = "Myxini", rank = "Class")$data$key
  keyFish5<-name_suggest(q = "Sarcopterygii", rank = "Class")$data$key
  
  Fetch_fish1<-mvt_fetch(taxonKey = keyFish1, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  Fetch_fish2<-mvt_fetch(taxonKey = keyFish2, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  Fetch_fish3a<-mvt_fetch(taxonKey = keyFish3a, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  Fetch_fish3b<-mvt_fetch(taxonKey = keyFish3b, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  Fetch_fish4<-mvt_fetch(taxonKey = keyFish4, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  Fetch_fish5<-mvt_fetch(taxonKey = keyFish5, srs = CRS_chosen, format=FORMAT_chosen) %>% st_transform(., st_crs(empty_grid)) 
  
  grid_group<-empty_grid %>% st_as_sf(.)
  
  inter_fish1<-st_join(grid_group, Fetch_fish1, join=st_intersects)
  inter_fish2<-st_join(grid_group, Fetch_fish2, join=st_intersects)
  inter_fish3a<-st_join(grid_group, Fetch_fish3a, join=st_intersects)
  inter_fish3b<-st_join(grid_group, Fetch_fish3b, join=st_intersects)
  inter_fish4<-st_join(grid_group, Fetch_fish4, join=st_intersects)
  inter_fish5<-st_join(grid_group, Fetch_fish5, join=st_intersects)
  
  sampling_fish1<-fasterize(inter_fish1, r, field="total", fun="sum")
  sampling_fish2<-fasterize(inter_fish2, r, field="total", fun="sum")
  sampling_fish3a<-fasterize(inter_fish3a, r, field="total", fun="sum")
  sampling_fish3b<-fasterize(inter_fish3b, r, field="total", fun="sum")
  sampling_fish4<-fasterize(inter_fish4, r, field="total", fun="sum")
  sampling_fish5<-fasterize(inter_fish5, r, field="total", fun="sum")
  
  st <- stack(sampling_fish1, sampling_fish2, sampling_fish3a, sampling_fish3b, sampling_fish4, sampling_fish5)
  sampling_effort <- sum(st, na.rm=T)
}



# Extract sampling effort for each species
if(nlevels(droplevels(as.factor(species$Vgroup)))==1){
  
  distributions$Sampling_MEAN<-exact_extract(sampling_effort, distributions, "mean")
  distributions[, c("Sampling_MED", "Sampling_Q10", "Sampling_Q90")]<-exact_extract(sampling_effort, distributions, "quantile", quantiles=c(0.5, 0.1, 0.9))  
  
  species$Sampling_MEAN<-distributions$Sampling_MEAN[match(species$scientific_name, distributions$binomial)]
  species$Sampling_MED<-distributions$Sampling_MED[match(species$scientific_name, distributions$binomial)]
  species$Sampling_Q10<-distributions$Sampling_Q10[match(species$scientific_name, distributions$binomial)]
  species$Sampling_Q90<-distributions$Sampling_Q90[match(species$scientific_name, distributions$binomial)]
  
}


write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1b.(TEMPO1)", Group, ".csv"), row.names=FALSE)







### GBIF coverage of range
cat(paste("\n", "\n", "Starting GBIF coverage of range", Sys.time(), "\n"))

species$gbifSP_check<-species$gbifSP_propempty<-species$gbifSP_propemptyASS<-NA



for(i in 1:nrow(species)){ 
  
  if(species$nb_GBIFgeo[i]==0){species$gbifSP_check[i]<-"No Observation"} else {
  
  keySP<-name_backbone(name=species$scientific_name[i])$usageKey
  range_SP<-subset(distributions, distributions$binomial==species$scientific_name[i]) %>% st_transform(., st_crs(empty_grid))
  Y_ass<-as.numeric(format(as.Date(species$assessment_date[i], "%Y-%m-%d"), "%Y"))

  if(nrow(range_SP)==0){species$gbifSP_check[i]<-"No Range"}else{
    
  Fetch<-mvt_fetch(taxonKey = keySP, srs = CRS_chosen, format=FORMAT_chosen) 
  if(is.null(Fetch)){species$gbifSP_check[i]<-"Fetch Null"}else{
  
    Fetch<-Fetch %>% st_transform(., st_crs(range_SP)) 

    # Calculate coverage 
    grid_sp<-st_crop(empty_grid, range_SP) %>% st_as_sf(.)
    inter<-st_join(grid_sp, range_SP, join=st_intersects) %>% subset(., is.na(N)==F)
    inter<-st_join(inter, Fetch, join=st_intersects)
    #Plot to check: ggplot()+ geom_sf(data=inter, aes(fill=is.na(total)))+ geom_sf(data=range_SP, col="darkred", fill=NA)+ geom_sf(data=Fetch, aes(col=total))


    # Extract values
    species$gbifSP_propempty[i]<-as.numeric(prop.table(table(factor(is.na(inter$total), c("TRUE", "FALSE"))))["TRUE"])
    interASS<-as.data.frame(inter[,which(names(inter)<Y_ass)]) ; interASS$geometry<-NULL
    interASS$totalASS<-rowSums(interASS, na.rm=T)
    species$gbifSP_propemptyASS[i]<-as.numeric(prop.table(table(factor(interASS$totalASS>0, c("TRUE", "FALSE"))))["FALSE"])

    # Proportion of data in cells outside the range
    species$gbif_outside_range[i]<- 1 - (sum(inter$total, na.rm=T) / sum(Fetch$total, na.rm=T))
    FetchASS<-as.data.frame(Fetch[,which(names(Fetch)<Y_ass)]) ; FetchASS$geometry<-NULL
    FetchASS$totalASS<-rowSums(FetchASS, na.rm=T)
    species$gbif_outside_rangeASS[i]<- 1 - (sum(interASS$totalASS, na.rm=T) / sum(FetchASS$totalASS, na.rm=T))

  }}} 
  print_loop(i, 100)
}





### How many DD species in range in the broad taxonomic group
cat(paste("\n", "\n", "Starting DDness", Sys.time(), "\n"))

# Create a dataframe with each species in a grid of 100km
grid<-st_make_grid(distributions, cellsize=100000, square=F) %>% st_transform(., st_crs(distributions))

intersection_all<-c()
for (i in 1:nrow(distributions)){
  print_loop(i, 100)
  
  intersection_sps_in_loop <- st_intersects(distributions[i,], grid) %>% as.data.frame(.)
  intersection_sps_in_loop$row.id=distributions$binomial[i]
  intersection_all<-rbind(intersection_all,intersection_sps_in_loop)
}  

saveRDS(intersection_all, file=paste0("1.Tables/intersection_all_gridDDness_TEMPO.", Group, ".rds"))


# Calculate for each species
intersection_all$row.id<-droplevels(as.factor(intersection_all$row.id))
list_DDs<-droplevels(as.factor(species$scientific_name[species$category=="DD"]))

for(i in 1:nlevels(intersection_all$row.id)){
  SP=levels(intersection_all$row.id)[i]
  Grid<-subset(intersection_all, intersection_all$row.id != SP & intersection_all$col.id %in% intersection_all$col.id[intersection_all$row.id==SP]) 
  spcs<-unique(droplevels(Grid$row.id))
  
  species$N_tot[species$scientific_name==SP]<-length(spcs)
  species$N_DDs[species$scientific_name==SP]<-length(spcs[spcs %in% list_DDs])
  species$DDness[species$scientific_name==SP]<- nrow(Grid[Grid$row.id %in% list_DDs,])/nrow(Grid)
  print_loop(i, 100)
}

species$propDD<-species$N_DDs/species$N_tot





### Save
write.csv(species, paste0("1.Tables/Species.characteristics.datasufficiency.Script1b.", Group, ".csv"), row.names=FALSE)

