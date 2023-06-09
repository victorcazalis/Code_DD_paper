library(sf) ; library(exactextractr) ; library(gfcanalysis)

CRS_forest<-"+init=epsg:4326"


### Keep only forest specialists (i.e., the only terrestrial habitat is forest)
speciesDDFor<-subset(speciesDD, Bin_forest==T & 
                       Bin_savannah==F & Bin_shrubland==F & Bin_grassland==F & Bin_rocky==F & Bin_cave==F & Bin_desert==F & Bin_artiTerr==F & Bin_introVeg==F)


distDDfor<-subset(distDD, distDD$binomial %in% speciesDDFor$scientific_name) %>%
  st_transform(., st_crs(CRS_forest))
rm(distDD)

sf::sf_use_s2(FALSE)



### List needed tiles
tilesSP<-calc_gfc_tiles(distDDfor)
tilesSP$Code<-paste0("Code", 1:length(tilesSP))
tiles<-st_as_sf(tilesSP)
intersect<-st_intersection(tiles, distDDfor)
tiles<-subset(tiles, tiles$Code %in% intersect$Code)
tilesSP<-subset(tilesSP, tilesSP$Code %in% intersect$Code)

### Format
tile_crs<-raster(paste0(Main_dir, "Trends EoN/Analyses/Global maps/Data/Hansen_tiles/Hansen_GFC-2019-v1.7_treecover2000_60N_000E.tif"))
distDDfor<-st_transform(distDDfor, crs(tile_crs))
intersect<-st_transform(intersect, crs(tile_crs))
intersect$Code<-as.factor(intersect$Code)
intersect$loss3GL<-intersect$loss<-intersect$area<-intersect$gain<-intersect$Forest2000<-NA

# Add GL
intersect$GL<-GL$GL_estimate[match(intersect$binomial, GL$internal_taxon_name)] ; intersect$GL[is.na(intersect$GL)==T]<-1
for(i in 1:nrow(intersect)){intersect$Year_GL[i] <- (2021-max(round(3*intersect$GL[i]), 10)) %>% max(., 2000)}



### EXTRACT VALUES
for(TILE in 1:nrow(tiles)){
  # Clean Temp
  unlink(tempdir(), recursive=T)
  gc()
  
  
  EXIST=0 # Always put back to 0, tile will be downloaded only it EXIST stays at 0
  
  interSUB<-subset(intersect, intersect$Code==levels(intersect$Code)[TILE])
  cat("Tile n", TILE, " (", length(which(is.na(interSUB$Forest2000))), ")", "\n") ; rm(tile)
  if(length(which(is.na(interSUB$Forest2000))) >0){
  
  ### Charge the tile
  # If it's already downloaded in Hansen_tiles, charge it
  tryCatch({tile<-extract_gfc(tiles[tiles$Code==levels(intersect$Code)[TILE],], data_folder="0.Data/Hansen_tiles", dataset="GFC-2021-v1.9") ; EXIST=1}, error=function(e){cat("Not in data sufficiency", "\n")})
  
  # If not, download
  if(EXIST==0){
    cat("Downloads")
    download_tiles(
      tilesSP[tiles$Code==levels(intersect$Code)[TILE],],
      output_folder="0.Data/Hansen_tiles",
      images = c("treecover2000", "lossyear", "gain", "datamask"),
      dataset = "GFC-2021-v1.9")
    tile<-extract_gfc(tiles[tiles$Code==levels(intersect$Code)[TILE],], data_folder="0.Data/Hansen_tiles", dataset="GFC-2021-v1.9")
  }
  
  
  
  for(i in 1:nrow(interSUB)){
  cat(i)
    tryCatch({
      tile_crop<-crop(tile, interSUB[i,])
      
      tile_thresh<-threshold_gfc(tile_crop, forest_threshold = 25) # Default
      
      intersect$Forest2000[intersect$Code==interSUB$Code[i] & intersect$binomial==interSUB$binomial[i]]<-exact_extract(tile_thresh$forest2000, interSUB[i,], fun="sum")
      intersect$gain[intersect$Code==interSUB$Code[i] & intersect$binomial==interSUB$binomial[i]]<-exact_extract(tile_thresh$gain, interSUB[i,], fun="sum")
      intersect$area[intersect$Code==interSUB$Code[i] & intersect$binomial==interSUB$binomial[i]]<-exact_extract(replace(tile_thresh$datamask, tile_thresh$datamask>1, 0), interSUB[i,], fun="sum")
      intersect$loss[intersect$Code==interSUB$Code[i] & intersect$binomial==interSUB$binomial[i]]<-exact_extract(replace(tile_thresh$lossyear, tile_thresh$lossyear>0, 1), interSUB[i,], fun="sum")
      intersect$loss3GL[intersect$Code==interSUB$Code[i] & intersect$binomial==interSUB$binomial[i]]<-exact_extract(tile_thresh$lossyear>=(interSUB$Year_GL[i]-2000+1), interSUB[i,], fun="sum") # Everything lost from the year after Year_GL
      
    } ,error=function(e){cat(paste("Bug at TILE=", TILE, " and i=", i, "\n"))})
  }
  
  cat("\n", TILE, "\n");
  if((TILE/10)==round(TILE/10)){
    DF<-data.frame(Code=intersect$Code, binomial=intersect$binomial, Forest2000=intersect$Forest2000, gain=intersect$gain, area=intersect$area, loss=intersect$loss, loss3GL=intersect$loss3GL)
    write.csv(DF, paste0("1.Tables/AOH.TEMPOFOR", Group, ".", TILE, ".csv"), row.names=F)
    }
}
}


# ### Fix the bugs
# for(i in which(is.na(intersect$Forest2000))){
#   
#   tryCatch({
#     
#     urban_sub<-st_buffer(intersect[i,],0)
#     tile<-extract_gfc(urban_sub, data_folder="H:/Trends EoN/Analyses/Global maps/Data/Hansen_tiles", dataset="GFC-2021-v1.9")
#     
#     tile_thresh<-threshold_gfc(tile, forest_threshold = 25) # Default
#     
#     intersect$Forest2000[i]<-exact_extract(tile_thresh$forest2000, urban_sub, fun="sum")
#     intersect$gain[i]<-exact_extract(tile_thresh$gain, urban_sub, fun="sum")
#     intersect$area[i]<-exact_extract(replace(tile_thresh$datamask, tile_thresh$datamask>1, 0), urban_sub, fun="sum")
#     intersect$loss[i]<-exact_extract(replace(tile_thresh$lossyear, tile_thresh$lossyear>1, 1), urban_sub, fun="sum")
#     
#     cat(paste0(i, " "))
#     
#   } ,error=function(e){cat(paste("Bug at i=", i, "\n"))})
#   
# }

# Save raw extraction
intersectDF<-as.data.frame(intersect)
intersectDF$geometry<-NULL
write.csv(intersectDF, paste0("1.tables/Hansen.extract.raw.", Group, ".csv"))

### SUM
interSUM<-ddply(intersectDF, .(binomial), function(x){data.frame(
  Forest2000=sum(x$Forest2000, na.rm=T),
  gain=sum(x$gain, na.rm=T),
  area=sum(x$area, na.rm=T),
  loss=sum(x$loss, na.rm=T),
  loss3GL=sum(x$loss3GL, na.rm=T)
)})

interSUM$Forest3GL<-interSUM$Forest2000 - (interSUM$loss - interSUM$loss3GL)
interSUM$Forest2021<-interSUM$Forest2000 - interSUM$loss
interSUM$ForestLoss<- (interSUM$loss3GL / interSUM$Forest3GL) # loss3GL corresponds exactly to the difference between FOrest3GL and Forest2021

speciesDD$ForestLoss<-interSUM$ForestLoss[match(speciesDD$scientific_name, interSUM$binomial)]
speciesDD$Forest2021<-interSUM$Forest2021[match(speciesDD$scientific_name, interSUM$binomial)]

### Save
write.csv(speciesDD, paste0("1.Tables/AOH.", Group, ".csv"), row.names=F)




