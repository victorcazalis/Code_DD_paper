

### Charge all distribution shapefiles
dist1<-st_read("0.Data/Distributions/FISH/ANGELFISH/ANGELFISH.shp")
dist2<-st_read("0.Data/Distributions/FISH/BLENNIES/BLENNIES.shp")
dist3<-st_read("0.Data/Distributions/FISH/BONEFISH_TARPONS/BONEFISH_TARPONS.shp")
dist4<-st_read("0.Data/Distributions/FISH/BUTTERFLYFISH/BUTTERFLYFISH.shp")
dist5<-st_read("0.Data/Distributions/FISH/CLUPEIFORMES/CLUPEIFORMES.shp")
dist6<-st_read("0.Data/Distributions/FISH/GROUPERS/GROUPERS.shp")
dist7<-st_read("0.Data/Distributions/FISH/HAGFISH/HAGFISH.shp")
dist8<-st_read("0.Data/Distributions/FISH/PUFFERFISH/PUFFERFISH.shp")
dist9<-st_read("0.Data/Distributions/FISH/SEABREAMS_PORGIES_PICARELS/SEABREAMS_PORGIES_PICARELS.shp")
dist10<-st_read("0.Data/Distributions/FISH/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp")
dist11<-st_read("0.Data/Distributions/FISH/SURGEONFISHES_TANGS_UNICORNFISHES/SURGEONFISHES_TANGS_UNICORNFISHES.shp")
dist12<-st_read("0.Data/Distributions/FISH/SYNGNATHIFORM_FISHES/SYNGNATHIFORM_FISHES.shp")
dist13<-st_read("0.Data/Distributions/FISH/TUNAS_BILLFISHES/TUNAS_BILLFISHES.shp")
dist14<-st_read("0.Data/Distributions/FISH/WRASSES_PARROTFISHES/WRASSES_PARROTFISHES.shp")

distMAR1<-st_read("0.Data/Distributions/FISH/MARINEFISH/MARINEFISH_PART1.shp")
distMAR2<-st_read("0.Data/Distributions/FISH/MARINEFISH/MARINEFISH_PART2.shp")
distMAR3<-st_read("0.Data/Distributions/FISH/MARINEFISH/MARINEFISH_PART3.shp")
distFW1<-st_read("0.Data/Distributions/FISH/FW_FISH/FW_FISH_PART1.shp")
distFW2<-st_read("0.Data/Distributions/FISH/FW_FISH/FW_FISH_PART2.shp")



### Compare species in groups and in Marine / Freshwaters
group_names<-c(dist1$binomial, dist2$binomial, dist3$binomial, dist4$binomial, dist5$binomial, dist6$binomial, dist7$binomial, dist8$binomial, dist9$binomial, dist10$binomial, dist11$binomial, dist12$binomial, dist13$binomial, dist14$binomial)
mar_names<-c(distMAR1$binomial, distMAR2$binomial, distMAR3$binomial)
fw_names<-c(distFW1$binomial, distFW2$binomial)


table(group_names %in% mar_names)
table(group_names %in% fw_names)
table(mar_names %in% group_names)
table(fw_names %in% group_names)

table(fw_names %in% mar_names)

table(group_names %in% c(fw_names, mar_names))
table(c(fw_names, mar_names) %in% group_names)


### Transform to character before merging (otherwise transforms to NA)
distMAR1<-distMAR1 %>% mutate_if(is.factor, as.character)
distMAR2<-distMAR2 %>% mutate_if(is.factor, as.character)
distMAR3<-distMAR3 %>% mutate_if(is.factor, as.character)
distFW1<-distFW1 %>% mutate_if(is.factor, as.character)
distFW2<-distFW2 %>% mutate_if(is.factor, as.character)
dist1<-dist1 %>% mutate_if(is.factor, as.character)
dist2<-dist2 %>% mutate_if(is.factor, as.character)
dist3<-dist3 %>% mutate_if(is.factor, as.character)
dist4<-dist4 %>% mutate_if(is.factor, as.character)
dist5<-dist5 %>% mutate_if(is.factor, as.character)
dist6<-dist6 %>% mutate_if(is.factor, as.character)
dist7<-dist7 %>% mutate_if(is.factor, as.character)
dist8<-dist8 %>% mutate_if(is.factor, as.character)
dist9<-dist9 %>% mutate_if(is.factor, as.character)
dist10<-dist10 %>% mutate_if(is.factor, as.character)
dist11<-dist11 %>% mutate_if(is.factor, as.character)
dist12<-dist12 %>% mutate_if(is.factor, as.character)
dist13<-dist13 %>% mutate_if(is.factor, as.character)
dist14<-dist14 %>% mutate_if(is.factor, as.character)




## Merge tables: first keep all marines, then add all freshwaters that are not in marine, then add all others
dist_merged<-distMAR1
dist_merged[(nrow(dist_merged)+1):(nrow(dist_merged)+nrow(distMAR2)),]<-distMAR2
dist_merged[(nrow(dist_merged)+1):(nrow(dist_merged)+nrow(distMAR3)),]<-distMAR3

FW1_toadd<-subset(distFW1, distFW1$binomial %not in% dist_merged$binomial) ; table(FW1_toadd$binomial %in% dist_merged$binomial) # Should be 100% FALSE
FW2_toadd<-subset(distFW2, distFW2$binomial %not in% dist_merged$binomial) ; table(FW2_toadd$binomial %in% dist_merged$binomial) # Should be 100% FALSE
dist_merged[(nrow(dist_merged)+1):(nrow(dist_merged)+nrow(FW1_toadd)),]<-FW1_toadd
dist_merged[(nrow(dist_merged)+1):(nrow(dist_merged)+nrow(FW2_toadd)),]<-FW2_toadd

table(dist10$binomial %in% dist_merged$binomial) # Check all tables, currently all species not in Marine / Freshwaters are in 10 (Chondrichthyes)

Chondr_toadd<-subset(dist10, dist10$binomial %not in% dist_merged$binomial) ; table(Chondr_toadd$binomial %in% dist_merged$binomial) # Should be 100% FALSE
dist_merged[(nrow(dist_merged)+1):(nrow(dist_merged)+nrow(Chondr_toadd)),]<-Chondr_toadd


### Final checks
table(group_names %in% dist_merged$binomial) # Should be 100% TRUE
table(mar_names %in% dist_merged$binomial) # Should be 100% TRUE
table(fw_names %in% dist_merged$binomial) # Should be 100% TRUE

### Save distributions (did not work in a single one, I save in two and merge later)
st_write(dist_merged[1:6034,], "0.Data/Distributions/FISH/FISHa.shp")
st_write(dist_merged[6035:20063,], "0.Data/Distributions/FISH/FISHb.shp")






