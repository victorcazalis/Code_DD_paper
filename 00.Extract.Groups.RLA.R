

species$RLA<-species$Dist_group<-species$Vgroup<-NA



# ### Red List Authority
RLAs<-read.csv("0.Data/SSC_Groups/Published_scc_groups_v1.csv")
species_raw$RLA<-RLAs$name[match(species_raw$scientific_name, RLAs$scientific_name)]


### Large group
species$Vgroup[species$order=="ODONATA"]<-"Odonata"
species$Vgroup[species$class %in% c("ACTINOPTERYGII", "CEPHALASPIDOMORPHI", "CHONDRICHTHYES", "MYXINI", "SARCOPTERYGII")]<-"Fish"
species$Vgroup[species$class=="AMPHIBIA"]<-"Amphibian"
species$Vgroup[species$class=="REPTILIA"]<-"Reptile"
species$Vgroup[species$class=="AVES"]<-"Bird"
species$Vgroup[species$class=="MAMMALIA"]<-"Mammal"

table(is.na(species$Vgroup))


### Distribution name
species$Dist_group<-revalue(species$Vgroup, c("Odonata"="FW_ODONATA", "Fish"="FISH", "Amphibian"="AMPHIBIANS", "Reptile"="REPTILES", "Bird"="BIRDS", "Mammal"="MAMMALS"))

