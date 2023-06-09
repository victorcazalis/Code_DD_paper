

### Extract the id of the most recent assessment
CorrAss<-readRDS("0.Data/SIS/CorrAss_adaptedfromSIS_PUBLICATIONDATE.rds")
species$assessmentid_V<-CorrAss$Assessment[match(species$taxonid, CorrAss$taxonid)]



### Number of GBIF records (in 2022 and before last assessment)
cat("\n", "Start GBIF extract")

species$nb_GBIFgeo<-species$nb_GBIF<-species_GBIFgeoASS<-species$nb_GBIFASS<-NA
Now<-as.numeric(format(Sys.Date(), "%Y"))

for(i in 1:nrow(species)){
  Y_sp<-as.numeric(format(as.Date(species$assessment_date[i], "%Y-%m-%d"), "%Y"))
  
  key<-name_backbone(name=species$scientific_name[i])$usageKey
  species$nb_GBIFgeo[i]<-occ_count(taxonKey=key, georeferenced = TRUE)
  
  N_geo<-0
  for(j in Y_sp:Now){
    N_geo<-N_geo+occ_count(taxonKey=key, year=j, georeferenced=TRUE)
  }
  species$nb_GBIFgeoASS[i]<-species$nb_GBIFgeo[i]-N_geo ; rm(N_geo)
  
  print_loop(i, 100)
}



### Number of data in OBIS (for fishes only)
species$nb_OBIS<-species$nb_OBISASS<-NA

for(i in which(species$Vgroup=="Fish")){
  dat_sp<-occurrence(scientificname=species$scientific_name[i])
  dat_sp<-subset(dat_sp, is.na(dat_sp$decimalLatitude)==F) # To keep only georeferenced data (to be able to compare with rgbif)
  species$nb_OBIS[i]<-nrow(dat_sp)
  
  Y_sp<-as.numeric(format(as.Date(species$assessment_date[i], "%Y-%m-%d"), "%Y"))
  if("date_year" %in% names(dat_sp)){
    dat_spASS<-subset(dat_sp, dat_sp$date_year < Y_sp | is.na(dat_sp$date_year)==T)
    species$nb_OBISASS[i]<-nrow(dat_spASS)} else {
      species$nb_OBISASS[i]<-nrow(dat_sp)
    }
  print_loop(i, 100)
}




### Number of articles in WOS
species$WOSASS<-species$WOS<-NA

Date1<-"1500-01-01" # First date for research (past limit)
Date3<-as.character(Sys.Date()) # Third date for research (current)


for(sp in 1:nrow(species)){
  Date2<-species$assessment_date[sp] # Second date for research (time of last assessment)
  
  QUERY_WOS <- paste0('TS=("', species$scientific_name[sp], '")') # Important to have real quote marks (" is ok ' is not)

  # Extract the number of articles between Date1 and Date3 (i.e., total)
  resp1 <- httr::GET('https://api.clarivate.com/api/woslite', 
                     httr::add_headers(accept = 'application/json', `X-APIKey` = KEY_WOS), # A token to the WOS lite API is required for this step (stored as the object KEY_WOS)
                     query = list(databaseId = 'WOK', 
                                  usrQuery = QUERY_WOS, 
                                  count = 100,  # Can be up to 100 for the list of papers
                                  firstRecord = 1, 
                                  publishTimeSpan=paste(Date1, Date3, sep="+")))
  
  data1 = fromJSON(rawToChar(resp1$content))
  species$WOS[sp]<-data1$QueryResult$RecordsFound
  
  # Extract the number of articles between Date1 and Date2 (i.e., number of articles available at last assessment)
  resp2 <- httr::GET('https://api.clarivate.com/api/woslite', 
                     httr::add_headers(accept = 'application/json', `X-APIKey` = KEY_WOS), 
                     query = list(databaseId = 'WOK', 
                                  usrQuery = QUERY_WOS, 
                                  count = 1, # Can be up to 100 for the list of papers
                                  firstRecord = 1,
                                  publishTimeSpan=paste(Date1, Date2, sep="+")))
  
  data2 = fromJSON(rawToChar(resp2$content))
  species$WOSASS[sp]<-data2$QueryResult$RecordsFound
  
  print_loop(sp, 50)
}





### Time since description 
cat("\n", "Start Time since description")

species$description_year<-NA

for(i in 1:nrow(species)){
  v<-unlist(strsplit(species$authority[i]," "))
  species$description_year[i]<-substr(v[length(v)], 1,4)
}

species$description_year<-as.numeric(as.character(species$description_year))

table(species$description_year[species$description_year %not in% 1500:2022])
species$description_year<-revalue(as.character(species$description_year), c('201'=2010)) %>% as.numeric(.) # Manual fix for a spelling mistake in the Red List data




### Habitat
cat("\n", "Start Habitat extract")

Habs_raw<-readRDS("0.Data/SIS/Habitats_adaptedfromSIS_GENERALHABITATSSUBFIELD.rds")


# Prepare depth for marine species
Habs_raw$depth_index<-replace(Habs_raw$Code0, Habs_raw$Code0 %in% c(1, 2, 3, 4, 6, 7, 8, 14, 15, 16, 17, 18), NA) %>% 
  revalue(., c("5"="1", "9"="1", "12"="1", "13"="1")) %>% 
  as.character(.)

Habs_raw$depth_index[Habs_raw$Code0 %in% c(10,11)]<-revalue(Habs_raw$Code[Habs_raw$Code0 %in% c(10,11)], c(
  "10.1"="1", 
  "10.2"="2", "11.1"="2", "11.1.1"="2", "11.1.2"="2", "10"="2",
  "10.3"="3", "11.5"="3", "11.6"="3", "11"="3",
  "10.4"="4", "11.2"="4", "11.3"="4", "11.4"="4")) 

Habs_raw$depth_index<-as.numeric(Habs_raw$depth_index)

# Extract per species
habitats<-ddply(Habs_raw, .(taxonid), function(x){data.frame(
  N.hab1=length(unique(x$Code0)),
  Major=subset(x, x$majorimportance_value=="Yes")$Code0 %>% unique(.) %>% paste(., collapse=" "),
  Bin_forest= 1 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_savannah= 2 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_shrubland= 3 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_grassland= 4 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_wetland= 5 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_rocky= 6 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_cave= 7 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_desert= 8 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_neretic= 9 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_oceanic= 10 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_deep= 11 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_intertidal= 12 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_coastal= 13 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_artiTerr= 14 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_artiAqua= 15 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  Bin_introVeg= 16 %in% subset(x, x$suitability_value=="Suitable")$Code0,
  
  Depth_min=min(x$depth_index, na.rm=T),
  Depth_max=max(x$depth_index, na.rm=T)
)})

species<-left_join(species, habitats, by="taxonid")
species$Depth_min[species$Depth_min==Inf]<-NA
species$Depth_max[species$Depth_max==-Inf]<-NA

# Fill in Depth with the median in the family when unknown (1 for Arctic fox)
depth_fam<-ddply(species[species$marine_system==T,], .(family), function(x){data.frame(Depth_min_fam=median(x$Depth_min, na.rm=T), Depth_max_fam=median(x$Depth_max, na.rm=T))})
species$Depth_min[is.na(species$Depth_min)==T & species$marine_system==T]<-depth_fam$Depth_min_fam[match(species$family[is.na(species$Depth_min)==T & species$marine_system==T], depth_fam$family)]
species$Depth_max[is.na(species$Depth_max)==T & species$marine_system==T]<-depth_fam$Depth_max_fam[match(species$family[is.na(species$Depth_max)==T & species$marine_system==T], depth_fam$family)]





### ADD COUNTRIES STATs (GDP, Corruption, Research capacity)
cat("\n", "Start Countries extract")

countriesTOT<-readRDS("0.Data/SIS/Countries_adaptedfromSIS_COUNTRYOCCURRENCESUBFIELD.rds")
table(species$taxonid %in% countriesTOT$taxonid) # The false are because no Country of Occurrence was mentioned (some fishes with a range map in the middle of the ocean, reptiles with no range)

### Prepare GDP data (take the youngest GDP for each country, see e.g. Syria)
gdp<-read.csv("0.Data/GDP/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_3159040.csv", sep="\t")

gdp$GDP<-NA
for(i in 1:nrow(gdp)){
  v<-gdp[i,grepl("X", names(gdp))]
  v<-v[is.na(v)==F]
  if(length(v)>0){gdp$GDP[i]<-v[length(v)]}
}

gdp<-gdp[, c("Country.Name", "Country.Code", "GDP")]


### Prepare Research capacity data
ResCap<-read.csv("0.Data/Research capacity/scimagojr country rank 1996-2020.csv", sep=",")


### Prepare Corruption data
Corruption<-read.csv("0.Data/Corruption/CPI2020_GlobalTablesTS_210125.csv", sep=",") ; names(Corruption)[1]<-"Country"


### Prepare conflicts data
conflicts<-openxlsx::read.xlsx("H:/Postdoc/Data sufficiency/1.Analyses/0.Data/Conflicts/ucdp-prio-acd-211.xlsx")
conflicts<-subset(conflicts, conflicts$year>=2000 & conflicts$cumulative_intensity==1)

# Assign manually location for some missing data
conflicts$location[conflicts$conflict_id=="13692"]<-"Afghanistan" # War in Afghanistan only located in Afghanistan
conflicts$location[conflicts$conflict_id=="420"]<-"Iraq" # War in Iraq only located in Iraq
conflicts$location[conflicts$conflict_id=="274"]<-"India" # Conflict in Ladaq
conflicts$location[conflicts$conflict_id=="218"]<-"India" # In reality both in India and Pakistan (but both are already in the list)
conflicts$location[conflicts$conflict_id=="409"]<-"Eritrea" # In reality both in Eritrea and Ethiopia (but Ethiopia is already in the list)

conflicts_countries<-ddply(conflicts, .(location), function(x){data.frame(Freq=nlevels(as.factor(x$year)))})




### Match 4 variables with IUCN countries 
ctr<-read.csv("0.Data/Countries_crosswalk.csv")

# GDP
ctr$GDP<-gdp$GDP[match(ctr$Countries_WB, gdp$Country.Name)]
ctr$GDP[ctr$Countries_WB=="Korea, Dem. People's Rep."]<-18000000000/25778810 # Only country with missing data (for Gibraltar and St Martin I took UK and France), estimate taken from this site: https://tradingeconomics.com/north-korea/gdp and population size here: https://data.worldbank.org/indicator/SP.POP.TOTL?view=map
table(is.na(ctr$GDP)) # Antarctica is NA

# Research capacity
ctr$res.cap<-ResCap$Documents[match(ctr$Countries_SCIMAJ, ResCap$Country)]
table(is.na(ctr$res.cap)) # Antarctica is NA

# Corruption
ctr$corruption<-Corruption$CPI.score.2020[match(ctr$Countries_Corruption, Corruption$Country)]
table(is.na(ctr$corruption)) # There are some NAs (eg independent Caribbean islands, so there are some missing values)

# Conflicts
ctr$conflict<-conflicts_countries$Freq[match(ctr$Countries_conflicts, conflicts_countries$location)]
ctr$conflict<-replace(ctr$conflict, is.na(ctr$conflict), 0)
table(conflicts_countries$location %in% ctr$Countries_conflicts) # Should be 100% TRUE


### Create empty variables in species data frame
species$gdpOK<-species$nb_countries<-species$med_GDP<-NA
species$rescapOK<-species$res.cap<-NA
species$corruptionOK<-species$corruption<-NA
species$conflictMEANYEAR<-species$conflictMEDYEAR<-species$conflictONLY<-species$conflictEXIST<-NA


### Loop to calculate the 4 variables
for(i in 1:nrow(species)){
  
# Create countries table
  countries<-subset(countriesTOT, countriesTOT$taxonid==species$taxonid[i])
  countries<-subset(countries, countries$origin %not in% c("Introduced", "Vagrant")) # Keep native, reintroduced, uncertain, assisted colonisation (remove introduced and vagrant)
  countriesUNIQ<-ddply(countries, .(Countries_IUCN), function(x){data.frame(N=nrow(x))})

# GDP
  countriesUNIQ$gdp<-ctr$GDP[match(countriesUNIQ$Countries_IUCN, ctr$Countries_IUCN)]
  
  species$med_GDP[i]<-median(countriesUNIQ$gdp, na.rm=T)
  species$nb_countries[i]<-nrow(countriesUNIQ)
  species$gdpOK[i]<-table(factor(is.na(countriesUNIQ$gdp), levels=c("TRUE", "FALSE")))["FALSE"]

# Research capacity
  countriesUNIQ$res.cap<-ctr$res.cap[match(countriesUNIQ$Countries_IUCN, ctr$Countries_IUCN)]
  
  species$res.cap[i]<-median(countriesUNIQ$res.cap, na.rm=T)
  species$rescapOK[i]<-table(factor(is.na(countriesUNIQ$res.cap), levels=c("TRUE", "FALSE")))["FALSE"]
  
  
# Corruption
  countriesUNIQ$corruption<-ctr$corruption[match(countriesUNIQ$Countries_IUCN, ctr$Countries_IUCN)]
  
  species$corruption[i]<-median(countriesUNIQ$corruption, na.rm=T)
  species$corruptionOK[i]<-table(factor(is.na(countriesUNIQ$corruption), levels=c("TRUE", "FALSE")))["FALSE"]

# Conflicts
  countriesUNIQ$conflicts<-ctr$conflict[match(countriesUNIQ$Countries_IUCN, ctr$Countries_IUCN)]
  
  species$conflictEXIST[i]<- sum(countriesUNIQ$conflicts)>0 # If sum>0, it means there is at least one country with conflict
  species$conflictONLY[i]<- 0 %not in%  countriesUNIQ$conflicts # If 0 not in the vector, it means there are only countries with conflict
  species$conflictMEDYEAR[i]<- median(countriesUNIQ$conflicts) # Takes the median number of year in conflicts
  species$conflictMEANYEAR[i]<- mean(countriesUNIQ$conflicts) # Takes the mean number of year in conflicts
  
}
  




### Present in zoos and aquarium
aqua<-read.table("0.Data/ZoosAquariums/datListIUCNRL_fishfromDalia.txt", header=T)
zoo<-read.csv("0.Data/ZoosAquariums/ZIMSData_v1.csv") ; names(zoo)[1]<-"species"
species$zoos<- (species$scientific_name %in% aqua$binSpecies | species$scientific_name %in% zoo$species)





### Trait data for mammals
cat("\n", "Start Traits extract")

combine_rep<-read.csv("0.Data/COMBINE/COMBINE_archives/trait_data_reported.csv")
combine_imp<-read.csv("0.Data/COMBINE/COMBINE_archives/trait_data_imputed.csv")

table(species$scientific_name[species$Vgroup =="Mammal"] %in% combine_imp$iucn2020_binomial)

# Mass (if NA, I take the average of the family)
species$bodymass_mammals<-combine_imp$adult_mass_g[match(species$scientific_name, combine_imp$iucn2020_binomial)] #table(species$category[species$Vgroup =="Mammal"], is.na(species$adult_mass_g[species$Vgroup =="Mammal"]))
mass_fam<-ddply(species[species$Vgroup=="Mammal",], .(family), function(x){data.frame(Mass_fam=mean(x$bodymass_mammals, na.rm=T))})
species$bodymass_mammals[is.na(species$bodymass_mammals) & species$Vgroup=="Mammal"]<-mass_fam$Mass_fam[match(species$family[is.na(species$bodymass_mammals) & species$Vgroup=="Mammal"], mass_fam$family)]

# Nocturnality (I assume mixed if NA)
species$nocturnal_mammals<-combine_imp$activity_cycle[match(species$scientific_name, combine_imp$iucn2020_binomial)]==1 # I put 2 and 3 together (partial nocturnal and full diurnal, as full nocturnal is majority and it means you can observe 2 and 3 in the day)
species$nocturnal_mammals[species$Vgroup=="Mammal"]<-replace(species$nocturnal_mammals[species$Vgroup=="Mammal"], is.na(species$nocturnal_mammals[species$Vgroup=="Mammal"])==T, "FALSE")

# Trait availability in COMBINE with only 1 of the binary diet variables
combine_rep$trait_availability<-rowSums(is.na(combine_rep[,c(7:34,45:60)])==F)
species$trait_availability_mammals<-combine_rep$trait_availability[match(species$scientific_name, combine_rep$iucn2020_binomial)]




### Trait data for amphibians
load("0.Data/Pablo_amphib/merged_all_databases_mod_V0.RData") # Data before imputation
load("0.Data/Pablo_amphib/trait_data_imp.RData") # Imputed data

# Take imputed body size 
trait_data_imp$ximp$scientific_name<-merged_all_databases_mod_V0$scientificName # The imputed data don't have the species name, but it's in the same order than in the other file, proof: plot(trait_data_imp$ximp$log10_Size_combo, merged_all_databases_mod_V0$log10_Size_combo)
species$SVL_amphibians<-trait_data_imp$ximp$log10_Size_combo[match(species$scientific_name, trait_data_imp$ximp$scientific_name)]


# Take trait availability (I chose one trait per category)
for(i in 1:ncol(merged_all_databases_mod_V0)) {merged_all_databases_mod_V0[,i]<-replace(merged_all_databases_mod_V0[,i], is.nan(merged_all_databases_mod_V0[,i]), NA)}
merged_all_databases_mod_V0$trait_availability<-rowSums(is.na(merged_all_databases_mod_V0[,c(62,63,67,70,73,80,81,83,84,87,88,90,92,94)])==F)
species$trait_availability_amphibians<-merged_all_databases_mod_V0$trait_availability[match(species$scientific_name, merged_all_databases_mod_V0$scientificName)]




############
### Save ###
############
write.csv(species, "1.Tables/Species.characteristics.datasufficiency.Script1a.csv", row.names=FALSE)

