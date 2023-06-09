

### Merge SIS files
RegionLS<-list()
HabLS<-list()
CountryLS<-list()
CritLS<-list()
DatesLS<-list()

groups<-c("actinopterygii", "actinopterygii2", "amphibia", "cephalaspidomorphi", "chondrichthyes", "mammal", "myxinidae", "odonata", "reptilia", "sarcopterygii")

for(GR in 1:length(groups)){
  RegionLS[[GR]]<-read.xlsx(xlsxFile=paste0("0.Data/SIS/1.Extracted files/REGIONINFORMATION_", groups[GR], ".xlsx"))
  HabLS[[GR]]<-read.xlsx(xlsxFile=paste0("0.Data/SIS/1.Extracted files/GENERALHABITATSSUBFIELD_", groups[GR], ".xlsx"))
  CountryLS[[GR]]<-read.xlsx(xlsxFile=paste0("0.Data/SIS/1.Extracted files/COUNTRYOCCURRENCESUBFIELD_", groups[GR], ".xlsx"))
  CritLS[[GR]]<-read.xlsx(xlsxFile=paste0("0.Data/SIS/1.Extracted files/REDLISTCRITERIA_", groups[GR], ".xlsx"))
  DatesLS[[GR]]<-read.xlsx(xlsxFile=paste0("0.Data/SIS/1.Extracted files/REDLISTASSESSMENTDATE_", groups[GR], ".xlsx"), detectDates=TRUE)
}

Region<-rbind(RegionLS[[1]], RegionLS[[2]], RegionLS[[3]], RegionLS[[4]], RegionLS[[5]], RegionLS[[6]], RegionLS[[7]], RegionLS[[8]], RegionLS[[9]], RegionLS[[10]])
Habs<-rbind(HabLS[[1]], HabLS[[2]], HabLS[[3]], HabLS[[4]], HabLS[[5]], HabLS[[6]], HabLS[[7]], HabLS[[8]], HabLS[[9]], HabLS[[10]])
Country<-rbind(CountryLS[[1]], CountryLS[[2]], CountryLS[[3]], CountryLS[[4]], CountryLS[[5]], CountryLS[[6]], CountryLS[[7]], CountryLS[[8]], CountryLS[[9]], CountryLS[[10]])
Crit<-rbind(CritLS[[1]], CritLS[[2]], CritLS[[3]], CritLS[[4]], CritLS[[5]], CritLS[[6]], CritLS[[7]], CritLS[[8]], CritLS[[9]], CritLS[[10]])
Dates<-rbind(DatesLS[[1]], DatesLS[[2]], DatesLS[[3]], DatesLS[[4]], DatesLS[[5]], DatesLS[[6]], DatesLS[[7]], DatesLS[[8]], DatesLS[[9]], DatesLS[[10]])

table(Crit$taxonid %in% species$taxonid) # The "FALSE" correspond to subspecies (e.g., "nov" or "subspecies")
table(species$taxonid %in% Crit$taxonid)
table(species$taxonid %in% Dates$taxonid)
table(species$taxonid %in% Country$taxonid)


### Are all assessments in this file global assessments (i.e., no assessments that are only regional)?
test<-as.data.frame.matrix(table(Region$assessmentid, Region$regions_value=="Global"))
table(test$'TRUE') # Should not be any 0



### Select latest assessment
Dates$day<-julian(Dates$value,as.Date("1900-01-01")) # Calculate number of days since 1900 (just to make sure that the "max" in the next line really captures the max)
CorrAssessment<-ddply(Dates, .(taxonid), function(x){
  v<-unique(x$assessmentid[x$day==max(x$day)])
  data.frame(Assessment=v[1], N_assess=length(v), Tot=paste0(v, collapse="/"))}) # Calculates the most recent assessment in the published assessments
CorrAssessment<-subset(CorrAssessment, CorrAssessment$taxonid %in% species$taxonid)

# Manual fix for 3 species with duplicate assessments
CorrAssessment$Assessment[CorrAssessment$taxonid=="42462"]="125997433" ; CorrAssessment$Assessment[CorrAssessment$taxonid=="178508"]="1537197" ; CorrAssessment$Assessment[CorrAssessment$taxonid=="4840"]="134828519"

species$Correct_Assessmentid<-CorrAssessment$Assessment[match(species$taxonid, CorrAssessment$taxonid)]
table(is.na(species$Correct_Assessmentid))


# Keep only most recent assessments
Habs2<-subset(Habs, Habs$assessmentid %in% species$Correct_Assessmentid)
Country2<-subset(Country, (Country$assessmentid %in% species$Correct_Assessmentid) & is.na(Country$countryoccurrencelookup_value)==F)
Crit2<-subset(Crit, Crit$assessmentid %in% species$Correct_Assessmentid) # Note that there can be several 2 lines for a species when 2 DD reasons

# Adapt habitat codes
Habs2$Code0<-Habs2$Code<-NA
for(HH in 1:nrow(Habs2)){Habs2$Code[HH]<-unlist(strsplit(Habs2$generalhabitatslookup_value[HH], " "))[1] ; Habs2$Code0[HH]<-unlist(strsplit(Habs2$Code[HH], "[.]"))[1]}

# Adapt countries codes
coun_classif<-read.csv("0.Data/Countries_classification.csv")
table(Country2$countryoccurrencelookup_value %in% coun_classif$Var1) # Should be 100% TRUE

Country2$Countries_IUCN<-coun_classif$Countries_IUCN[match(Country2$countryoccurrencelookup_value, coun_classif$Var1)]
table(is.na(Country2$Countries_IUCN)) # Should be 100% TRUE


# Save
saveRDS(Crit2, "0.Data/SIS/Habitats_adaptedfromSIS_REDLISTCRITERIA.rds")
saveRDS(Habs2, "0.Data/SIS/Habitats_adaptedfromSIS_GENERALHABITATSSUBFIELD.rds")
saveRDS(Country2, "0.Data/SIS/Countries_adaptedfromSIS_COUNTRYOCCURRENCESUBFIELD.rds")
saveRDS(CorrAssessment, "0.Data/SIS/CorrAss_adaptedfromSIS_PUBLICATIONDATE.rds")





#########################################
### MERGE SYNONYMOUS FROM SIS CONNECT ###
#########################################

# Charge and merge data 
syn1<-read.csv("0.Data/SIS/2.SISConnect/Amphibian/synonyms.csv")
syn2<-read.csv("0.Data/SIS/2.SISConnect/Reptile/synonyms.csv")
syn4<-read.csv("0.Data/SIS/2.SISConnect/Cephalaspidomorphi/synonyms.csv")
syn5<-read.csv("0.Data/SIS/2.SISConnect/Mammal/synonyms.csv")
syn6<-read.csv("0.Data/SIS/2.SISConnect/Myxini/synonyms.csv")
syn7<-read.csv("0.Data/SIS/2.SISConnect/Odonata/synonyms.csv")
syn8<-read.csv("0.Data/SIS/2.SISConnect/Sarcopterygii/synonyms.csv")
syn9<-read.csv("0.Data/SIS/2.SISConnect/Chondrichtyes/synonyms.csv")

synonyms<-rbind(syn1, syn2, syn4, syn5, syn6, syn7, syn8, syn9)

# Prepare year of publication
synonyms$year2<-synonyms$year<-NA
for(i in 1:nrow(synonyms)){
  v<-unlist(strsplit(synonyms$name[i]," "))
  synonyms$year[i]<-v[length(v)] %>% gsub("[^0-9]", "", .) %>% substr(., 1,4) 
  synonyms$year2[i]<-gsub("[^0-9]", "", synonyms$name[i])
}

synonyms$year[synonyms$year==""]<-synonyms$year2[synonyms$year==""] # Those where the first trick did not work, I use the second (ie taking all numbers from the name and checking manually if some did not work)
table(synonyms$year[synonyms$year %not in% 1500:2022])

# I remove the few that are not in 1500:2022, I will have to go check species one by one if I end up keeping this variable
synonyms$year[synonyms$year %not in% 1500:2022 | synonyms$year =='']<-NA 
hist(as.numeric(synonyms$year))

# Save
write.csv(synonyms, "1.Tables/Synonymous.processed.table.csv", row.names=F)



