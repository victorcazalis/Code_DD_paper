
### Charge functions that are needed in the script
source("2.Scripts/03.Functions.analyses.sufficiency.R")

  
  
library(plyr) ; library(dplyr)  ; library(ggcorrplot) ; library(rfUtilities) ; library(gridExtra) ; library(grid) ; library(ggpubr) ; library(stringr)

# Prepare table for validation
perf<-data.frame(expand.grid(Metric=c("Accuracy", "Sensitivity", "Specificity", "TSS"), Group=c("Odonata", "Fish", "Amphibian", "Reptile", "Mammal"), Domain=c("Marine", "Terrestrial"), Threshold=c(0.5,0.65,"EQUAL", "MAXTSS", "SENS90"), Threshold.QTT=NA)) ; perf$Time<-perf$Value<-NA

# Prepare empty list for plots
Plot_tot<-Imp_plot<-Cov_plot<-G_Cor<-Imp.var.plot<-Plot_valBin<-M0_pred<-SpMOD_stored<-vector(mode="list", length=6)





### LOOP WITH ALL GROUPS
for(Domain in c("Marine", "Terrestrial")){
  if(Domain=="Marine"){GROUPS<-c("Fish")}
  if(Domain=="Terrestrial"){GROUPS<-c("Fish", "Reptile", "Mammal", "Odonata", "Amphibian")}
  
  for(Group in GROUPS){
    
    cat(paste0("\n", "\n", "\n", "Starting analysis with: ", Group, " ", Domain, "\n"))


# For Odonata I created some range maps (Hydro-to-create) but decided not to keep them because they are not of good quality, so I eventually exclude those species
if(Group=="Odonata"){
  species<-readRDS("1.Tables/Species.characteristics.Script1d.FW_ODONATA.rds")
  species<-subset(species, is.na(species$Dist_quality) | species$Dist_quality != "Hydro_to_create")
}
    
if(Group=="Reptile"){species<-readRDS("1.Tables/Species.characteristics.Script1d.REPTILES.rds")}
if(Group=="Amphibian"){species<-readRDS("1.Tables/Species.characteristics.Script1d.AMPHIBIANS.rds")}
if(Group=="Fish"){species<-readRDS("1.Tables/Species.characteristics.Script1d.FISH.rds")}
if(Group=="Mammal"){species<-readRDS("1.Tables/Species.characteristics.Script1d.MAMMALS.rds")}


# Index to order plots
IND<-revalue(Group, c("Mammal"="1", "Reptile"="2", "Amphibian"="3", "Fish"="4", "Odonata"="6")) %>% as.numeric(.)
if(Domain=="Marine"){IND<-5}

# Remove DD species with no distribution and remove extinct species
species<-subset(species, species$Dist_exist==T & !species$category %in% c("EX", "EW"))

# Keep only domain (marine or terrestrial + freshwater)
if(Domain=="Marine"){species<-subset(species, species$marine_system==TRUE)}
if(Domain=="Terrestrial" & Group!="Reptile"){species<-subset(species, species$terrestrial_system==TRUE | species$freshwater_system==TRUE)}



################################
### FINETUNE SOME COVARIATES ###
################################

# When no data within range Sampling NA, so I fix to 0
species$Sampling_MED[which(is.na(species$Sampling_MED))]<-0

# When no GBIF data or no data with correct coordinates for fetch function (Fetch Null), propempty is fixed to 1
species$gbifSP_propempty[which(species$gbifSP_check=="Fetch Null" | species$nb_GBIFgeo==0)]<-1
species$gbifSP_propemptyASS[which(species$gbifSP_check=="Fetch Null" | species$nb_GBIFgeoASS==0)]<-1

# Transform propempty into coverage
species$gbifSP_coverage<-1-species$gbifSP_propempty
species$gbifSP_coverageASS<-1-species$gbifSP_propemptyASS

# Keep magnitude of range size only
species$range_magnitude<-trunc(log10(species$range_km2))

# Create a DD variable
species$DD<-species$category=="DD" ; table(species$DD)

# Time since description
species$time_desc<-2022-species$description_year
species$time_descASS<-species$published_year-species$description_year

# Relative rural
species$rur.prop<-species$hrur/species$hpop
species$rur.prop[which(species$hpop==0)]<-1

# For fishes, GBIF data is the maximum between GBIF and OBIS
if(Group == "Fish"){
  species$nb_GBIFgeo<-apply(species[,c("nb_GBIFgeo", "nb_OBIS")], 1, "max", na.rm=T)
  species$nb_GBIFgeoASS<-apply(species[,c("nb_GBIFgeoASS", "nb_OBISASS")], 1, "max", na.rm=T)
}

# Taxonomy
if(Group %in% c("Reptile", "Amphibian", "Mammal")){species$taxo<-species$order; species$taxo_valid<-species$family}
if(Group %in% c("Odonata")){species$taxo<-species$order; species$taxo_valid<-species$family} # Only one order so no control for taxonomy here
if(Group %in% c("Fish")){species$taxo<-species$class; species$taxo_valid<-species$family} 

# Habitat
if(Group %in% c("Amphibian", "Reptile", "Mammal", "Odonata")){
  species$Habitat<-revalue(
    as.factor(paste0(species$Bin_forest, species$N.hab1==1)), 
    c("FALSEFALSE"="NF_generalist", "FALSETRUE"="NF_specialist", "TRUEFALSE"="Forest_generalist", "TRUETRUE"="Forest_specialist")
  ) %>% factor(., c(levels(.), "RockCave"))
  species$CaveRocky<-as.numeric(species$Bin_cave)+as.numeric(species$Bin_rocky)
  species$Habitat[species$CaveRocky>0 & species$CaveRocky==species$N.hab1]<-"RockCave"  # Species with some rocky habitats + no other habitats are considered "RockCave"
  species$CaveRocky<-NULL
}




###########################
### RANDOM FOREST MODEL ###
###########################

library(ranger)


### Fix the four temporal covariates to their value by the time of assessment
species$time_descMOD<-species$time_descASS
species$nb_GBIFgeoMOD<-species$nb_GBIFgeoASS
species$gbifSP_coverageMOD<-species$gbifSP_coverageASS
species$WOSMOD<-species$WOSASS


### Create set of variables (general + terrestrial / marine)
VARS<-c(
  "nb_GBIFgeoMOD",
  "time_descMOD",
  "med_GDP",
  "conflictMEDYEAR",
  "range_magnitude", 
  "Sampling_MED", "gbifSP_coverageMOD",
  "DDness", 
  "WOSMOD",
  "zoos",
  "taxo"
)

if(Domain=="Terrestrial"){VARS<-c(VARS,
   "realmsTERR_mode",
   "acc_terrMED",
   "roadsQ90",
   "hpop", "rur.prop",
   "alt_med"
)}


if(Domain=="Marine"){VARS<-c(VARS,
  "realmsMAR_mode", 
  "acc_marMED",
  "fishingMED",
  "Depth_min"
)}

if(Group=="Mammal"){VARS<-c(VARS, "trait_availability_mammals", "bodymass_mammals", "nocturnal_mammals", "RLA")}  
if(Group=="Amphibian"){VARS<-c(VARS, "trait_availability_amphibians", "SVL_amphibians")}
if(Group!="Fish" & Domain=="Terrestrial"){VARS<-c(VARS, "Habitat")}
if(Group %in% c("Reptile", "Fish")){VARS<-c(VARS, "RLA")}


### Create complete dataframe with only variables of interest (DD and covariates)
species_MOD<-species[,names(species) %in% c("DD", "taxonid", "nb_GBIFgeo", "nb_GBIFgeoASS", "time_desc", "time_descASS", "taxo_valid", "gbifSP_coverage", "gbifSP_coverageASS", "WOS", "WOSASS", VARS)] %>% .[complete.cases(.),]
cat("Completeness: ", nrow(species_MOD), " / ", nrow(species))
if("taxo" %in% names(species_MOD)){species_MOD$taxo<-droplevels(as.factor(species_MOD$taxo))}



## Check covariate correlation
VARS_qtt<-names(unlist(lapply(species_MOD[,names(species_MOD) %in% VARS], is.numeric)) %>% .[.==TRUE])
df_cor<-species_MOD[,names(species_MOD) %in% VARS_qtt] ; names(df_cor)<-Pretty_name(names(df_cor))
Cor<-cor(df_cor)
G_Cor[[IND]]<-ggcorrplot(Cor, type="upper", tl.srt=60, colors=c("red", "white", "blue"))+ggtitle(paste0(Group, " (", ifelse(Domain=="Terrestrial" & Group=="Fish", "Freshwater", Domain), "), max = ", round(max(abs(Cor[Cor!=1])),2)))
cowplot::save_plot(paste0("3.Figures/Correlation/Correlation.plot.", Group, ".", Domain, ".png"), G_Cor[[IND]], base_width=8, base_height=8)






### Prepare DD variable
species_MOD$DD<-factor(revalue(as.character(species_MOD$DD), c("TRUE"="DD", "FALSE"="DS")), levels=c("DS", "DD")) # The model uses the first factor as reference

# Reformat categorical variables as factor (ranger does something weird with characters)
if(Domain=="Marine"){species_MOD$realmsMAR_mode <- factor(species_MOD$realmsMAR_mode)}
if(Domain=="Terrestrial"){species_MOD$realmsTERR_mode <- factor(species_MOD$realmsTERR_mode)}
species_MOD$taxo <- factor(species_MOD$taxo)
if("Habitat" %in% VARS){species_MOD$Habitat <- factor(species_MOD$Habitat)}
species_MOD$taxo_valid <- factor(species_MOD$taxo_valid)
species_MOD$taxo <- factor(species_MOD$taxo)
if("RLA" %in% VARS){species_MOD$RLA <- factor(species_MOD$RLA)}
SpMOD_stored[[IND]]<-species_MOD


### Run model (twice because for predictions it's recommended to run without impurity_corrected)
eval(parse(text= 
             paste0("M0<-ranger(DD ~", paste(VARS, collapse="+"), ", data=species_MOD, class.weights=prop.table(table(species_MOD$DD))[c(2,1)], num.trees=1000, importance='impurity_corrected', probability=T) ") 
)) 
eval(parse(text= 
             paste0("M0_pred[[IND]]<-ranger(DD ~", paste(VARS, collapse="+"), ", data=species_MOD, class.weights=prop.table(table(species_MOD$DD))[c(2,1)], num.trees=1000, probability=T) ") 
)) 


## Run Validation
species_MOD$taxo_valid<-droplevels(as.factor(species_MOD$taxo_valid))
species_MOD$Predict_validation<-NA

for(FAM in levels(species_MOD$taxo_valid)){
  tryCatch({
    cat(paste0("Starting family: ", FAM, "(N=", nrow(species_MOD[species_MOD$taxo_valid==FAM,]), ")", "\n"))
    species_val<-subset(species_MOD, species_MOD$taxo_valid != FAM)

    eval(parse(text=paste0("M_val<-ranger(DD ~", paste(VARS, collapse="+"), ", data=species_val, class.weights=prop.table(table(species_val$DD))[c(2,1)], num.trees=1000, importance='impurity_corrected', probability=T) ") ))

    species_MOD$Predict_validation[species_MOD$taxo_valid==FAM]<-predict(M_val, data=species_MOD[species_MOD$taxo_valid==FAM,], type="response")$predictions[,"DS"]

  } ,error=function(e){cat(paste0("Bug at family", FAM))})
}

# Calculate performance and plot
perf<-DD_perf(THRE=0.5, perf.df=perf)
perf<-DD_perf(THRE=0.65, perf.df=perf)
perf<-DD_perf(THRE="EQUAL", perf.df=perf)
perf<-DD_perf(THRE="MAXTSS", perf.df=perf)
perf<-DD_perf(THRE="SENS90", perf.df=perf)

Bin1<-ggplot()+
  geom_histogram(data=species_MOD, aes(x=Predict_validation, fill=DD))+
  geom_vline(data=perf[perf$Group==Group & perf$Domain==Domain & perf$Threshold %in% c("MAXTSS", "SENS90"),], aes(xintercept=Threshold.QTT, col=Threshold), size=1.5, show.legend=F)+
  scale_colour_manual(values=c("#fdae61", "#d73027"))+
  scale_fill_manual(values=c("#abd9e9", "#4575b4"), name="")+
  xlab("Probability of being DS")+
  theme_minimal()


Thre_MAXTSS<-perf$Threshold.QTT[perf$Group==Group & perf$Domain==Domain & perf$Metric=="Accuracy" & perf$Threshold=="MAXTSS"]
Thre_SENS90<-perf$Threshold.QTT[perf$Group==Group & perf$Domain==Domain & perf$Metric=="Accuracy" & perf$Threshold=="SENS90"]
tab_MAXTSS<-as.data.frame.matrix(table(factor(species_MOD$DD, c("DD", "DS")), revalue(as.character(species_MOD$Predict_validation>Thre_MAXTSS), c("FALSE"="Predicted_DD", "TRUE"="Predicted_DS"))))
confusion_MAXTSS<-ggtexttable(tab_MAXTSS, theme=ttheme("mOrange")) %>%
  tab_add_title(., text="Maximised TSS", hjust=-0.8,  padding = unit(1, "line")) %>%
  table_cell_bg(., row=3, column=2, fill="#fdc79aff", col="white") %>%
  table_cell_bg(., row=4, column=3, fill="#fdc79aff", col="white") %>%
  table_cell_bg(., row=4, column=2, fill="#ffe9d7ff", col="white") %>%
  table_cell_bg(., row=3, column=3, fill="#ffe9d7ff", col="white")

tab_SENS90<-as.data.frame.matrix(table(factor(species_MOD$DD, c("DD", "DS")), revalue(as.character(species_MOD$Predict_validation>Thre_SENS90), c("FALSE"="Predicted_DD", "TRUE"="Predicted_DS"))))
confusion_SENS90<-ggtexttable(tab_SENS90, theme=ttheme("mRed")) %>%
  tab_add_title(., text="High sensitivity", hjust=-0.8,  padding = unit(1, "line")) %>%
  table_cell_bg(., row=3, column=2, fill="#ef9596ff", col="white") %>%
  table_cell_bg(., row=4, column=3, fill="#ef9596ff", col="white") %>%
  table_cell_bg(., row=4, column=2, fill="#f7dbdaff", col="white") %>%
  table_cell_bg(., row=3, column=3, fill="#f7dbdaff", col="white")

Plot_valBin[[IND]]<-gridExtra::grid.arrange(confusion_SENS90, confusion_MAXTSS, Bin1,
                        layout_matrix=matrix(c(4,1,4,2,4,3,3,3,3,3,3,3,3,3,3), ncol=5, byrow=T),
                        top=textGrob(paste0(Group, ifelse(Domain=="Terrestrial", "", "-Marine")),
                                     gp=gpar(fontsize=20,font=2)))

cowplot::save_plot(paste0("3.Figures/Plot.Validation.Binarisation.", Group, ".", Domain, ".png"), Plot_valBin[[IND]], base_width=7.5, base_height=5)


# Save
write.csv(species_MOD, paste0("1.Tables/Validation/Validation.table.", Group, ".", Domain, ".csv"), row.names=F)
write.csv(perf, "1.Tables/Performance.saved.csv", row.names = F)



####################
### PLOTS MODELS ###
####################

### Plot variable importance
Imp<-ranger::importance(M0)
Imp.Vars<-data.frame(Variable=factor(names(Imp), levels=names(Imp)[order(Imp, decreasing=T)]), Importance=Imp, row.names=1:length(Imp))
Imp.Vars<-Imp.Vars[order(Imp.Vars$Importance, decreasing=T),]

Imp.var.plot[[IND]]<-ggplot(Imp.Vars)+
  geom_bar(aes(x=Pretty_name(Variable), y=Importance), stat="identity")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle(paste0("Variable importance for ", Group, " (", Domain, ")"))

cowplot::save_plot(paste0("3.Figures/Variable.Importance.", Group, ".", Domain, ".png"), Imp.var.plot[[IND]], base_width=8, base_height=4)



### Plot partial dependence
par.dep<-ParDepCalc(mod=M0, df=species_MOD)

Cov_plot[[IND]]<-grid.arrange(
  perfectPartialPlot(Var2plotNUM=1, par.dep=par.dep, df=species_MOD),
  perfectPartialPlot(Var2plotNUM=2, par.dep=par.dep, df=species_MOD),
  perfectPartialPlot(Var2plotNUM=3, par.dep=par.dep, df=species_MOD),
  perfectPartialPlot(Var2plotNUM=4, par.dep=par.dep, df=species_MOD),
  nrow=1, top=textGrob(paste0(Group, ifelse(Domain=="Terrestrial", "", "-Marine")), gp=gpar(fontsize=25,font=1)))

Imp_plot[[IND]]<-ggplot(Imp.Vars[1:4,])+
  geom_segment(x=1, xend=4, y=1, yend=1, col="gray70")+
  geom_point(aes(y=1, x=Variable, size=Importance), show.legend=F, fill="gray70", shape=21)+
  ylim(c(0.5,1.2))+
  scale_size_continuous(range=c(1,12))+
  theme_void()

Plot_tot[[IND]]<-grid.arrange(Cov_plot[[IND]], Imp_plot[[IND]], layout_matrix=matrix(c(1,1,1,2), ncol=1))

cowplot::save_plot(paste0("3.Figures/Covariate.effects.", Group, ".", Domain, ".png"), Plot_tot[[IND]], base_width=12, base_height=4)






#######################################
### PREDICT DATA SUFFICIENT SPECIES ###
#######################################

DD<-subset(species_MOD, species_MOD$DD == "DD")

# At assessment
DD$Prob.DS.ASS<- predict(M0_pred[[IND]], data=DD, type="response", num.trees=1000)$predictions[,"DS"]

# Replace by current values
DD$nb_GBIFgeoMOD<-DD$nb_GBIFgeo
DD$time_descMOD<-DD$time_desc
DD$gbifSP_coverageMOD<-DD$gbifSP_coverage
DD$WOSMOD<-DD$WOS

# Calculate pDS and iDS
DD$Prob.DS<- predict(M0_pred[[IND]], data=DD, type="response", num.trees=1000)$predictions[,"DS"]
DD$Prob.increase<- (DD$Prob.DS-DD$Prob.DS.ASS) %>% replace(., .<0, 0)


#DD<-DD[order(DD$Prob.increase, decreasing=T),]
DD$Vgroup<-Group ; DD$Domain<-Domain
write.csv(DD, paste0("1.Tables/DD_predictions/DD.proba.predicted.", Group, ".", Domain, ".csv"), row.names=F)
}
}













#########################
### PLOT MULTISPECIES ###
#########################

### Model performance
perf$Value<-as.numeric(perf$Value)
perf_plot<-subset(perf, is.na(perf$Value)==F & perf$Threshold %in% c("MAXTSS", "SENS90"))

# Prepare xaxis
perf_plot$Group_to_plot<-ifelse(perf_plot$Group=="Fish", paste0(perf_plot$Group, " (", perf_plot$Domain, ")"), paste0(perf_plot$Group)) %>%
  replace(., .=="Fish (Terrestrial)", "Fish (Freshwater)") %>%
  factor(., c("Mammal", "Reptile", "Amphibian", "Fish (Freshwater)", "Fish (Marine)", "Odonata"))


plot_perf<-function(MET, YMIN, YMAX){
  ggplot(perf_plot[perf_plot$Metric==MET,])+
  geom_point(aes(x=Group_to_plot, y=Value, col=Threshold), size=4, show.legend=F)+
  ylim(c(YMIN, YMAX))+
  ggtitle(revalue(MET, c("Sensitivity"="(a) Sensitivity", "Specificity"="(b) Specificity", "TSS"="(c) TSS")))+
  xlab("")+ylab("")+
  theme_minimal() %+replace%  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_manual(values=c("#fdae61", "#d73027"))}


Perf<-grid.arrange(
  plot_perf("Sensitivity", 0, 1),
  plot_perf("Specificity", 0, 1),
  plot_perf("TSS", -1, 1),
  ncol=3)

cowplot::save_plot("3.Figures/0.Main/Fig.1.Performance.plot.svg", Perf, base_width=7.5, base_height=4.2)








### Covariate
Cov_plot_glob<-grid.arrange(
  Plot_tot[[1]], 
  Plot_tot[[2]], 
  Plot_tot[[3]], 
  Plot_tot[[4]], 
  Plot_tot[[5]], 
  Plot_tot[[6]],
  layout_matrix=matrix(c(1,1,1,1,1,7,2,2,2,2,2,3,3,3,3,3,7,4,4,4,4,4,5,5,5,5,5,7,6,6,6,6,6), nrow=3, byrow=T)
)
cowplot::save_plot(paste0("3.Figures/0.Main/Covariate.effects.AllGroups.svg"), Cov_plot_glob, base_width=15, base_height=9.5)







### Predict
DD_tot_raw<-read.csv(paste0("1.Tables/DD_predictions/", list.files(paste0("1.Tables/DD_predictions"))[1]))
for(i in 2:length(list.files(paste0("1.Tables/DD_predictions")))){DD_tot_raw<-rbind.fill(DD_tot_raw, read.csv(paste0("1.Tables/DD_predictions/", list.files(paste0("1.Tables/DD_predictions"))[i])) )}

# Merging fishes together
DD_tot<-ddply(DD_tot_raw, .(taxonid), function(x){data.frame(
  Prob.DS=max(x$Prob.DS),
  Prob.increase=x$Prob.increase[which(x$Prob.DS==max(x$Prob.DS))],
  Vgroup=x$Vgroup[1]
)})

# Charge AOH results
AOH_tot<-read.csv(paste0("1.Tables/", list.files("1.Tables")[grepl("AOH.", list.files("1.Tables"))][1]))
for(i in 2:4){AOH_tot<-rbind.fill(AOH_tot, read.csv(paste0("1.Tables/", list.files("1.Tables")[grepl("AOH.", list.files("1.Tables"))][i])) )}

# Extrapolate if 3 generations time brings before 2000 (Forest Loss) or 1992 (CCI)
AOH_tot$AOHlost_extrapol<-AOH_tot$AOHlost ; AOH_tot$ForestLoss_extrapol<-AOH_tot$ForestLoss
for(i in 1:nrow(AOH_tot)){
  AOH_tot$ref_year[i]<- 2020-max(10, round(3* replace(AOH_tot$GL[i], is.na(AOH_tot$GL[i]), 1) ))
  
  #CCI
  if(AOH_tot$ref_year[i]<1992){AOH_tot$AOHlost_extrapol[i]<-AOH_tot$AOHlost[i] * (1+(1992-AOH_tot$ref_year[i]) / (2020-1992))}
  
  #Forest
  if(AOH_tot$ref_year[i]<2000){AOH_tot$ForestLoss_extrapol[i]<-AOH_tot$ForestLoss[i] * (1+(2000-AOH_tot$ref_year[i]) / (2020-2000))}
  
} 

table(AOH_tot$ref_year, AOH_tot$GL)

# Bring in results dataset
DD_tot$AOHlost<-AOH_tot$AOHlost_extrapol[match(DD_tot$taxonid, AOH_tot$taxonid)]
DD_tot$Forestloss<-AOH_tot$ForestLoss_extrapol[match(DD_tot$taxonid, AOH_tot$taxonid)]

# Calculate PrioDS
DD_tot$Dist_11<-sqrt((1-DD_tot$Prob.DS)^2 + (1-DD_tot$Prob.increase)^2) 
DD_tot$PrioDS<- 1- sqrt(0.5)*DD_tot$Dist_11
DD_tot$PrioDS<-replace(DD_tot$PrioDS, (DD_tot$AOHlost>0.2 | DD_tot$Forestloss>0.2) , 1)



# Plot function
plot_isoc<-function(GR){
DD_GR<-subset(DD_tot, DD_tot$Vgroup==GR)
DD_GR<-subset(DD_GR, is.na(PrioDS)==F)
QUANT<-quantile(DD_GR$PrioDS, probs=c(0,0.5,0.75,0.9,1))
DD_GR$Dist_cat<-cut(DD_GR$PrioDS[DD_GR$Vgroup==GR], breaks=QUANT, include.lowest=T)

Isoclines<-data.frame(X=seq(0,1,0.01))
Isoclines$Y1<- 1-sqrt(2*(1-QUANT[2])^2- (1-Isoclines$X)^2)
Isoclines$Y2<- 1-sqrt(2*(1-QUANT[3])^2- (1-Isoclines$X)^2)
Isoclines$Y3<- 1-sqrt(2*(1-QUANT[4])^2- (1-Isoclines$X)^2)
Isoclines[, c("Y1", "Y2", "Y3")]<-replace(Isoclines[, c("Y1", "Y2", "Y3")], Isoclines[, c("Y1", "Y2", "Y3")]<0, NA)


ggplot()+
  geom_line(data=Isoclines, aes(x=X, y=Y1), col="#de77ae", linetype="dashed")+
  geom_line(data=Isoclines, aes(x=X, y=Y2), col="#c51b7d", linetype="dashed")+
  geom_line(data=Isoclines, aes(x=X, y=Y3), col="#8e0152", linetype="dashed")+
  geom_point(data=DD_GR, aes(x=Prob.DS, y=Prob.increase, col=Dist_cat), show.legend=F)+
  geom_point(data=DD_GR[DD_GR$AOHlost>0.2 | DD_GR$Forestloss>0.2,], aes(x=Prob.DS, y=Prob.increase), col="black", fill="white", shape=21, size=2, show.legend=F)+
  geom_point(data=DD_GR[DD_GR$AOHlost>0.3 | DD_GR$Forestloss>0.3,], aes(x=Prob.DS, y=Prob.increase), col="black", size=2)+
  scale_colour_manual(values=(c("gray75", "#de77ae", "#c51b7d", "#8e0152")), drop=FALSE)+
  scale_fill_manual(values=(c("gray75", "#de77ae", "#c51b7d", "#8e0152")), drop=FALSE)+
  ylim(c(min(DD_GR$Prob.increase), max(DD_tot$Prob.increase)))+
  xlab("")+ylab("")+
  theme_minimal() %+replace% theme(plot.title=element_text(hjust=0.5, size=12, face="bold"))+ 
  ggtitle(GR)
}


Isoc_Mam<-plot_isoc(GR="Mammal")
Isoc_Rept<-plot_isoc(GR="Reptile")
Isoc_Amph<-plot_isoc(GR="Amphibian")
Isoc_Fish<-plot_isoc(GR="Fish")
Isoc_Odon<-plot_isoc(GR="Odonata")

 
Plot_Predict<-grid.arrange(
  Isoc_Mam,
  Isoc_Rept,
  Isoc_Amph,
  Isoc_Fish,
  Isoc_Odon,
  ncol=3, left=" ", bottom=" ")

cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.svg"), Plot_Predict, base_width=8.5, base_height=5.1)

cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.Mammal.svg"), Isoc_Mam, base_width=3.3, base_height=3)
cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.Reptile.svg"), Isoc_Rept, base_width=3.3, base_height=3)
cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.Amphibian.svg"), Isoc_Amph, base_width=3.3, base_height=3)
cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.Fish.svg"), Isoc_Fish, base_width=3.3, base_height=3)
cowplot::save_plot(paste0("3.Figures/0.Main/Fig.3.Predicted.isoclines.Odonata.svg"), Isoc_Odon, base_width=3.3, base_height=3)






### Save all data (covariates + predictions) in Extended_data_S1.csv
# Charge groups
dfTM<-read.csv("1.Tables/Validation/Validation.table.Mammal.Terrestrial.csv") ; dfTM$Domain="Terrestrial" ; dfTM$Group="Mammal" 
dfTR<-read.csv("1.Tables/Validation/Validation.table.Reptile.Terrestrial.csv") ; dfTM$Domain="Terrestrial" ; dfTM$Group="Reptile" 
dfTA<-read.csv("1.Tables/Validation/Validation.table.Amphibian.Terrestrial.csv") ; dfTM$Domain="Terrestrial" ; dfTM$Group="Amphibian" 
dfTF<-read.csv("1.Tables/Validation/Validation.table.Fish.Terrestrial.csv") ; dfTM$Domain="Terrestrial" ; dfTM$Group="Fish" 
dfMF<-read.csv("1.Tables/Validation/Validation.table.Fish.Marine.csv") ; dfTM$Domain="Marine" ; dfTM$Group="Fish" 
dfTO<-read.csv("1.Tables/Validation/Validation.table.Odonata.Terrestrial.csv") ; dfTM$Domain="Terrestrial" ; dfTM$Group="Odonata" 

# Merge groups
df<-rbind.fill(dfTM, dfTR, dfTA, dfTF, dfMF, dfTO)
df$Predict_validation<-NULL

# Add scientific name
search<-readRDS("1.Tables/Species.merged.april22.rds")
df$scientific_name<-search$scientific_name[match(df$taxonid, search$taxonid)]

# Order columns
df<-df[, c(which(names(df)=="Domain"), which(names(df)=="Group"), which(names(df)=="taxonid"), which(names(df)=="scientific_name"), which(names(df)=="DD"), which(!names(df) %in% c("taxonid", "Domain", "Group", "DD", "scientific_name")))]

# Add predictions from DD_tot
df$pDS<-DD_tot$Prob.DS[match(df$taxonid, DD_tot$taxonid)]
df$iDS<-DD_tot$Prob.increase[match(df$taxonid, DD_tot$taxonid)]
df$AOHloss_CCI<- (-1)*DD_tot$AOHlost[match(df$taxonid, DD_tot$taxonid)] # -1 to come from a loss (what I calculate in the code) to a change (what I talk about in the paper)
df$AOHloss_Forest<- (-1)*DD_tot$Forestloss[match(df$taxonid, DD_tot$taxonid)]# -1 to come from a loss (what I calculate in the code) to a change (what I talk about in the paper)
df$PrioDS<-DD_tot$PrioDS[match(df$taxonid, DD_tot$taxonid)]

# Save
write.csv(df, "3.Figures/0.Supp/Extended_data_S1.csv", row.names=F)








### Correlation
Corr_plot_glob<-grid.arrange(
  G_Cor[[1]], 
  G_Cor[[2]], 
  G_Cor[[3]], 
  G_Cor[[4]], 
  G_Cor[[5]], 
  G_Cor[[6]],
  ncol=3
)
cowplot::save_plot(paste0("3.Figures/0.Supp/Correlation.checks.AllGroups.png"), Corr_plot_glob, base_width=20, base_height=13)





### Variable importance
Imp.var_glob<-grid.arrange(
  Imp.var.plot[[1]], 
  Imp.var.plot[[2]], 
  Imp.var.plot[[3]], 
  Imp.var.plot[[4]], 
  Imp.var.plot[[5]], 
  Imp.var.plot[[6]],
  ncol=3
)
cowplot::save_plot(paste0("3.Figures/0.Supp/Variable.importance.AllGroups.svg"), Imp.var_glob, base_width=15, base_height=10)





### Plot validation binarisation
Plot_valbin_glob<-grid.arrange(
  Plot_valBin[[1]], 
  Plot_valBin[[2]], 
  Plot_valBin[[3]], 
  Plot_valBin[[4]], 
  Plot_valBin[[5]], 
  Plot_valBin[[6]],
  ncol=2
)
cowplot::save_plot(paste0("3.Figures/0.Supp/Plot_binarisation.AllGroups.svg"), Plot_valbin_glob, base_width=15, base_height=20)



DD_tot$Dist_11<-sqrt((1-DD_tot$Prob.DS)^2 + (1-DD_tot$Prob.increase)^2) 
DD_tot$PrioDS<- 1- sqrt(0.5)*DD_tot$Dist_11


### Compare predictions fishes
Fterr<-read.csv("1.Tables/DD_predictions/DD.proba.predicted.Fish.Terrestrial.csv")
Fmar<-read.csv("1.Tables/DD_predictions/DD.proba.predicted.Fish.Marine.csv")

Fterr<-subset(Fterr, Fterr$taxonid %in% Fmar$taxonid)
Fterr$Prob.DS.MARINE<-Fmar$Prob.DS[match(Fterr$taxonid, Fmar$taxonid)]
Fterr$Prob.increase.MARINE<-Fmar$Prob.increase[match(Fterr$taxonid, Fmar$taxonid)]

Fterr$PrioDS<- 1- sqrt(0.5*((1-Fterr$Prob.DS)^2 + (1-Fterr$Prob.increase)^2)) 
Fterr$PrioDS.MARINE<- 1- sqrt(0.5*((1-Fterr$Prob.DS.MARINE)^2 + (1-Fterr$Prob.increase.MARINE)^2)) 


Gcomp0<-ggplot(Fterr)+
  geom_point(aes(x=PrioDS, y=PrioDS.MARINE))+
  xlab("According to the terrestrial model")+ylab("According to the marine model")+
  ggtitle("PrioDS")+
  theme_bw()

Gcomp1<-ggplot(Fterr)+
  geom_point(aes(x=Prob.DS, y=Prob.DS.MARINE))+
  xlab("According to the terrestrial model")+ylab("According to the marine model")+
  ggtitle("pDS")+
  theme_bw()

Gcomp2<-ggplot(Fterr)+
  geom_point(aes(x=Prob.increase, y=Prob.increase.MARINE))+
  xlab("According to the terrestrial model")+ylab("According to the marine model")+
  ggtitle("δpDS")+
  theme_bw()

cowplot::save_plot(paste0("3.Figures/0.Supp/Fish_comparison_terrestrial_marine.png"), grid.arrange(Gcomp0, Gcomp1, Gcomp2, ncol=3), base_width=11, base_height = 5)






### Map results
source("2.Scripts/03.Mapping.results.R")







### GBIF knee
Plot_knee=list()

for(IND in 1:6){
  
  species_MOD=SpMOD_stored[[IND]]
  tryCatch({species_MOD$realmsMAR_mode <- factor(species_MOD$realmsMAR_mode)}, error=function(e){cat("No var")})
  tryCatch({species_MOD$realmsTERR_mode <- factor(species_MOD$realmsTERR_mode)}, error=function(e){cat("No var")})
  tryCatch({species_MOD$Habitat <- factor(species_MOD$Habitat)}, error=function(e){cat("No var")})
  tryCatch({species_MOD$taxo <- factor(species_MOD$taxo)}, error=function(e){cat("No var")})
  tryCatch({species_MOD$RLA <- factor(species_MOD$RLA)}, error=function(e){cat("No var")})

  par.dep<-partial(M0_pred[[IND]], pred.data = species_MOD, pred.var="nb_GBIFgeoMOD", prob=T, quantiles=T, probs=seq(0,0.95,0.01)) 
  names(par.dep)<-c("x", "y")

  mod<-chngpt::chngptm(formula.1=y ~ 1, formula.2= ~x, par.dep, type="M10", family="gaussian") 
  par.dep$Pred=predict(mod, par.dep)

  Plot_knee[[IND]]=ggplot(par.dep[par.dep$x<100,])+
    geom_line(aes(x=x, y=Pred), col="black", size=2)+
    geom_line(aes(x=x, y=y), col="#33a02c", size=1)+
    xlab("Number of GBIF records")+ylab("Probability of being DS")+
    ggtitle(c("Mammal", "Reptile", "Amphibian", "Fish", "Fish - marine", "Odonata")[IND])+
    labs(subtitle=paste0("Knee: ", mod$coefficients["chngpt"]))+
    theme_minimal()
  
  cat(IND)
}

Plot_knee_glob<-grid.arrange(
  Plot_knee[[1]], Plot_knee[[2]], Plot_knee[[3]], Plot_knee[[4]], Plot_knee[[5]], Plot_knee[[6]], ncol=3
)
cowplot::save_plot(paste0("3.Figures/0.Supp/Plot_knee_GBIF.svg"), Plot_knee_glob, base_width=8, base_height=4.5)








### Local characterisation for a species (for details in Fig. 4)
predict.function <- function(model, new_observation) predict(model, new_observation, type="response")$predictions[,"DS"]

Brok_SP="Zamenis lineatus"
Brok_TAXONID<-species$taxonid[species$scientific_name==Brok_SP]

predict.function(M0_pred[[IND]], DD[DD$taxonid==Brok_TAXONID,names(DD) %in% VARS])
DD$Prob.DS[DD$taxonid==Brok_TAXONID]

Res_Broken<- breakDown::broken(M0_pred[[IND]], DD[DD$taxonid==Brok_TAXONID,names(DD) %in% VARS], data = DD,
                    predict.function = predict.function, 
                    direction = "down")

Res_BrokenDF<-data.frame(Var=Res_Broken$variable, Contribution=Res_Broken$contribution) %>% .[order(.$Contribution, decreasing=T),] %>% subset(., .$Var %not in% c("(Intercept)", "final_prognosis"))

ggplot(Res_BrokenDF)+
  geom_bar(aes(x=Var, y=Contribution), stat="identity")+
  theme_minimal() %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 0.5))





### GBIF map species (for details in Fig. 4)
library(rgbif)

NAME="Zamenis lineatus"

distSP<-distributions[distributions$binomial==NAME,] %>% st_transform(., "+init=epsg:4326")

# Download and remove records with no spatial coordinates
GBIF <- rgbif::occ_search(scientificName = NAME, hasCoordinate = T, limit=10000)$data 
GBIF <- GBIF %>% filter(!is.na(decimalLongitude)) %>% filter(!is.na(decimalLatitude))

# Year column
GBIF$New_data<- GBIF$year > as.numeric(species$published_year[species$scientific_name==NAME])
GBIF<-subset(GBIF, is.na(New_data)==F)

countries<-read_sf("0.Data/Administrative boundaries/Admin_dissolved_by.country.shp") %>% st_transform(., "+init=epsg:4326")


### Extent GBIF (choose LIM depending on what is more widespread between the range and the GBIF data)
gbif_shp<-SpatialPoints(cbind(GBIF$decimalLongitude, GBIF$decimalLatitude), proj4string=CRS("+init=epsg:4326")) %>% st_as_sf(.)
LIM<-1.2*raster::extent(st_union(distSP, gbif_shp)) ; LIM@xmin<-0.95*LIM@xmin ; LIM@xmax<-1.05*LIM@xmax

# Plot
ggplot()+
  geom_sf(data=st_crop(countries, LIM), fill="gray85", col=NA)+
  geom_sf(data=distSP, fill="#ffffbf", col="gray50")+
  geom_point(data=GBIF, aes(decimalLongitude, decimalLatitude, col=New_data), show.legend=F, size=4)+
  scale_colour_manual(values=c("#91bfdb", "#d73027"))+
  theme_void()



### WOS species (for details in Fig. 4)
library(httr) ; library(jsonlite)

QUERY_WOS <- paste0('TS=("', NAME, '")') # Important to have real quote marks (" is ok ' is not)
  
resp1 <- httr::GET('https://api.clarivate.com/api/woslite', 
                     httr::add_headers(accept = 'application/json', `X-APIKey` = KEY_WOS), 
                     query = list(databaseId = 'WOK', 
                                  usrQuery = QUERY_WOS, 
                                  count = 100,  # Can be up to 100 for the list of papers
                                  firstRecord = 1, # To check when I extract the list of papers on DD, but I think 1 is the most recent
                                  publishTimeSpan=paste("1500-01-01", Sys.Date(), sep="+")))
  
data1 = fromJSON(rawToChar(resp1$content))
data1$Data



