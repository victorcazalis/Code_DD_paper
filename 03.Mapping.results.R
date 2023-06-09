

countries<-read_sf("0.Data/Administrative boundaries/Admin_dissolved_by_country_Simplif_0.005MOLLWEIDE.shp")
Map_plot<-list()


for(GR in 1:5){

Group=c("MAMMALS", "REPTILES", "AMPHIBIANS", "FISH", "FW_ODONATA")[GR]
Group2<-revalue(Group, c("AMPHIBIANS"="Amphibian", "REPTILES"="Reptile", "FISH"="Fish", "MAMMALS"="Mammal", "FW_ODONATA"="Odonata"))



### Recharge distributions (needed to recreate the grid with the exact same dimensions)
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

if(Group=="REPTILES"){
  distributions$tobuffer<-distributions$binomial %in% species$scientific_name[species$marine_system==T & species$terrestrial_system==F & species$freshwater_system==F]
  distributions[distributions$tobuffer==T,]<-st_buffer(distributions[distributions$tobuffer==T,], 5000)
}

gc()


### Calculate species of priority for reassessment
species<-readRDS("1.Tables/Species.merged.april22.rds")
res<-DD_tot[DD_tot$Vgroup==Group2,]
QUANT<-quantile(res$PrioDS, probs=c(0,0.75,1))
res$Dist_cat<-cut(res$PrioDS ,breaks=QUANT, labels=c("NotPrio", "Prio"))
res$scientific_name<-species$scientific_name[match(res$taxonid, species$taxonid)]


### Prepare grid
grid_raw<-readRDS(paste0("1.Tables/intersection_all_gridDDness_TEMPO.", Group, ".rds")) ; names(grid_raw)<-c("Species", "Cell_id")
grid_DD<-subset(grid_raw, grid_raw$Species %in% res$scientific_name)
grid_DD$Dist_cat<-res$Dist_cat[match(grid_DD$Species, res$scientific_name)]

DD_density<-ddply(grid_DD, .(Cell_id), function(x){data.frame(
  N_DD=nrow(x),
  N_DD25=nrow(x[x$Dist_cat=="Prio",])
)})

grid<-st_make_grid(distributions, cellsize=100000, square=F) %>% st_transform(., st_crs(distributions)) %>% st_as_sf(.)
grid$N_DD<-DD_density$N_DD[match(rownames(grid), DD_density$Cell_id)]
grid$N_DD25<-DD_density$N_DD25[match(rownames(grid), DD_density$Cell_id)]

### Plot
Plot_N_DD<-ggplot()+
  geom_sf(data=countries, fill="gray90", col=NA)+
  geom_sf(data=grid[is.na(grid$N_DD)==F,], aes(fill=N_DD), col=NA)+
  scale_fill_gradient(low="#fed976", high="#b10026", trans="log10", limits=c(1,max(grid$N_DD, na.rm=T)), name="")+
  ggtitle(ifelse(GR==1, "Number of DD species", ""))+
  theme_void() %+replace%   theme(plot.title=element_text(hjust=0.5, size=12))

Plot_N_DD25<-ggplot()+
  geom_sf(data=countries, fill="gray90", col=NA)+
  geom_sf(data=grid[is.na(grid$N_DD25)==F & grid$N_DD25>0,], aes(fill=N_DD25), col=NA)+
  scale_fill_gradient(low="#fed976", high="#b10026", trans="log10", limits=c(1,max(grid$N_DD, na.rm=T)), name="")+
  ggtitle(ifelse(GR==1, "Number of species to reassess", ""))+
  theme_void() %+replace%   theme(plot.title=element_text(hjust=0.5, size=12))

Map_plot[[GR]]<-grid.arrange(Plot_N_DD, Plot_N_DD25, ncol=2, top=grid::textGrob(Group2,gp=grid::gpar(fontsize=30,font=2)))

}



Map_tot<-grid.arrange(Map_plot[[1]], Map_plot[[2]], Map_plot[[3]], Map_plot[[4]], Map_plot[[5]], ncol=1)

cowplot::save_plot(paste0("3.Figures/0.Supp/Maps.AllGroups.png"), Map_tot, base_width=15, base_height=20)


