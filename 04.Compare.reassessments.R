

### Charge results
DD_tot<-read.csv("3.Figures/0.Supp/Extended_data_S1.csv")
DD_tot<-subset(DD_tot, DD_tot$DD==T)





### Download current version of the Red List
search<-as.data.frame(rl_search("Lophornis brachylophus", key=Vapi)$result) ; search[1,]<-NA ; search$Group<-NA

for(S in 1:nrow(DD_tot)){

    # Search TO ADD
    tryCatch({Ssearch<-as.data.frame(rl_search(DD_tot$scientific_name[S], key=Vapi)$result)} ,error=function(e){ # I have to do that because the connection bugs sometimes, so I re-run the function after 1/2 sec
      cat(paste0("\n Bug Search with ", DD_tot$scientific_name[S]))
      Sys.sleep(15)
      cat(" slept")
      Ssearch<-as.data.frame(rl_search(DD_tot$scientific_name[S], key=Vapi)$result)
      cat(paste0(" Solved: N=", nrow(Ssearch), "\n"))})

    if(nrow(Ssearch)>0){
      search<-rbind.fill(search, Ssearch)} else {cat("No row in Ssearch for ", DD_tot$scientific_name[S], "\n")}

    if((S/100)==round(S/100)){cat(paste0("S=",S, "\n")) ; saveRDS(search, "1.Tables/Search_DD_2022.1TEMP.rds")}

  }

saveRDS(search, "1.Tables/Search_DD_2022.1.rds")
#search<-readRDS("1.Tables/Search_DD_2022.1.rds")
search<-readRDS("C:/Users/Victor/Documents/sRedList/Platform/InProgress/sredlist-server-develop/Species/species-all-page.rds")

### Compare
DD_tot$New_Cat<-search$category[match(DD_tot$scientific_name, search$scientific_name)]
DD_tot$New_Year<-search$published_year[match(DD_tot$scientific_name, search$scientific_name)]

DD_reass<-subset(DD_tot, DD_tot$New_Year==2022)
DD_reass$New_Cat2<-ifelse(DD_reass$New_Cat=="DD", "DD", "DS")

table(DD_reass$New_Cat)
table(DD_reass$Group, DD_reass$New_Cat2)




### Plot main (priority score)
DD_reass$Vgroup<-factor(DD_reass$Group, c("Mammal", "Reptile", "Amphibian", "Fish", "Odonata"))

G3=ggplot(DD_reass)+
  geom_violin(aes(x=New_Cat2, y=PrioDS), fill="#c994c7ff")+
  geom_point(aes(x=New_Cat2, y=PrioDS))+
  facet_wrap(~Vgroup, nrow=1, drop=T)+
  xlab("Reassessment status")+
  ylab("Priority score")+
  ylim(c(0,1))+
  theme_minimal()

cowplot::save_plot("3.Figures/0.Main/Fig.5.validation.reassessment.svg", G3, base_height = 4, base_width = 12)

t.test(DD_reass$PrioDS ~ DD_reass$New_Cat2, alternative="less")


mean(DD_reass$PrioDS[DD_reass$New_Cat2=="DD"])
sd(DD_reass$PrioDS[DD_reass$New_Cat2=="DD"])
mean(DD_reass$PrioDS[DD_reass$New_Cat2=="DS"])
sd(DD_reass$PrioDS[DD_reass$New_Cat2=="DS"])

GR="Fish"
mean(DD_reass$PrioDS[DD_reass$New_Cat2=="DD" & DD_reass$Group==GR])
sd(DD_reass$PrioDS[DD_reass$New_Cat2=="DD" & DD_reass$Group==GR])
mean(DD_reass$PrioDS[DD_reass$New_Cat2=="DS" & DD_reass$Group==GR])
sd(DD_reass$PrioDS[DD_reass$New_Cat2=="DS" & DD_reass$Group==GR])
t.test(DD_reass$PrioDS[DD_reass$Group==GR] ~ DD_reass$New_Cat2[DD_reass$Group==GR], alternative="less")



### Plot Supplementary (pDS and iDS separately)
library(gridExtra)

Violin<-grid.arrange(
ggplot(DD_reass)+
  geom_violin(aes(x=New_Cat2, y=pDS), fill="#c994c7", draw_quantiles = 0.5)+
  geom_point(aes(x=New_Cat2, y=pDS))+
  xlab("Reassessment category")+
  ylab("pDS")+
  facet_wrap(~Vgroup, ncol=5)+
  theme_minimal(),

ggplot(DD_reass)+
  geom_violin(aes(x=New_Cat2, y=iDS), fill="#c994c7", draw_quantiles = 0.5)+
  geom_point(aes(x=New_Cat2, y=iDS))+
  xlab("Reassessment category")+
  ylab("ΔpDS")+
  facet_wrap(~Vgroup, ncol=5)+
  theme_minimal(),

nrow=2)

cowplot::save_plot("3.Figures/0.Supp/Violin_plot_NewRLVersion.SUPP.png", Violin, base_height=4.5, base_width = 10)







