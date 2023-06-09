library(pdp)

Pretty_name<-function(x){
  revalue(x, c(
    "res.cap"="Research capacity",
    "nb_GBIFgeoMOD"="# GBIF records",
    "time_descMOD"="Time since description",
    "N.hab1"="# habitats",
    "med_GDP"="Median GDP",
    "res.cap"= "Research capacity",
    "corruption"="Corruption",
    "range_magnitude"="Range size magnitude", 
    "realmsTERR_mode"="Realms",
    "Sampling_MED"="GBIF sampling effort", 
    "DDness"="Overlap with DD species",
    "acc_terrMED"="Travel time to cities (median)", 
    "acc_terrQ10"="Travel time to cities (Q10)",
    "taxo"="Taxonomical order",
    "alt_med"="Altitude",
    "hpop"="Human population", "hrur"="Rural population",
    "realmsMAR_mode"="Realms", "realmsMAR_var"="# realms",
    "acc_marMED"="Distance to ports (median)", "acc_marQ10"="Distance to ports (Q10)",
    "hf_mar"="Human footprint",
    "depth_med"="Depth",
    "RLA"="Red List Authority",
    "trait_availability_mammals"="Trait availability", 
    "gbifSP_coverageMOD"="Coverage of GBIF records", 
    "trait_availability_amphibians"="Trait availability",
    "WOSMOD"="# articles (WOS)",
    "conflictMEDYEAR"="Armed conflicts",
    "bodymass_mammals"="Body mass",
    "SVL_amphibians"="Body SVL",
    "roadsQ90"="Road density",
    "rur.prop"="Rural proprotion",
    "fishingMED"="Fishing intensity",
    "Depth_min"="Depth",
    "zoos"="Zoos and Aquaria",
    "nocturnal_mammals"="Nocturnality"
    
  ), warn_missing = F)
}





### Function to calculate the partial dependence of all plots, I have to do all at once to be able to have the same yaxis on 4 facets
ParDepCalc<-function(mod, df){
  
  Var1<-as.character(Imp.Vars$Variable[1])
  Var2<-as.character(Imp.Vars$Variable[2])
  Var3<-as.character(Imp.Vars$Variable[3])
  Var4<-as.character(Imp.Vars$Variable[4])
  
  # Create dataframe (removing above 95% for the plot)
  par.dep1<-partial(mod, pred.data = df, pred.var=Var1, prob=T, quantiles=T, probs=seq(0,0.95,0.01)) ; names(par.dep1)<-c("x", "y") ; par.dep1$Var<-"V1"
  par.dep2<-partial(mod, pred.data = df, pred.var=Var2, prob=T, quantiles=T, probs=seq(0,0.95,0.01)) ; names(par.dep2)<-c("x", "y") ; par.dep2$Var<-"V2"
  par.dep3<-partial(mod, pred.data = df, pred.var=Var3, prob=T, quantiles=T, probs=seq(0,0.95,0.01)) ; names(par.dep3)<-c("x", "y") ; par.dep3$Var<-"V3"
  par.dep4<-partial(mod, pred.data = df, pred.var=Var4, prob=T, quantiles=T, probs=seq(0,0.95,0.01)) ; names(par.dep4)<-c("x", "y") ; par.dep4$Var<-"V4"
  par.dep<-rbind(par.dep1, par.dep2, par.dep3, par.dep4)
  
  return(par.dep)
}


### Plot the partial dependence (four individual plots to be more customable)
perfectPartialPlot <- function(Var2plotNUM, par.dep, df){
  
  Var2plot<-as.character(Imp.Vars$Variable[Var2plotNUM])
  
  # Use partial dependence calculated in ParDepCalc
  par.depSUB<-subset(par.dep, par.dep$Var==paste0("V", Var2plotNUM))
  
  # Choose xlab names
  xlab_manu<-Pretty_name(Var2plot)
  if(nchar(xlab_manu)<16){xlab_manu<-paste(xlab_manu, " XXXXXXXXXX")} # To ensure all xlab are on 2 lines
  
  
  # Sqrt if the max is above five times the median
  TRANS<-"identity" ; if(max(par.depSUB$x) > 5*par.depSUB$x[round(nrow(par.depSUB)/2)]){TRANS<-"sqrt"}
  
  # Plot
  BREAKS=seq(par.depSUB[1,1],par.depSUB[nrow(par.depSUB),1], length.out=30)
  if(TRANS=="sqrt"){BREAKS<-(seq(sqrt(par.depSUB[1,1]), sqrt(par.depSUB[nrow(par.depSUB),1])^2, length.out=30))}
  
  Hist<-ggplot(df)+
    geom_histogram(aes(x=df[,Var2plot]), breaks=BREAKS)+
    scale_x_continuous(limits=c(par.depSUB[1,1],par.depSUB[nrow(par.depSUB),1]), expand=c(0,0), trans=TRANS)+
    scale_y_continuous(expand=c(0,0))+
    theme_void()
  
  YMIN<-min(par.dep$y) - 0.15*(max(par.dep$y)-min(par.dep$y))
  YMAX<-max(par.dep$y) + 0.05*(max(par.dep$y)-min(par.dep$y))
  
  if(Var2plotNUM>1){YLAB_lab<-""} else {YLAB_lab<-"Probability of being data sufficient"}
  if(TRANS=="sqrt"){HistXmin=sqrt(par.depSUB[1,1]); HistXmax=sqrt(par.depSUB[nrow(par.depSUB),1])} else {HistXmin=par.depSUB[1,1]; HistXmax=par.depSUB[nrow(par.depSUB),1]}
  
  
  GG<-ggplot(par.depSUB)+
    geom_line(aes(x=x, y=y), col="#33a02c", size=2)+
    ylab(YLAB_lab)+
    xlab(str_wrap(xlab_manu, width=15))+ # Sets the number of characters that can appear on one line
    scale_y_continuous(limits=c(YMIN, YMAX), expand=c(0,0))+
    scale_x_continuous(trans=TRANS)+ 
    annotation_custom(ggplotGrob(Hist), xmin=HistXmin, xmax=HistXmax, ymin=YMIN, ymax=min(par.depSUB$y))
  
  if(Var2plotNUM==1){return(GG+theme_minimal() %+replace% theme(plot.margin=unit(rep(0,4), units="points")))} else{
    return(GG+theme_minimal() %+replace% theme(plot.margin=unit(rep(0,4), units="points"), axis.text.y=element_text(colour="white"), axis.ticks.y=element_line(colour="white")))
  }
  
}







### Calculate performance

DD_perf<-function(THRE, perf.df){
  
  # Calculate the threshold that gives equal sensitivity and specificity or highest TSS (otherwise assign the quantitative value for 0.5 and 0.65)
  if(THRE %in% c("EQUAL", "MAXTSS", "SENS90")){
    test.thre<-data.frame(Threshold=seq(0,1,0.01), Sensitivity=NA, Specificity=NA, TSS=NA)
    
    for(i in 1:nrow(test.thre)){
      comp.thresh<-table(species_MOD$DD, factor(species_MOD$Predict_validation<test.thre$Threshold[i], c('FALSE', 'TRUE')))
      test.thre$Specificity[i]<-(comp.thresh["DD", "TRUE"]/sum(comp.thresh["DD",]))
      test.thre$Sensitivity[i]<-(comp.thresh["DS", "FALSE"]/sum(comp.thresh["DS",]))
    }
    
    test.thre$TSS<-test.thre$Sensitivity + test.thre$Specificity -1
    test.thre$Diff.SensSpe<-abs(test.thre$Sensitivity-test.thre$Specificity)
    
    if(THRE=="MAXTSS"){THRE.QTT<-test.thre$Threshold[which(test.thre$TSS==max(test.thre$TSS))][1]} # Selects the threshold that provides the highest TSS (the lower value in case there are several with the exact same TSS)
    
    if(THRE=="EQUAL"){THRE.QTT<-test.thre$Threshold[which(test.thre$Diff.SensSpe==min(test.thre$Diff.SensSpe))][1]} # Selects the threshold that provides the closest sensitivity and specificity (the lower value in case there are several with the exact same difference)
    
    if(THRE=="SENS90"){test.thre<-subset(test.thre, test.thre$Sensitivity>0.9) ; THRE.QTT<-test.thre$Threshold[which(test.thre$TSS==max(test.thre$TSS))][1]} # # Selects the threshold that provides the highest TSS but with specificity staying above 0.9 (the lower value in case there are several with the exact same TSS)
    
  } else {THRE.QTT<-THRE}
  

  # Calculate performance metrics
  comp<-table(species_MOD$DD, factor(species_MOD$Predict_validation<THRE.QTT, c('FALSE', 'TRUE')))
  
  # Save metrics
  perf.df[perf.df$Metric=="Accuracy" & perf.df$Group==Group & perf.df$Domain==Domain & perf.df$Threshold==THRE, c("Value", "Time")]<-c( (comp["DS","FALSE"]+comp["DD", "TRUE"])/sum(comp) , as.character(Sys.time()))
  perf.df[perf.df$Metric=="Sensitivity" & perf.df$Group==Group & perf.df$Domain==Domain & perf.df$Threshold==THRE, c("Value", "Time")]<-c( (comp["DS", "FALSE"]/sum(comp["DS",])), as.character(Sys.time()))
  perf.df[perf.df$Metric=="Specificity" & perf.df$Group==Group & perf.df$Domain==Domain & perf.df$Threshold==THRE, c("Value", "Time")]<-c( (comp["DD", "TRUE"]/sum(comp["DD",])), as.character(Sys.time()))
  perf.df[perf.df$Metric=="TSS" & perf.df$Group==Group & perf.df$Domain==Domain & perf.df$Threshold==THRE, c("Value", "Time")]<-c(((comp["DS", "FALSE"]/sum(comp["DS",])) + (comp["DD", "TRUE"]/sum(comp["DD",])) -1), as.character(Sys.time()))
 
  # Save used threshold
  perf.df$Threshold.QTT[perf.df$Group==Group & perf.df$Domain==Domain & perf.df$Threshold==THRE]<-THRE.QTT
  
  
  return(perf.df) 
}

