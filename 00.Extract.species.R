


### Build data set to fill
search<-as.data.frame(rl_search("Lophornis brachylophus", key=Vapi)$result) ; search[1,]<-NA ; search$Group<-NA
Groups<-rep(NA, 14)



for(G in 0:14){ ## CAREFUL ! THERE IS A PAGE 0!!
  
  Gtable<-data.frame(rl_sp(page=G, key=Vapi)$result)
  Gtable<-subset(Gtable, is.na(Gtable$infra_rank)==TRUE) 
  Groups[(G+1)]<-paste0("Page.", G)
  
  cat(paste0("\n Starting a new group: ", Groups[(G+1)], " with ", nrow(Gtable), " species \n"))
  
  for(S in 1:nrow(Gtable)){
   
    # Search TO ADD
    tryCatch({Ssearch<-as.data.frame(rl_search(Gtable$scientific_name[S], key=Vapi)$result)} ,error=function(e){ # I have to do that because the connection bugs sometimes, so I re-run the function after 1/2 sec
      cat(paste0("\n Bug Search with ", Gtable$scientific_name[S])) ; beep(4)
      Sys.sleep(15)
      cat(" slept")
      Ssearch<-as.data.frame(rl_search(Gtable$scientific_name[S], key=Vapi)$result)
      cat(paste0(" Solved: N=", nrow(Ssearch), "\n"))})
    
    if(nrow(Ssearch)>0){
    Ssearch$Group<-Groups[(G+1)]
    
    search<-rbind.fill(search, Ssearch)} else {cat("No row in Ssearch for ", Gtable$scientific_name[S], "\n")}
    
    if((S/100)==round(S/100)){cat(paste0("S=",S, "\n"))}
    
  }
  
  saveRDS(search, paste0("Species.download.", G, ".rds"))
  
  cat(paste0("Group completed: ", G)) ; beep(1)
}


# Not included: Tachigali pilgeriana (G1), Coussarea cepha?loides (G2), Tricholoma borgsjoe?nse (G9)
search2<- search %>% distinct(taxonid, .keep_all=T)

saveRDS(search2, paste0("1.Tables/Species.april22.merged.rds"))








