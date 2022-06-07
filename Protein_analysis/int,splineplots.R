
Int_plot<-function(x,Peptide=FALSE,raw=FALSE){
  if(any(x$temp_ref=="131C")){
    x<-x %>% dplyr::filter(!temp_ref=="131C")
  }
  x<-data.frame(x)
  x$treatment<-as.factor(x$treatment)
  if(!any(stringr::str_detect(names(x),"uniqueID"))){
    x<-x %>% dplyr::rename("uniqueID"="Accession")
    x$uniqueID<-as.factor(x$uniqueID)
  }
  if(any(stringr::str_detect(x$sample_name,"S"))&any(stringr::str_detect(x$Spectrum.File,"FAIMS"))){
    x$sample_name<-stringr::str_replace(x$sample_name,"S","\u03A6")
  }
  if(any(names(x)=="I3")){
    x<-x %>% dplyr::rename("I"="I3")
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x<-x %>% dplyr::rename("I"="value")
  }
  x<-x %>% dplyr::filter(!is.na(x$I))
  if(!any(names(x)=="C")){
    x<-x %>% dplyr::rename("C"="temperature")
  }
  x$C<-as.factor(x$C)
  x$I<-as.numeric(x$I)
  if(!isTRUE(Peptide)){#if this is a protein file
    if(!isTRUE(raw)){
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Normalized intensity protein")+
        ggtitle(x$sample_name[1])+
        ylim(-0.1,5)+
        theme(legend.position="bottom")+ labs(colour = "Temperature (\u00B0C)")
    }else{
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,10000000)
    }
  }else{#if this is a peptide file
    if(!isTRUE(raw)){
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Normalized intensity")+
        ggtitle(x$sample_name[1])+
        ylim(-0.1,2)+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")
    }else if(max(x$I,na.rm=TRUE)>900000){#if this is a PSM file
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,10000000)
    }else{#if this is a PG file
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,2)
    }
  }
  
}
PlotTrilinear<-function(df_norm,target,df.temps,Ft,filt,Peptide=FALSE,show_results=FALSE){
  df_norm$CC<-ifelse(df_norm$treatment=="vehicle",0,1)
  if(isTRUE(Peptide) & any(names(df_norm)=="XCorr")){
    #remove columns not needed for curve fitting 
    df_norm<-df_norm %>% dplyr::select(-XCorr,-temp_ref,-Modifications,-MissedCleavages,-DeltaM,-"Annotated_Sequence",-tidyr::contains("PEP"),-Charge)
    if(any(names(df_norm)=="I_Interference")){
      df_norm<-df_norm %>% dplyr::select(-I_Interference,-IonInjTime,-S_N,-Spectrum.File)
    }
    
  }
  if(any(names(df_norm)=="I3")){
    df_norm<-df_norm %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
  }
  ##SCRIPT STARTS HERE
  DF<-df_norm %>% dplyr::group_split(uniqueID) #split null treatment only by protein ID
  d_<-df_norm %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,treatment) #split vehicle treatment
  d_1<-df_norm %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,treatment) #split treated treatment
  if(length(d_1)==0){
    d_1<-d_
  }
  if(length(d_)==0){
    warning(paste0("No vehicle curves found for ",df_norm$sample_name[1]))
  }
  #convert to data frame for uniqueID presence
  DF<-dplyr::bind_rows(DF)
  d_<-dplyr::bind_rows(d_)
  d_1<-dplyr::bind_rows(d_1)#1 
  #make sure uniqueIDs are present for treated and vehicle
  CID<-intersect(d_1$uniqueID,d_$uniqueID)
  CID<-intersect(CID,DF$uniqueID)
  DF<-DF %>% subset(uniqueID %in% CID)
  DF$LineRegion<-1
  d_<-d_%>% subset(uniqueID %in% CID)
  d_$LineRegion<-1
  d_1<-d_1%>% subset(uniqueID %in% CID)
  d_1$LineRegion<-1
  #split treatment into equal-sized lists
  DF<-DF %>%  dplyr::group_split(uniqueID) 
  d_<-d_ %>% dplyr::group_split(uniqueID) 
  d_1<-d_1 %>% dplyr::group_split(uniqueID) 
  
  #preallocate list
  results<-vector(mode = "list", length(d_))
  results_t<-vector(mode = "list",length(d_1))
  results_n<-vector(mode = "list",length(DF))
  
  
  results<-suppressWarnings(DLR(d_))#First guess at line regions
  results_t<-suppressWarnings(DLR(d_1))
  results_n<-suppressWarnings(DLR(DF))
  
  
  
  #reassign shared points between line regions
  df_<-suppressWarnings(purrr::map2(results,d_,function(x,y) CP(x,y,PSM=Peptide)))
  df_1<-suppressWarnings(purrr::map2(results_t,d_1,function(x,y) CP(x,y,PSM=Peptide)))
  DFN<-suppressWarnings(purrr::map2(results_n,DF,function(x,y) CP(x,y,PSM=Peptide)))
  
  
  #remove results to save space 
  rm(results,results_t,results_n,d_,d_1,DF)#10
  df_<-dplyr::bind_rows(df_) %>% dplyr::group_split(uniqueID)
  df_1<-dplyr::bind_rows(df_1)%>% dplyr::group_split(uniqueID)
  DFN<-dplyr::bind_rows(DFN)%>% dplyr::group_split(uniqueID)
  
  df_<-lapply(df_,function(x)x[order(x$C),])
  df_1<-lapply(df_1,function(x)x[order(x$C),])
  DFN<-lapply(DFN,function(x)x[order(x$C),])
  #get # identifier
  df_<-purrr::map2(df_,seq(df_),function(x,y)x %>% dplyr::mutate(N=y))
  df_1<-purrr::map2(df_1,seq(df_1),function(x,y)x %>% dplyr::mutate(N=y))
  DFN<-purrr::map2(DFN,seq(DFN),function(x,y)x %>% dplyr::mutate(N=y))
  
  #prealloate variables
  tlresults<-list()
  tlresults_PI<-list()
  #split data 
  df_<-dplyr::bind_rows(df_)
  df_1<-dplyr::bind_rows(df_1)
  DFN<-dplyr::bind_rows(DFN)
  #confidence intervals
  
  tlresults<-tlstat(DFN,df_,df_1,norm=FALSE,Filters=filt,Ftest=Ft,show_results=show_results)
  if(isTRUE(show_results)){
    
    return(tlresults)
  }
  #return filtered lists
  res<-tlf(tlresults,DFN,APfilt=FALSE,PF=FALSE)
  
  i=which(res[[1]]$uniqueID %in% target)
  plotTL1<-tlCI(i,res[[1]],res[[2]],res[[3]],overlay=TRUE,residuals=FALSE,df.temps=df.temps,PSMs=Peptide)
  
  return(plotTL1)
}


###############################
plot_Splines<-function(x,Protein="Q02750",df.temps,Filters=FALSE,fT=TRUE,show_results=TRUE,Peptide=TRUE,CI=TRUE,simulations=FALSE,CARRIER=FALSE,Frac=TRUE,raw=FALSE){
  Filters=Filters
  fT=fT
  if(any(names(x)=="Charge")){
    x$Charge<-as.factor(x$Charge)
  }
  if(any(names(x)=="value")&!any(names(x)=="C")){
    x<-x %>% dplyr::rename("C"="temperature")
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x<-x %>% dplyr::rename("I"="value")
    
  }
  if(any(names(x)=="Fraction")){
    x<-x %>% dplyr::select(-Fraction) %>% distinct(.)
    x<-x %>% dplyr::mutate(Fraction=1)
    
  }
  if(any(names(x)=="Spectrum.File")){
    x<-x %>% dplyr::select(-Spectrum.File) %>% distinct(.)
    
  }
  if(any(names(x)=="Accession")&!any(names(x)=="uniqueID")){
    x<-x %>% dplyr::rename("uniqueID"="Accession")
    
  }
  if(!any(stringr::str_detect(names(x),"replicate"))&any(names(x)=="Charge")){
    x$replicate<-x$Charge
  }
  if(length(unique(x$treatment))==1){
    return(warning("Please check that you have vehicle and treated samples in your data.  Only one of the conditions read."))
  }
  if(!any(names(x)=="replicate")&any(names(x)=="Fraction")){
    x<-x %>% dplyr::mutate(replicate=Fraction)
  }
  if(isTRUE(Peptide)){
    x<-x %>% dplyr::group_split(uniqueID,Annotated_Sequence,treatment,temp_ref) 
    x<-purrr::map(x,function(y) y%>% distinct(.) %>% dplyr::mutate(replicate=row.names(.)))
    x<-dplyr::bind_rows(x)
  }
  # if(any(stringr::str_detect(names(x),"replicate") & isFALSE(Peptide))){
  #   if(x$replicate[1]==""){
  #     x<-x %>% dplyr::mutate(replicate=stringr::str_sub(x$id, -1))
  #     
  #   }
  # }
  
  if(isTRUE(CARRIER)){
    df.temps<-df.temps %>% dplyr::filter(temperature<68)
  }
  
  if(!isTRUE(Peptide)){#if this is a protein file
    DFN<-x 
    df_<-x %>% dplyr::filter(treatment=="vehicle")
    df_1<-x %>% dplyr::filter(treatment=="treated")
    
    if(nrow(df_1)==0){
      df_1<-df_
    }
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    if(isTRUE(show_results)){
      spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
      return(spresults)
    }else{
      spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
      
    }
    if(any(class(spresults)=="list")){
      spresults<-dplyr::bind_rows(spresults)
      spresults<-spresults[!is.na(spresults$C),]
      res_sp<-spf(spresults,DFN,filters=Filters)
    }else{
      res_sp<-spf(spresults,DFN,filters=FALSE)
    }
    
    if(isTRUE(simulations)){
      ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
      return(ROCs)
    }
    
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
    
    return(Pred1)
  }else{#if this is a peptide file
    DFN<-dplyr::bind_rows(x)
    df_<-dplyr::bind_rows(x) %>% dplyr::filter(treatment=="vehicle")
    df_1<-dplyr::bind_rows(x) %>% dplyr::filter(treatment=="treated")
    
    if(nrow(df_1)==0){#if there are no treated files, make one up to avoid errors in curve fitting
      df_1<-df_
    }
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    spresults<-tryCatch(spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters),
                        error = function(cond){
                          message("Error found in fitting and summarise function")
                          message(cond)
                        })
    
    if(isTRUE(show_results)){
      return(spresults)
    }
    if(any(class(spresults)=="list")){
      spresults<-dplyr::bind_rows(spresults)
      spresults<-spresults[!is.na(spresults$C),]
      res_sp<-spf(spresults,DFN,filters=Filters)
    }else{
      
      DFN$temperature<-DFN$C
      res_sp<-spf(spresults,DFN,filters=FALSE)
    }
    
    if(isTRUE(simulations)){
      ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
      return(ROCs)
    }
    
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
    
    return(Pred1)
  }
  
  DFN<-x %>% dplyr::filter(uniqueID %in% as.character(Protein))
  df_<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),treatment=="vehicle")
  df_1<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),treatment=="treated")
  #get spline results
  spresults<-list()
  spresults_PI<-list()
  
  spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
  if(isTRUE(show_results)){
    return(spresults)
  }
  
  if(any(class(spresults)=="list")){
    spresults<-dplyr::bind_rows(spresults)
    spresults<-spresults[!is.na(spresults$C),]
    res_sp<-spf(spresults,DFN,filters=Filters)
  }else{
    DFN$temperature<-DFN$C
    res_sp<-spf(spresults,DFN,filters=Filters)
  }
  if(isTRUE(simulations)){
    ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
    return(ROCs)
  }
  #saveIDs filtered
  i<-which(res_sp[[1]]$uniqueID %in% Protein)
  #generate 95%CI for splines
  Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
  return(Pred1)
  
  
}


plot_Sigmoidal<-function(df_,Protein,Peptide=FALSE){
  #remove duplicated intensity values
  df_<-df_ %>% distinct(.)
  df_$sample_name<-str_replace(df_$sample_name,"S","\u03A6")
  
  PlSig<-try(sigC(df_,Protein,Peptide=Peptide))
  
  ID<-unique(PlSig$uniqueID)
  
  sig<- try(sigfit(PlSig,Peptide=Peptide))
  return(sig)
}
