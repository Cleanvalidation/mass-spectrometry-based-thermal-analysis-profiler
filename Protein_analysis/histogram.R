#plot global histograms for variability 
his_sp<-function(Df1,df.temps,MD=FALSE){
  if(any(names(Df1)=="C")){
    Df1<-dplyr::bind_rows(Df1) 
    Df1<-Df1%>% dplyr::rename("temperature"="C")
    
  }
  if(any(names(dplyr::bind_rows(Df1))=="sample_name.x")){
    test<-dplyr::bind_rows(Df1) %>% dplyr::rename("sample_name"=ifelse(any(names(.)=="sample_name.x"),"sample_name.x","sample_name.y"))
    test<-test%>% dplyr::select(treatment,CV_pct,treatment,C,sample_name) %>% unique(.)
    test<-test %>% dplyr::rename("temperature"="C")
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
    
    
    test<-test %>% dplyr::mutate(carrier=ifelse(str_detect(test$sample_name,"NOcarrier"),"+ Carrier","- Carrier"),
                                 FAIMS = ifelse(str_detect(test$sample_name,"NO_FAIMS"),"+ FAIMS","- FAIMS"),
                                 Phi = ifelse(str_detect(test$sample_name,"PhiSDM"),"+ PhiSDM","- PhiSDM"))
  }else{
    test<-dplyr::bind_rows(Df1)
    test<-test %>% unique(.)
    if (any(names(test)=="C")){
      test<-test %>% dplyr::rename("temperature"="C")
    }
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
  }
  
  
  test<-test %>% dplyr::mutate(carrier=ifelse(str_detect(test$sample_name,"NOcarrier"),"+ Carrier","- Carrier"),
                               FAIMS = ifelse(str_detect(test$sample_name,"NO_FAIMS"),"+ FAIMS","- FAIMS"),
                               Phi = ifelse(str_detect(test$sample_name,"PhiSDM"),"+ PhiSDM","- PhiSDM"))
  
  if(isTRUE(MD)){
    test$treatment<-as.factor(test$treatment)
    
    ggplot(test,aes(y=CV_pct,x=treatment,fill=treatment))+
      facet_grid(c(~temperature,~carrier),scales="free_x")+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample_id")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,200)
  }else{  
    
    test$treatment<-test$treatment
    ggplot(test,aes(y=CV_pct,x=treatment,fill=treatment))+
      facet_grid(~temperature)+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample_id")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,200)
    
    
  }
}
#plot global histograms for variability 
Violin_panels<-function(df_raw,df.temps,MD=TRUE){
  if(any(names(df_raw)=="C")){
    Df1<-dplyr::bind_rows(df_raw) 
    Df1<-Df1%>% dplyr::rename("temperature"="C")
    
  }
  if(any(names(dplyr::bind_rows(Df1))=="sample_name.x")){
    test<-dplyr::bind_rows(Df1) %>% dplyr::rename("sample_name"=ifelse(any(names(.)=="sample_name.x"),"sample_name.x","sample_name.y"))
    test<-test%>% dplyr::select(treatment,CV_pct,treatment,C,sample_name) %>% unique(.)
    test<-test %>% dplyr::rename("temperature"="C")
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum.File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$treatment<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
    test<-test %>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                 "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                 "I_Interference"="Isolation_Interference_[%]")
    
    
  }else{
    test<-df_raw
    test<-test %>% dplyr::left_join(df.temps,by="temp_ref")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum.File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$treatment<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
    test<-test %>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                 "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                 "I_Interference"="Isolation_Interference_[%]")
    
    
    
  }
  
  
  if(isTRUE(MD)){
    test$sample_name<-as.factor(test$sample_name)
    
    P<-ggplot(test,aes(y=I_Interference,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("Isolation Interference (%)")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,100)+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=IonInjTime,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("Ion Injection Time (ms)")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,75)+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=PEP,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("PEP Score")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=DeltaM,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      scale_fill_brewer(palette = "Dark2") +
      ggplot2::ylab(expression(paste(Delta, " M (ppm)")))+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(-15,15)+
      guides(fill=guide_legend(title="Method"))
    
  }
}