###Upset plots#################
#To generate information on missing values
################################
upPSM_SN<-function(df_){
  
  df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                            "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                            "I_Interference"="Isolation_Interference_[%]")
  
  
  #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
  df_1<-dplyr::bind_rows(df_)
  df_1$uniqueID<-as.factor(df_1$uniqueID)
  df_1$treatment<-as.factor(df_1$treatment)
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1<-df_1%>% 
    dplyr::group_split(uniqueID,sample_id,Annotated_Sequence)
  df_1<-purrr::map(df_1,function(x) x %>% dplyr::mutate(missing_pct=(100*sum(is.na(x$I))/length(x$I))) %>% head(.,1)) 
  df_1<-dplyr::bind_rows(df_1)
  
  df_1$SNrank<-dplyr::ntile(df_1$S_N,3)
  minSNrank<-min(df_1[df_1$SNrank==3,'S_N'],na.rm=TRUE)
  maxSNrankL<-max(df_1[df_1$SNrank==1,'S_N'],na.rm=TRUE)
  low<-paste0("low"," < ",as.character(maxSNrankL))
  medium<-"medium"
  high<-paste0("high"," > ",as.character(minSNrank))
  df_1<-df_1 %>% dplyr::mutate(rank=ifelse(df_1$rank==1,"low",ifelse(df_1$rank==2,"medium","high")),
                               SNrank=ifelse(df_1$SNrank==3,"high",ifelse(df_1$SNrank==2,"medium","low")))
  
  df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
  df_1$SNrank<-factor(df_1$SNrank,levels=c("high","medium","low"))
  df_1$missing_pct<-round(df_1$missing_pct,1)
  df_1$DeltaM<-round(df_1$DeltaM,1)
  df_1<-df_1%>% 
    dplyr::group_split(sample_id)
  df_1<-purrr::map(df_1,function(x)
    pivot_wider(x,names_from=c(Charge),values_from=Charge,values_fill=NA) %>%
      dplyr::select(-uniqueID,-Spectrum.File,-Annotated_Sequence,-IonInjTime,-MissedCleavages,-S_N,-sample_id,-treatment,-C,-CC,-DeltaM,-PEP,-Modifications,-I_Interference,-missing_pct,-XCorr,-I,-Protein_value))
  
  df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(SNrank,sample_name))
  df_1<-purrr::map(df_1,function(x) 1+x[,sapply(x,class)=="numeric"])
  df<-lapply(df_1,function(x) colnames(x))
  df<-lapply(df,function(x) str_replace(x,x,
                                        paste0("Charge: ","+",x, sep = " ")))
  
  
  
  
  df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
  df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
  df_1<-purrr::map(df_1,function(x)x[!is.na(x$SNrank),])
  
  
  rating_scale = scale_fill_manual(name="Ranked S/N",
                                   values=c("low" ='#fee6ce', "medium" ='#fdae6b', "high"  = '#e6550d'))
  
  show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
  UPD<-purrr::map(df_1,function(x) upset_data(x,names(x),
                                              min_degree=1,
                                              n_intersections=N,
                                              with_sizes %>% dplyr::select(order, intersection) %>%
                                                distinct() %>% arrange(order) %>% pull("intersection"))
                  
  )
  check<-list()
  check<- purrr::map(UPD,function(x)upset(x,names(x),
                                          min_degree=1,
                                          sort_intersections="ascending",
                                          set_sizes=FALSE,
                                          guides='collect',
                                          n_intersections=N,
                                          height_ratio = 0.7,
                                          stripes='white',
                                          base_annotations=list(
                                            '# of PSMs'=intersection_size(
                                              counts=TRUE,
                                            )
                                          ),
                                          annotations =list(
                                            "Bin %"=list(
                                              aes=aes(x=intersection, fill=x$SNrank),
                                              
                                              geom=list(
                                                geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                
                                                geom_text(
                                                  aes(
                                                    label=!!aes_percentage(relative_to='intersection'),
                                                    color=ifelse(!is.na(x$SNrank), 'show', 'hide')
                                                  ),
                                                  stat='count',
                                                  position=position_fill(vjust = .5),
                                                  
                                                ),
                                                
                                                scale_y_continuous(labels=scales::percent_format()),
                                                show_hide_scale,
                                                rating_scale
                                                
                                              )
                                            )
                                          )
                                          
  )+ggtitle(str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
  check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                set_sizes=FALSE,
                guides='collect',
                n_intersections=N,
                
                stripes='white',
                base_annotations=list(
                  '# of PSMs'=intersection_size(
                    counts=TRUE,
                  )
                ),
                annotations =list(
                  "Ranked Intensity  "=list(
                    aes=aes(x=intersection, fill=df_1[[1]]$SNrank),
                    
                    geom=list(
                      geom_bar(stat='count', position='fill', na.rm=TRUE),
                      themes=theme(legend.position="bottom", legend.box = "horizontal"),
                      geom_text(
                        aes(
                          label=!!aes_percentage(relative_to='intersection'),
                          color=ifelse(!is.na(df_1[[1]]$SNrank), 'show', 'hide')
                        ),
                        stat='count',
                        position=position_fill(vjust = .5),
                        
                      ),
                      
                      scale_y_continuous(labels=scales::percent_format()),
                      show_hide_scale,
                      rating_scale
                      
                    )
                  )
                )
                
  )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
  y<-get_legend(check1$patches$plots[[1]])
  data<-unlist(lapply(check,function(x) x$labels$title))
  check<-check[order(data)]
  P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
  
  return(print(P))
  
}
#listUP<-upPSM_SN(df_raw)
#df_ is the input from df_norm which includes missing percentages
#returns segmented data for vehicle and treated
upMV <- function(df_,condition,plot_multiple=FALSE,PSMs=FALSE,Frac=TRUE){
  
  if(isTRUE(plot_multiple)){
    df_<-dplyr::bind_rows(df_)
    if(isTRUE(PSMs)){
      #df_<- df_ %>% dplyr::right_join(df.samples,by="sample_id")
      df_ <-df_ %>% dplyr::mutate(CC=ifelse(stringr::str_detect(treatment,"vehicle")==TRUE,0,1))#concentration values are defined in uM
      
      if(any(stringr::str_detect(names(df_),"Accession"))){
        df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                  "IonInjTime"="Ion_Inject_Time_[ms]","I_Interference"="Isolation_Interference_[%]")
      }else{
        df_<-df_%>% dplyr::rename("I"="value","C"="temperature","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                  "IonInjTime"="Ion_Inject_Time_[ms]","I_Interference"="Isolation_Interference_[%]")
      }
      
      #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
      df_1<-dplyr::bind_rows(df_)
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$treatment<-as.factor(df_1$treatment)
      df_1$sample_name<-as.factor(df_1$sample_name)
      if(isTRUE(Frac)){
        df_1<-dplyr::bind_rows(df_1)%>% 
          dplyr::group_split(uniqueID,sample_id,Annotated_Sequence,Fraction,treatment,sample_name)
      }else{
        df_1<-dplyr::bind_rows(df_1)%>% 
          dplyr::group_split(uniqueID,sample_id,Annotated_Sequence,treatment,sample_name)
      }
      
      df_1<-dplyr::bind_rows(df_1)
      if(any(is.na(df_1$rank))){
        if(any(stringr::str_detect(names(df_1),"Fraction"))){
          rank<-df_1 %>% dplyr::select(-rank) %>% dplyr::filter(C==min(.$C,na.rm=TRUE)) %>% dplyr::mutate(rank=dplyr::ntile(dplyr::desc(.$I),3),
                                                                                                          SNrank=dplyr::ntile(dplyr::desc(.$S_N),3)) %>% 
            dplyr::select("uniqueID","treatment","sample_name","sample_id","rank","Fraction","Annotated_Sequence","missing_pct","SNrank","Charge","S_N",-temp_ref)
          maxSNrank<-min(rank[rank$SNrank==1,'S_N'],na.rm=TRUE)
          minSNrank<-max(rank[rank$SNrank==3,'S_N'],na.rm=TRUE)
          low<-paste0("low"," < ",as.character(minSNrank))
          medium<-"medium"
          high<-paste0("high"," > ",as.character(maxSNrank))
          rank<-rank %>% dplyr::mutate(rank=ifelse(rank==1,"high",ifelse(rank==2,"medium","low")),
                                       SNrank=ifelse(SNrank==3,"low",ifelse(SNrank==2,"medium","high")))
          rank$SNrank<-factor(rank$SNrank)
          rank$rank<-factor(rank$rank)
          df_1<-dplyr::bind_rows(rank)
          
        }else{
          rank<-df_1 %>% dplyr::select(-rank) %>% dplyr::filter(C==min(.$C,na.rm=TRUE)) %>% dplyr::mutate(rank=dplyr::ntile(dplyr::desc(.$I),3),
                                                                                                          SNrank=dplyr::ntile(dplyr::desc(.$S_N),3)) %>% 
            dplyr::select("uniqueID","treatment","sample_name","sample_id","rank","Annotated_Sequence","missing_pct","SNrank","Charge","S_N",-temp_ref)
          maxSNrank<-min(rank[rank$SNrank==1,'S_N'],na.rm=TRUE)
          minSNrank<-max(rank[rank$SNrank==3,'S_N'],na.rm=TRUE)
          low<-paste0("low"," < ",as.character(minSNrank))
          medium<-"medium"
          high<-paste0("high"," > ",as.character(maxSNrank))
          rank<-rank %>% dplyr::mutate(rank=ifelse(rank==1,"high",ifelse(rank==2,"medium","low")),
                                       SNrank=ifelse(SNrank==3,"low",ifelse(SNrank==2,"medium","high")))
          rank$SNrank<-factor(rank$SNrank)
          rank$rank<-factor(rank$rank)
          df_1<-dplyr::bind_rows(rank)
        }
        #name<-dplyr::intersect(names(df_1),names(rank))
        #df_1<-df_1 %>% dplyr::right_join(rank,by=name)
        
      }
      #df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
      df_1$missing_pct<-100-round(df_1$missing_pct,1)
      df_1<-dplyr::bind_rows(df_1)%>% 
        dplyr::group_split(sample_name)
      #separate data by charge states
      Df_1<-purrr::map(df_1,function(x)x %>% dplyr::select(uniqueID,Annotated_Sequence,treatment,sample_id,Charge,rank,sample_name,SNrank) %>% distinct(.))
      
      hi<-dplyr::bind_rows(Df_1) %>% pivot_wider(names_from="Charge",values_from="Charge")
      #df_1<-purrr::map(hi,function(x)x %>% dplyr::rename_if(is.numeric))
      df_1<-hi %>% rename_at(names(.)[which(sapply(hi,class)=="numeric")],~paste0("Charge: ", .))
      df_1<-dplyr::bind_rows(df_1) %>% dplyr::group_split(sample_name)
      if(any(names(df_1[[1]])=="SNrank")){
        df_1<-purrr::map(df_1,function(x)x[!is.na(x$SNrank),])
      }
      
      rating_scale = scale_fill_manual(name="Ranked S/N",
                                       values=c("high" = '#e6550d', "medium" ='#fdae6b', "low"  = '#fee6ce'))
      
      show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
      
      check<-list()
      check<- purrr::map(df_1,function(x)upset(x,names(x)[which(stringr::str_detect(names(x),"Charge"))],
                                               min_degree=1,
                                               set_sizes=FALSE,
                                               guides='collect',
                                               height_ratio = 0.7,
                                               stripes='white',
                                               base_annotations=list(
                                                 '# of PSMs'=intersection_size(
                                                   counts=TRUE,
                                                 )
                                               ),
                                               annotations =list(
                                                 "Bin %"=list(
                                                   aes=aes(x=intersection, fill=x$SNrank),
                                                   
                                                   geom=list(
                                                     geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                     
                                                     geom_text(
                                                       aes(
                                                         label=!!aes_percentage(relative_to='intersection'),
                                                         color=ifelse(!is.na(x$SNrank), 'show', 'hide')
                                                       ),
                                                       stat='count',
                                                       position=position_fill(vjust = .5),
                                                       
                                                     ),
                                                     
                                                     scale_y_continuous(labels=scales::percent_format()),
                                                     show_hide_scale,
                                                     rating_scale
                                                     
                                                   )
                                                 )
                                               )
                                               
      )+ggtitle(stringr::str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
      check1<-upset(df_1[[1]],names(df_1[[1]])[which(stringr::str_detect(names(df_1[[1]]),"Charge"))],
                    min_degree=1,
                    set_sizes=FALSE,
                    guides='collect',
                    stripes='white',
                    base_annotations=list(
                      '# of PSMs'=intersection_size(
                        counts=TRUE,
                      )
                    ),
                    annotations =list(
                      "Ranked Intensity  "=list(
                        aes=aes(x=intersection, fill=df_1[[1]]$SNrank),
                        
                        geom=list(
                          geom_bar(stat='count', position='fill', na.rm=TRUE),
                          themes=theme(legend.position="bottom", legend.box = "horizontal"),
                          geom_text(
                            aes(
                              label=!!aes_percentage(relative_to='intersection'),
                              color=ifelse(!is.na(df_1[[1]]$SNrank), 'show', 'hide')
                            ),
                            stat='count',
                            position=position_fill(vjust = .5),
                            
                          ),
                          
                          scale_y_continuous(labels=scales::percent_format()),
                          show_hide_scale,
                          rating_scale
                          
                        )
                      )
                    )
                    
      )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
      y<-get_legend(check1$patches$plots[[1]])
      data<-unlist(lapply(check,function(x) x$labels$title))
      check<-check[order(data)]
      P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
      
      return(print(P))
      
    }else{#If this is a protein file
      df_1<-dplyr::bind_rows(df_)
      df_1<-df_1%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature")
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$treatment<-as.factor(df_1$treatment)
      df_1$sample_name<-as.factor(df_1$sample_name)
      df_1$rank<-ifelse(df_1$rank==3,"high",ifelse(df_1$rank==2,"medium","low"))
      df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
      
      df_1<-df_1%>% 
        dplyr::group_split(uniqueID,sample_id)
      df_1<-lapply(df_1,function(x) x[1,]) 
      df_1<-dplyr::bind_rows(df_1)
      
      df_1$missing_pct<-as.numeric(df_1$missing_pct)
      df_1$missing_pct<-round(df_1$missing_pct,0)
      
      
      df_1<-df_1 %>% dplyr::group_split(sample_name)
      df_1<-purrr::map(df_1,function(x)
        pivot_wider(x,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA) %>%
          dplyr::select(-C,-I,-missing,-treatment,-sample_id))
      
      df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(rank,sample_name))
      df_1<-purrr::map(df_1,function(x) 1+x[,lapply(x,class)=="numeric"])
      df<-lapply(df_1,function(x) colnames(x))
      df<-lapply(df,function(x) str_replace(x,x,
                                            paste0("missing ",x,"%", sep = " ")))
      
      
      
      
      df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
      df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
      df_1<-purrr::map(df_1,function(x)x[!is.na(x$rank),])
    }
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    check<-list()
    check<- purrr::map(df_1,function(x)upset(x,names(x)[!names(x)%in%c("rank","sample_name")],
                                             min_degree=1,
                                             set_sizes=FALSE,
                                             guides='collect',
                                             n_intersections=N,
                                             height_ratio = 0.7,
                                             stripes='white',
                                             base_annotations=list(
                                               '# of Replicates'=intersection_size(
                                                 counts=TRUE,
                                               )
                                             ),
                                             annotations =list(
                                               "Bin %"=list(
                                                 aes=aes(x=intersection, fill=x$rank),
                                                 
                                                 geom=list(
                                                   geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                   
                                                   geom_text(
                                                     aes(
                                                       label=!!aes_percentage(relative_to='intersection'),
                                                       color=ifelse(!is.na(rank), 'show', 'hide')
                                                     ),
                                                     stat='count',
                                                     position=position_fill(vjust = .5),
                                                     
                                                   ),
                                                   
                                                   scale_y_continuous(labels=scales::percent_format()),
                                                   show_hide_scale,
                                                   rating_scale
                                                   
                                                 )
                                               )
                                             )
                                             
    )+ggtitle(str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
    check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                  set_sizes=FALSE,
                  guides='collect',
                  n_intersections=N,
                  
                  stripes='white',
                  base_annotations=list(
                    '# of Replicates'=intersection_size(
                      counts=TRUE,
                    )
                  ),
                  annotations =list(
                    "Ranked Intensity  "=list(
                      aes=aes(x=intersection, fill=df_1[[1]]$rank),
                      
                      geom=list(
                        geom_bar(stat='count', position='fill', na.rm=TRUE),
                        themes=theme(legend.position="bottom", legend.box = "horizontal"),
                        geom_text(
                          aes(
                            label=!!aes_percentage(relative_to='intersection'),
                            color=ifelse(!is.na(rank), 'show', 'hide')
                          ),
                          stat='count',
                          position=position_fill(vjust = .5),
                          
                        ),
                        
                        scale_y_continuous(labels=scales::percent_format()),
                        show_hide_scale,
                        rating_scale
                        
                      )
                    )
                  )
                  
    )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
    y<-get_legend(check1$patches$plots[[1]])
    
    P<-ggarrange(plotlist=check,ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
    print(P)
  }else{
    ########Plot separately
    
    df_<-df_ %>% subset(sample_name== condition) 
    
    df_<-df_ %>% dplyr::rename("uniqueID"="Accession","C"="temperature","I"="value")
    
    df_$rank<-ifelse(df_$rank==1,"high",ifelse(df_$rank==2,"medium","low"))
    df_$rank<-factor(df_$rank,levels=c("high","medium","low"))
    df_<-df_ %>% dplyr::filter(sample_name==condition)%>%
      dplyr::group_split(uniqueID) 
    df_<-lapply(df_,function(x) x[1,]) %>% dplyr::bind_rows(.)
    df_$missing_pct<-as.numeric(df_$missing_pct)
    df_$missing_pct<-round(df_$missing_pct,0)
    df_$treatment<-as.factor(df_$treatment)
    
    d2<-dplyr::bind_rows(df_)%>% dplyr::filter(treatment=="vehicle") 
    d2$uniqueID<-as.factor(d2$uniqueID)
    d2$sample_name<-as.factor(d2$sample_name)
    d2$C<-as.factor(d2$C)
    d2$rank<-as.factor(d2$rank)
    
    
    d1<-dplyr::bind_rows(df_)%>% dplyr::filter(treatment=="treated") 
    d1$uniqueID<-as.factor(d1$uniqueID)
    d1$sample_name<-as.factor(d1$sample_name)
    d1$C<-as.factor(d1$C)
    d1$rank<-as.factor(d1$rank)
    
    
    d3<-rbind(d1,d2)
    
    d3<-tidyr::pivot_wider(d3,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d1<-pivot_wider(d1,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d2<-pivot_wider(d2,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    
    d1<-d1%>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    d2<-d2 %>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    d3<-d3 %>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    
    rd1<-d1$rank
    rd2<-d2$rank
    rd3<-d3$rank
    
    d1<-1+d1[,lapply(d1,class)=="numeric"]
    d2<-1+d2[,lapply(d2,class)=="numeric"]
    d3<-1+d3[,lapply(d3,class)=="numeric"]
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    colnames(d1) <- paste("missing ", colnames(d1),"%", sep = " ")
    colnames(d2) <- paste("missing ", colnames(d2),"%", sep = " ")
    colnames(d3) <- paste("missing ", colnames(d3),"%", sep = " ")
    
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    p<-upset(d2,names(d2),
             min_degree=1,
             set_sizes=FALSE,
             n_intersections=N,
             stripes='white',
             base_annotations=list(
               '# of Replicates'=intersection_size(
                 counts=TRUE,
               )
             ),
             annotations =list(
               "Ranked Intensity %"=list(
                 aes=aes(x=intersection, fill=d2$rank),
                 geom=list(
                   geom_bar(stat='count', position='fill', na.rm=TRUE),
                   geom_text(
                     aes(
                       label=!!aes_percentage(relative_to='intersection'),
                       color=ifelse(!is.na(rank), 'show', 'hide')
                     ),
                     stat='count',
                     position=position_fill(vjust = .5)
                   ),
                   scale_y_continuous(labels=scales::percent_format()),
                   
                   show_hide_scale,
                   rating_scale
                 )
               )
             )
             
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))
    
    p2<-upset(d1,names(d1),name="Sample bins",
              min_degree=1,
              n_intersections=N,
              set_sizes=FALSE,
              stripes='white',
              base_annotations=list(
                '# of Replicates'=intersection_size(
                  counts=TRUE,
                )
              ),
              annotations =list(
                "Ranked Intensity %"=list(
                  aes=aes(x=intersection, fill=d1$rank),
                  geom=list(
                    geom_bar(stat='count', position='fill', na.rm=TRUE),
                    
                    geom_text(
                      aes(
                        label=!!aes_percentage(relative_to='intersection'),
                        color=ifelse(!is.na(rank), 'show', 'hide')
                      ),
                      stat='count',
                      position=position_fill(vjust = .5)
                    ),
                    scale_y_continuous(labels=scales::percent_format()),
                    
                    show_hide_scale,
                    rating_scale
                  )
                )
              )
              
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    
    d<-names(d3)[which(!names(d3)=="rank")]
    dr<-which(names(d3)=="rank")
    d3<-d3 %>% dplyr::select(-rank)
    d<-as.numeric(unlist(str_extract_all(d,"[[:digit:]]+")))
    ints<-names(d3)[c(order(d,decreasing=FALSE))] %>% as.list(.)
    d3$rank<-rd3
    
    p3<-ComplexUpset::upset(d3,names(d3), 
                            n_intersections=N,
                            set_sizes=FALSE,
                            stripes='white',
                            base_annotations=list(
                              '# of Replicates'=intersection_size(
                                counts=TRUE,
                              )
                            ),
                            annotations =list(
                              "Ranked Intensity %"=list(
                                aes=aes(x=intersection, fill=d3$rank),
                                geom=list(
                                  geom_bar(stat='count', position='fill', na.rm=TRUE),
                                  
                                  geom_text(
                                    aes(
                                      label=!!aes_percentage(relative_to='intersection'),
                                      color=ifelse(!is.na(rank), 'show', 'hide')
                                    ),
                                    stat='count',
                                    position=position_fill(vjust = .5)
                                  ),
                                  scale_y_continuous(labels=scales::percent_format()),
                                  
                                  show_hide_scale,
                                  rating_scale
                                )
                              )
                            ) 
                            
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    return(list(plot(p),plot(p2),plot(p3)))
  }
}

