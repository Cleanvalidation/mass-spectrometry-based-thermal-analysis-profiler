UpSet_curves<-function(f,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE){
  f<-f %>% purrr::keep(function(x) !class(x)=='try-error')
  if(isTRUE(Peptide)){
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=as.factor(sample_name),treatment=as.factor(treatment)) %>%
      dplyr::select(-C,-I,-temp_ref,-CV_pct,-missing_pct) %>% dplyr::filter(!is.na(sample_name)) %>%
      dplyr::group_split(uniqueID,treatment,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>% group_by(sample_id,treatment) %>% dplyr::summarise(uniqueID=uniqueID,
                                                                                         treatment=treatment,
                                                                                         sample_name=sample_name,
                                                                                         p_dTm=p_dTm,
                                                                                         sample_id=sample_id,
                                                                                         Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                                                                         rss=mean(rss,na.rm=TRUE),
                                                                                         rsq=mean(rsq,na.rm=TRUE),
                                                                                         AUC=mean(AUC,na.rm=TRUE)) %>% 
                    ungroup(.) %>% distinct(.))
  }else if(!isTRUE(Peptide)){
    f<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=sample_name,treatment=as.factor(treatment)) %>% 
      dplyr::select(-rank,-C,-I,-temp_ref,-CV_pct,-missing_pct) %>%
      dplyr::group_split(uniqueID,treatment,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>%
                    group_by(sample_id,treatment) %>%
                    dplyr::summarise(uniqueID=uniqueID,
                                     treatment=treatment,
                                     sample_name=sample_name,
                                     sample_id=sample_id,
                                     p_dTm=p_dTm,
                                     Coverage=Coverage,
                                     MW_kDa=MW_kDa,
                                     Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                     rss=mean(rss,na.rm=TRUE),
                                     rsq=mean(rsq,na.rm=TRUE),
                                     AUC=mean(AUC,na.rm=TRUE)) %>%
                    ungroup(.) %>% distinct(.)) 
  }else if (isTRUE(Trilinear)){#if this is a trilinear result
    #f<-f %>% dplyr::group_split(uniqueID,treatment)
    # f<-lapply(f,function(x) dplyr::bind_rows(x))
    # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=sample_name,treatment=as.factor(treatment)) %>%
      dplyr::select(-rsq,-data) %>%
      dplyr::group_split(uniqueID,treatment,sample_name)
    
    
    f<-purrr::map(f,function(x) x %>% dplyr::summarise(uniqueID=uniqueID,
                                                       treatment=treatment,
                                                       M1=M1,
                                                       Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                                       rss=sum(rss,na.rm=TRUE),
                                                       rsq=mean(rsq,na.rm=TRUE),
                                                       AUC=mean(AUC,na.rm=TRUE),
                                                       rssDiff=rssDiff,
                                                       sample_name=x$sample_name,
                                                       Fvals=Fvals,
                                                       rssDiff=rssDiff,
                                                       pV=pV,
                                                       pAdj=pAdj) %>% 
                    ungroup(.) %>% distinct(.))
  }
  
  
  if(isTRUE(Trilinear)){
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="lm"))
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$treatment=="treated"][1]>x$Tm[x$treatment=="vehicle"][1],"Stabilized","Destabilized")),
                                                    dTm=x$Tm[x$treatment=="treated"][1]-x$Tm[x$treatment=="vehicle"][1]))
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    
    df_1<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=as.character(sample_name)) %>% dplyr::group_split(uniqueID,sample_name)
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,model_converged,stabilized,rsq_greater_than_0.8) %>% 
      pivot_wider(names_from=sample_name,values_from=c(model_converged)) %>% distinct(.)
    
  }else if(isTRUE(Splines) & !isTRUE(Peptide)){
    
    f<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(Coverage))%>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm>0 & p_dTm<0.05,"Stabilized","Destabilized")))
    
    f<-purrr::map(f,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    df_TPP<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(sample_name)) %>% 
      dplyr::select(uniqueID,sample_name,Tm,rss,AUC,p_dTm,stabilized,sample_id) %>% 
      distinct(.)
    
    df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,sample_id,Tm) %>% distinct(.)
    df_1<-dplyr::bind_rows(df_1)
    df_2<-df_1
    colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors1$sample_name<-as.factor(colors1$sample_name)
    colors1$sample_name<-levels(colors1$sample_name)
    colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
    colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
    colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
    df_1<-df_1[!duplicated(df_1),]
    
    if(length(unique(df_1$sample_name))==1){
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = 0) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }else{
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = ifelse(sample_name == "C_F_E", 0.2, 0)) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }
    
    # pdf("Number_of_curves_upset_splines_Protein.pdf",encoding="CP1253.enc",compress=FALSE,width=6.12,height=4.02)
    # check1
    # dev.off()
    
    df_1<-dplyr::bind_rows(df_2)
    
    df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
      pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
    
    
    
    
    #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
    
    IDs<-as.factor(df_1$uniqueID)
    if(!ncol(df_1)==2){
      df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
      df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
    }
    df_1<-cbind(IDs,df_1)
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    df_1<-na.omit(df_1)
    
    #keep the first row out of redundant peptide group data
    
    #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
    # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
    #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
    #                           RSQ=mean(Rsq,na.rm=TRUE)))
    # f<-dplyr::bind_rows(f) 
    # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
    #                        rss_n=ifelse(treatment=="null",RSS,NA))
    # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
    
    #select columns of interest
    
    
    #,rsq_greater_than_0.8,stabilized
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    
    if(isTRUE(filter)){
      check<-list()
      
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=FALSE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),     
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8
                   
      )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(filter=="TRUE","filtered)","unfiltered)")))
      
      level_data=rev(levels(check$data$intersection))
      #colors$sample_name<-as.character(levels(colors$sample_id))
      colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors$sample_name<-as.factor(colors$sample_name)
      colors$sample_name<-levels(colors$sample_name)
      colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')
      colors$sample_name<-levels(colors$sample_name)
      colors<-dplyr::bind_rows(colors) %>% dplyr::filter(sample_name %in% level_data)
      
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      #check upset plots for intersections
      check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
      
      if(str_count(check_intersections$IDs[1],"-")<7){
        queries<-queries
      }else{
        c(queries,list(
          upset_query(
            intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
            color='black',
            fill='black',
            only_components='# of fitted curves'
          )))
      }
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),               
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8,
                   queries<-queries
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines","trilinear")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      
      df_<-df_[apply(df_!=0, 1, all),]
      return(list(check,check1,df_))
    }else{#if the data isnt filtered
      
      f<-dplyr::bind_rows(f) %>% 
        distinct(.) %>% dplyr::group_split(uniqueID,sample_name,sample_id)
      
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm>0 & p_dTm<0.01,"Stabilized","Destabilized")))
      
      f<-purrr::map(f,function(x) x[1,])
      
      f<-dplyr::bind_rows(f) %>% dplyr::select(-sample_id) %>% 
        f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      df_TPP<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(sample_name)) %>% 
        dplyr::select(uniqueID,sample_name,Tm,rss,AUC,p_dTm,stabilized) %>% 
        distinct(.)
      
      df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,Tm) %>% distinct(.)
      df_1<-dplyr::bind_rows(df_1)
      df_2<-df_1
      colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors1$sample_name<-as.factor(colors1$sample_name)
      colors1$sample_name<-levels(colors1$sample_name)
      colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
      colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
      colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      
      df_1<-df_1[!duplicated(df_1),]
      if(length(unique(df_1$sample_name))==1){
        check1<-df_1 %>% count(sample_name) %>%
          mutate(focus = 0) %>%
          ggplot() +
          ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
          ggforce::theme_no_axes()+
          scale_fill_manual(values = colors1$hex,aesthetics="fill")+
          xlim(-1.1,1.45)+
          ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
          ggtitle("Number of fitted curves")+
          theme(legend.position="bottom", legend.box = "horizontal")
        return(check1)
      }else{
        check1<-df_1 %>% count(sample_name) %>%
          mutate(focus = ifelse(sample_name == "C_F_E", 0.2, 0)) %>%
          ggplot() +
          ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
          ggforce::theme_no_axes()+
          scale_fill_manual(values = colors1$hex,aesthetics="fill")+
          xlim(-1.1,1.45)+
          ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
          ggtitle("Number of fitted curves")+
          theme(legend.position="bottom", legend.box = "horizontal")
      }
      
      df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,Tm) %>% 
        pivot_wider(names_from=sample_name,values_from=c(Tm)) %>% distinct(.)
      
      df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
      df_1 <- df_1 %>%
        mutate_if(sapply(df_1, is.factor), as.numeric)
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1<-mutate_all(df_1, ~replace(., !is.na(.) & is.numeric(.), "TRUE"))
      df_1 <- mutate_all(df_1, ~replace(., is.na(.), "FALSE"))
      df_1 <- mutate_all(df_1, ~replace(., is.character(.), as.logical(.)))
      
      #keep the first row out of redundant peptide group data
      
      #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
      # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
      #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
      #                           RSQ=mean(Rsq,na.rm=TRUE)))
      # f<-dplyr::bind_rows(f) 
      # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
      #                        rss_n=ifelse(treatment=="null",RSS,NA))
      # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
      
      #select columns of interest
      
      
      #,rsq_greater_than_0.8,stabilized
      rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                       values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
      
      check<-list()
      
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% "uniqueID"],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=FALSE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),     
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8
                   
      )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(filter=="TRUE","filtered)","unfiltered)")))
      
      level_data=rev(levels(check$data$intersection))
      #colors$sample_name<-as.character(levels(colors$sample_id))
      colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors$sample_name<-as.factor(colors$sample_name)
      colors$sample_name<-levels(colors$sample_name)
      colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
      
      colors<-dplyr::bind_rows(colors) %>% dplyr::filter(sample_name %in% level_data)
      
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      #check upset plots for intersections
      check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
      
      if(str_count(check_intersections$IDs[1],"-")<7){
        queries<-queries
      }else{
        c(queries,list(
          upset_query(
            intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
            color='black',
            fill='black',
            only_components='# of fitted curves'
          )))
      }
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),               
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8,
                   queries<-queries
                   
                   
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines","trilinear")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      
      df_<-df_[apply(df_!=0, 1, all),]
    }
    return(list(check,check1,df_))
  }else{#if this is a peptide spline file
    
    f1<-dplyr::bind_rows(f) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_id,sample_name)
    
    
    f<-purrr::map(f1,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    df_TPP<-dplyr::bind_rows(f) %>%dplyr::filter(!is.na(sample_name)) %>% 
      dplyr::select(uniqueID,sample_name,rss,AUC,p_dTm,Tm) %>% 
      distinct(.)
    
    df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,Tm) %>% distinct(.)
    
    df_1<-dplyr::bind_rows(df_1)
    df_2<-df_1
    colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors1$sample_name<-as.factor(colors1$sample_name)
    colors1$sample_name<-levels(colors1$sample_name)
    colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
    colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
    colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
    
    
    
    df_1<-df_2[!duplicated(df_2),]
    if(length(unique(df_1$sample_name))==1){
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = 0) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
      return(check1)
    }else{
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }
    
    
    # pdf("Number of fitted_curves_Peptide_filt.pdf",encoding="CP1253.enc",compress=TRUE,width=5.31,height=4.02)
    # check1
    # dev.off()
    
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,Tm) %>% 
      pivot_wider(names_from=sample_name,values_from=c(Tm)) %>% distinct(.)
    
    df_<-df_ %>% dplyr::mutate(uniqueID=as.character(uniqueID))
    df_ <- df_ %>%
      mutate_if(sapply(df_, is.factor), as.numeric)
    df_$uniqueID<-as.factor(df_$uniqueID)
    df_<-mutate_all(df_, ~replace(., !is.na(.) & is.numeric(.), "TRUE"))
    df_ <- mutate_all(df_, ~replace(., is.na(.), "FALSE"))
    df_ <- mutate_all(df_, ~replace(., is.character(.), as.logical(.)))
    
    #keep the first row out of redundant peptide group data
    
    #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
    # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
    #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
    #                           RSQ=mean(Rsq,na.rm=TRUE)))
    # f<-dplyr::bind_rows(f) 
    # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
    #                        rss_n=ifelse(treatment=="null",RSS,NA))
    # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
    
    #select columns of interest
    
    
    #,rsq_greater_than_0.8,stabilized
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    
    check<-list()
    
    check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                 #min_degree=6,
                 set_sizes=FALSE,
                 n_intersections=10,
                 min_degree=1,
                 encode_sets=FALSE,
                 stripes='white',
                 # intersections=list(
                 #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                 #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                 #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                 #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                 #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                 #   
                 # ),
                 # mode = 'inclusive_union',
                 #mode=mode,
                 # intersections=list(
                 #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                 #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                 #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                 #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                 #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                 #   
                 # ),
                 #sort_intersections_by="cardinality",
                 base_annotations=list(
                   '# of fitted curves'=
                     intersection_size(
                       # mode='inclusive_intersection',
                       counts=TRUE
                       # mapping=aes(fill='inclusive_intersection')
                     )+
                     ggplot2::ylab("# of fitted curves")
                 ),     
                 # base_annotations=list(
                 #   'Number of protein events'=upset_annotate(
                 #     '..count..',
                 #     list(
                 #       geom_bar(aes(fill=df_$stabilized)),
                 #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                 #     )
                 #   )
                 # )#,
                 
                 # annotations =list(
                 #   'Stabilized Percentage'=list(
                 #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                 #     geom=list(
                 #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                 #       geom_text(
                 #         aes(
                 #           label=!!aes_percentage(relative_to='group'),
                 #           group=df_$stabilized
                 #         ),
                 #         stat='count',
                 #         position=position_fill(vjust = .5)
                 #       ),
                 #       scale_y_continuous(labels=scales::percent_format()),
                 #       rating_scale
                 #     )
                 #   )
                 # ,
                 # 'Tm'=(
                 #   # note that aes(x=intersection) is supplied by default and can be skipped
                 #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                 #     # checkout ggbeeswarm::geom_quasirandom for better results!
                 #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                 #     geom_violin(alpha=0.5, na.rm=TRUE)
                 # )
                 #),
                 width_ratio=0.1,
                 height_ratio=0.8
                 
    )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(isTRUE(filter),"filtered)","unfiltered)")))
    
    
    level_data=rev(levels(check$data$intersection))
    check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
    
    #colors$sample_name<-as.character(levels(colors$sample_id))
    colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors$sample_name<-levels(colors$sample_name)
    colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')
    colors$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors$sample_name))]
    colors$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors$sample_name))]
    
    colors<-colors %>% dplyr::filter(sample_name %in% check_intersections$IDs)
    queries=list(
      upset_query(
        intersect=colors$sample_name[1],
        color=colors$hex[1],
        fill=colors$hex[1],
        only_components='# of fitted curves'
      ),
      upset_query(
        intersect=colors$sample_name[2],
        color=colors$hex[2],
        fill=colors$hex[2],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=as.character(colors$sample_name[3]),
        color=colors$hex[3],
        fill=colors$hex[3],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[4],
        color=colors$hex[4],
        fill=colors$hex[4],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[5],
        color=colors$hex[5],
        fill=colors$hex[5],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[6],
        color=colors$hex[6],
        fill=colors$hex[6],
        only_components='# of fitted curves'
      ),
      upset_query(
        intersect= as.character(colors$sample_name[7]),
        color=as.character(colors$hex[7]),
        fill=as.character(colors$hex[7]),
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=as.character(colors$sample_name[8]),
        color=as.character(colors$hex[8]),
        fill=as.character(colors$hex[8]),
        only_components='# of fitted curves'
      ))[1:length(colors$sample_name)]
    #check upset plots for intersections
    
    if(str_count(check_intersections$IDs[1],"-")<7){
      queries<-queries
    }else{
      queries=c(queries,list(
        upset_query(
          intersect=c('C_F_E','C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
          color='black',
          fill='black',
          only_components='# of fitted curves'
        )))
    }
    if(!isTRUE(filter)){#if the data is not filtered
      check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),               
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8,
                   queries=queries
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines ","trilinear ")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      df_<-df_[apply(df_!=0, 1, all),]
      
      return(list(check,check1,df_))
      
    }else{#if the peptide data is unfiltered
      check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
                   ),               
                   # base_annotations=list(
                   #   'Number of protein events'=upset_annotate(
                   #     '..count..',
                   #     list(
                   #       geom_bar(aes(fill=df_$stabilized)),
                   #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
                   #     )
                   #   )
                   # )#,
                   
                   # annotations =list(
                   #   'Stabilized Percentage'=list(
                   #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
                   #     geom=list(
                   #       geom_bar(stat='count', position='fill', na.rm=TRUE),
                   #       geom_text(
                   #         aes(
                   #           label=!!aes_percentage(relative_to='group'),
                   #           group=df_$stabilized
                   #         ),
                   #         stat='count',
                   #         position=position_fill(vjust = .5)
                   #       ),
                   #       scale_y_continuous(labels=scales::percent_format()),
                   #       rating_scale
                   #     )
                   #   )
                   # ,
                   # 'Tm'=(
                   #   # note that aes(x=intersection) is supplied by default and can be skipped
                   #   ggplot(mapping=aes(y=log10(df_$dTm)))+
                   #     # checkout ggbeeswarm::geom_quasirandom for better results!
                   #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
                   #     geom_violin(alpha=0.5, na.rm=TRUE)
                   # )
                   #),
                   width_ratio=0.1,
                   height_ratio=0.8,
                   queries=queries
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines ","trilinear ")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      df_<-df_[apply(df_!=0, 1, all),]
      
      return(list(check,check1,df_))
    }
  } 
  
}