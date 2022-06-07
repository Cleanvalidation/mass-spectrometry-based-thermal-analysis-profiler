#Convert to MSStatsTMT
MSStats_converter<-function(df_raw,solvent,ref){
  editPSMs2<-data.frame()
  ref<-as.character(ref)
  editPSMs2<-df_raw %>% dplyr::rename("ProteinName"="Accession",
                                      "PeptideSequence"="Annotated_Sequence",
                                      "Run"="Spectrum.File",
                                      #"PSM"="PSMs_Peptide_ID",
                                      "Channel"="temp_ref",
                                      "Intensity"="value")
  #Condition, Bioreplicate and TechRepMixture need to be filled
  
  editPSMs2$Condition<-ifelse(editPSMs2$Channel==ref,"Norm",0)
  if(nchar(editPSMs2$Run[1])>8){
    editPSMs2$BioReplicate<-ifelse(stringr::str_detect(editPSMs2$Run,solvent)=="TRUE","vehicle","treated")
    editPSMs2<-editPSMs2 %>% 
      dplyr::mutate(Mixture=paste0(ifelse(stringr::str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))
  }else{
    editPSMs2$BioReplicate<-ifelse(stringr::str_detect(editPSMs2$Run,"01"),1,2)
    editPSMs2<-editPSMs2 %>% 
      dplyr::mutate(Mixture=paste0(ifelse(stringr::str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))
    
  }
  editPSMs2$TechRepMixture<-1
  editPSMs2<-editPSMs2 %>% dplyr::select(ProteinName,PeptideSequence,Charge,PSM,Mixture,TechRepMixture,Run,Channel,Condition,BioReplicate,Intensity)
  
  Annotation<-editPSMs2 %>% dplyr::select(Run,TechRepMixture,Channel,Condition,Mixture,BioReplicate)
  Annotation$Fraction<-1
  Annotation<-Annotation %>% dplyr::group_split(Run)
}



#TPP TR Reader
TPPbenchmark<-function(f,overlaps=NA,volcano=TRUE,filters="TPP"){
  f<-f
  f<-list.files(f)
  #read data
  df_TPP<-lapply(f,function(x) read_excel(x,.name_repair = "unique"))
  #extract experiment names
  f<-str_extract_all(f,c("C_F_E","C_F_S","C_nF_E","C_nF_S","nC_F_E","nC_F_S","nC_nF_E","nC_nF_S"))
  f<-lapply(f,function(x) data.frame(sample_name=as.factor(x)))
  
  #join data
  df_TPP<-purrr::map2(df_TPP,f,function(x,y)cbind(x,y))
  
  
  #select columns of interest
  if(filters=="HQ"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(!is.na(meltPoint_MEKi_2) & !is.na(meltPoint_MEKi_2) & !is.na(meltPoint_DMSO_2) &!is.na(meltPoint_DMSO_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_MEKi_2 >0.8 & R_sq_MEKi_1 >0.8 & R_sq_DMSO_2>0.8 &R_sq_DMSO_1>0.8 &
                      plateau_MEKi_2<0.3 & plateau_MEKi_1<0.3 & plateau_DMSO_2<0.3 & plateau_DMSO_1<0.3)
  }else if(filters=="Sigmoidal"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    model_converged_MEKi_2,model_converged_MEKi_1,model_converged_DMSO_2,model_converged_DMSO_1,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(model_converged_DMSO_1=="Yes",model_converged_DMSO_2=="Yes",model_converged_MEKi_1=="Yes",model_converged_MEKi_2=="Yes")
  }else{
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(!is.na(meltPoint_MEKi_2) & !is.na(meltPoint_MEKi_2) & !is.na(meltPoint_DMSO_2) &!is.na(meltPoint_DMSO_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_MEKi_2 >0.8 & R_sq_MEKi_1 >0.8 & R_sq_DMSO_2>0.8 &R_sq_DMSO_1>0.8 &
                      plateau_MEKi_2<0.3 & plateau_MEKi_1<0.3 & plateau_DMSO_2<0.3 & plateau_DMSO_1<0.3 & fulfills_all_4_requirements=="Yes")
  }
  df_TPP$Protein_ID<-as.factor(df_TPP$Protein_ID)
  df_TPP$sample_name<-str_replace(df_TPP$sample_name,"S","\u03A6")
  #get names
  
  df_TPP<-df_TPP %>% dplyr::group_split(sample_name) 
  df_TPP1<-purrr::map(df_TPP,function(x)x %>% dplyr::select(Protein_ID,sample_name,pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2) %>% 
                        pivot_longer(c(pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2),
                                     names_to = c("hi","dTm"),
                                     names_pattern = c("(.+)pVal_(.+)"),
                        ) %>% dplyr::rename("p_dTm"="value"))
  df_TPP1<-purrr::map(df_TPP1,function(x)x %>% dplyr::select(Protein_ID,sample_name,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,p_dTm) %>% 
                        pivot_longer(c(diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2),
                                     names_to = c("hi","set"),
                                     names_pattern = "(.+)diff_(.+)"
                        ) %>% dplyr::rename("dTm"="value") %>% dplyr::select(-hi,-set))
  df_TPP1<-purrr::map(df_TPP1,function(x) x[!is.na(x$dTm),])
  if(isTRUE(volcano)){
    
    df_TPP2<-dplyr::bind_rows(df_TPP1) 
    df_TPP2$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_TPP2$diffexpressed[df_TPP2$dTm > 2 & df_TPP2$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_TPP2$diffexpressed[df_TPP2$dTm < -2 & df_TPP2$p_dTm < 0.05] <- "Destabilized"
    df_TPP2$delabel <- NA
    df_TPP2$delabel[df_TPP2$Protein_ID %in%c("P36507","Q02750")] <- as.character(df_TPP2$Protein_ID[df_TPP2$Protein_ID %in%c("P36507","Q02750")])
    #get overlaps in data
    #df_TPP2$Stroke <- NA
    #df_TPP2$Stroke<-ifelse(df_TPP2$diffexpressed=="Stabilized" | df_TPP2$diffexpressed=="Destabilized",ifelse(df_TPP2$Protein_ID %in% overlaps,TRUE,FALSE),FALSE)
    
    
    
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% distinct(.) %>%  dplyr::group_split(sample_name,Protein_ID)
    df_TPP2<-purrr::map(df_TPP2,function(x) x[1,])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name)
    flevels<-data.frame(colors=c("green","black","red"))
    df_TPP3<-purrr::map(df_TPP2,function(x)
      ggplot(data=x,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+ 
        labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
        xlim(-15,15)+
        theme(legend.position="bottom", legend.box = "horizontal")+
        geom_label(aes(16,y=40),label=paste0("S = ",nrow(x[x$diffexpressed=="Stabilized",])),show.legend=FALSE)+
        geom_label(aes(-16,y=40),label=paste0("DS = ",nrow(x[x$diffexpressed=="Destabilized",])),show.legend=FALSE)+
        scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
        #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
        geom_label_repel(aes(dTm, -log10(p_dTm),label=delabel),color="black",  nudge_x = 1,
                         force = 3,
                         box.padding = 3,
                         segment.alpha = .5)+
        ggtitle(x$sample_name[1])+ylim(-0.1,55))#+
    #geom_point(data=df_TPP2[df_TPP2$Stroke==1,],
    #pch=21, fill=NA, size=4, colour="black", stroke=1)
    
    
    
    return(df_TPP3)
  }
  
  #for number of curves
  
  df_1<-dplyr::bind_rows(df_TPP1)
  
  df_1$sample_name<-str_replace(df_1$sample_name,"S","\u03A6")
  df_1<-dplyr::bind_rows(df_1) 
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-df_1 %>% dplyr::select(p_dTm,uniqueID,sample_name,dTm) %>% distinct(.)
  df_1<-df_1[!duplicated(df_1),]
  df_1<-df_1[!is.na(df_1$dTm),]
  df_1$p_dTm<-base::round(df_1$p_dTm,3) 
  df_1<-dplyr::bind_rows(df_1) %>% distinct(.)%>% dplyr::group_split(uniqueID,sample_name)
  df_1<-purrr::map(df_1,function(x) x[1,])
  df_1<-dplyr::bind_rows(df_1)
  
  colors<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_1)$sample_name)))
  colors$sample_name<-as.character(colors$sample_name)
  colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
  
  check1<-df_1 %>% count(sample_name) %>%
    mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
    ggplot() +
    ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
    ggforce::theme_no_axes()+
    scale_fill_manual(values = colors$hex,aesthetics="fill")+
    xlim(-1.1,1.45)+
    ggplot2::geom_label(mapping=aes(x=colors$x,y=colors$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
    ggtitle("Number of fitted curves")
  if(!isTRUE(filters)){
    pdf("Number_of_fitted_curves_TPP_Protein_only_converging.pdf",encoding="CP1253.enc",compress=TRUE,width=6.12,height=4.02)
    check1
    dev.off()
  }else{
    pdf("Number_of_fitted_curves_TPP_iMAATSA_4_req_Filters.pdf",encoding="CP1253.enc",compress=TRUE,width=6.12,height=4.02)
    check1
    dev.off()
  }
  # 
  # 'Size'=(
  #   intersection_size(
  #     mode=mode,
  #     mapping=aes(fill=exclusive_intersection),
  #     size=0,
  #     text=list(check_overlap=TRUE)
  #   ) + scale_fill_venn_mix(
  #     data=abc_data,
  #     guide=FALSE,
  #     colors=c('A'='red', 'B'='blue', 'C'='green3')
  # 
  df_1<-dplyr::bind_rows(df_TPP)
  df_1$uniqueID<-df_1$Protein_ID
  df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
    pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
  
  
  
  
  #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
  # df_ <- df_ %>%
  #   mutate_if(sapply(df_, is.factor), as.numeric)
  IDs<-df_1$uniqueID
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
  df_1<-cbind(IDs,df_1)
  colors<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
  colors$sample_name<-as.character(colors$sample_name)
  colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
  rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                   values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
  df_1<-na.omit(df_1)
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
                     counts=TRUE,
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   #scale_color_brewer(palette="Dark2")+
                   theme(legend.position = 'none',
                         # axis.text=element_text(size=20, face='bold'),
                         # axis.title=element_text(size=20, face='bold')
                   )+
                   ggplot2::ylab("# of fitted curves")
               ),               # base_annotations=list(
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
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) ",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP" | filters=="HQ","filtered","unfiltered")))
  
  #check upset plots for intersections
  check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
  #filter exclusive intersection colors
  colors<-colors %>% dplyr::filter(sample_name %in% check_intersections$IDs)
  
  #add query
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
               height_ratio=0.8,
               queries<-list(
                 upset_query(
                   intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
                   color='black',
                   fill='black',
                   only_components='# of fitted curves'
                 ),
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
                 # ,
                 # upset_query(
                 #   intersect=colors$sample_name[3],
                 #   color=colors$hex[3],
                 #   fill=colors$hex[3],
                 #   only_components='# of fitted curves'
                 # )
                 # ,
                 # upset_query(
                 #   intersect=colors$sample_name[4],
                 #   color=colors$hex[4],
                 #   fill=colors$hex[4],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect=colors$sample_name[5],
                 #   color=colors$hex[5],
                 #   fill=colors$hex[5],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect=colors$sample_name[6],
                 #   color=colors$hex[6],
                 #   fill=colors$hex[6],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect= as.character(colors$sample_name[7]),
                 #   color=as.character(colors$hex[7]),
                 #   fill=as.character(colors$hex[7]),
                 #   only_components='# of fitted curves'
                 # )
                 # ,
                 # upset_query(
                 #   intersect=as.character(colors$sample_name[8]),
                 #   color=as.character(colors$hex[8]),
                 #   fill=as.character(colors$hex[8]),
                 #   only_components='# of fitted curves'
                 # )
               )
               
               
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) (",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP","additional filters)",ifelse(filters=="HQ","high-quality filtered)","unfiltered)"))))
  if(!isTRUE(filters)){
    pdf("Upset_fitted_curves_TPP_Protein_only_converging.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
    check
    dev.off()
  }else{
    pdf("Upset_fitted_curves_TPP_iMAATSA_4_req_additional_Filters.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
    check
    dev.off()
  }
  return(list(check,check1))
  
  
}
TPPbenchmark_generic<-function(f,overlaps=NA,volcano=TRUE,filters="TPP",Peptide=FALSE,filter=TRUE){
  f<-f
  f<-list.files(f)
  #read data
  df_TPP<-lapply(f,function(x) read_excel(x,.name_repair = "unique"))
  #extract experiment names
  f<-str_extract_all(f,c("C_F_E","C_F_S","C_nF_E","C_nF_S","nC_F_E","nC_F_S","nC_nF_E","nC_nF_S"))
  
  f<-lapply(f,function(x) data.frame(sample_name=as.factor(x)))
  f<-f %>% purrr::keep(function(x) nrow(x)>0)
  #join data
  df_TPP<-purrr::map2(df_TPP,f,function(x,y)cbind(x,y))
  
  
  #select columns of interest
  if(filters=="HQ"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Vehicle_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(!is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) &!is.na(meltPoint_Vehicle_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_Treatment_2 >0.8 & R_sq_Treatment_1 >0.8 & R_sq_Treatment_2>0.8 &R_sq_Vehicle_1>0.8 &
                      plateau_Treatment_2<0.3 & plateau_Treatment_1<0.3 & plateau_Treatment_2<0.3 & plateau_Vehicle_1<0.3)
    
  }else if(filters=="Sigmoidal"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    model_converged_Treatment_2,model_converged_Treatment_1,model_converged_Treatment_2,model_converged_Vehicle_1,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Treatment_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(model_converged_Vehicle_1=="Yes",model_converged_Treatment_2=="Yes",model_converged_Treatment_1=="Yes",model_converged_Treatment_2=="Yes")
    
  }else{
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Treatment_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(!is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) &!is.na(meltPoint_Vehicle_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_Treatment_2 >0.8 & R_sq_Treatment_1 >0.8 & R_sq_Treatment_2>0.8 &R_sq_Vehicle_1>0.8 &
                      plateau_Treatment_2<0.3 & plateau_Treatment_1<0.3 & plateau_Treatment_2<0.3 & plateau_Vehicle_1<0.3 & fulfills_all_4_requirements=="Yes")
  }
  df_TPP$Protein_ID<-as.factor(df_TPP$Protein_ID)
  df_TPP$sample_name<-str_replace(df_TPP$sample_name,"S","\u03A6")
  #get names
  
  df_TPP<-df_TPP %>% dplyr::group_split(sample_name) 
  df_TPP1<-purrr::map(df_TPP,function(x)x %>% dplyr::select(Protein_ID,sample_name,pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2) %>% 
                        pivot_longer(c(pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2),
                                     names_to = c("hi","dTm"),
                                     names_pattern = c("(.+)pVal_(.+)"),
                        ) %>% dplyr::rename("p_dTm"="value"))
  df_TPP1<-purrr::map(df_TPP1,function(x)x %>% dplyr::select(Protein_ID,sample_name,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,p_dTm) %>% 
                        pivot_longer(c(diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2),
                                     names_to = c("hi","set"),
                                     names_pattern = "(.+)diff_(.+)"
                        ) %>% dplyr::rename("dTm"="value") %>% dplyr::select(-hi,-set))
  df_TPP1<-purrr::map(df_TPP1,function(x) x[!is.na(x$dTm),])
  if(isTRUE(volcano)){
    
    df_TPP2<-dplyr::bind_rows(df_TPP1) 
    df_TPP2$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_TPP2$diffexpressed[df_TPP2$dTm > 2 & df_TPP2$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_TPP2$diffexpressed[df_TPP2$dTm < -2 & df_TPP2$p_dTm < 0.05] <- "Destabilized"
    df_TPP2$delabel <- NA
    df_TPP2$delabel[df_TPP2$Protein_ID %in%c("P36507","Q02750")] <- as.character(df_TPP2$Protein_ID[df_TPP2$Protein_ID %in%c("P36507","Q02750")])
    #get overlaps in data
    #df_TPP2$Stroke <- NA
    #df_TPP2$Stroke<-ifelse(df_TPP2$diffexpressed=="Stabilized" | df_TPP2$diffexpressed=="Destabilized",ifelse(df_TPP2$Protein_ID %in% overlaps,TRUE,FALSE),FALSE)
    
    
    
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% distinct(.) %>%  dplyr::group_split(sample_name,Protein_ID)
    df_TPP2<-purrr::map(df_TPP2,function(x) x[1,])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name)
    flevels<-data.frame(treatment=c("No","Destabilized","Stabilized"),
                        colors=c("black","green","purple"))
    
    df_TPP3<-purrr::map(df_TPP2,function(x)
      ggplot(data=x,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+ 
        labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
        xlim(-15,15)+
        theme(legend.position="bottom", legend.box = "horizontal")+
        geom_label(aes(16,y=40),label=paste0("S = ",nrow(x[x$diffexpressed=="Stabilized",])),show.legend=FALSE)+
        geom_label(aes(-16,y=40),label=paste0("DS = ",nrow(x[x$diffexpressed=="Destabilized",])),show.legend=FALSE)+
        scale_color_manual("Stabilization",values=flevels$colors[which(flevels$treatment %in% x$diffexpressed)],labels = flevels$treatment[which(flevels$treatment %in% x$diffexpressed)])+
        #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
        geom_label_repel(aes(dTm, -log10(p_dTm),label=delabel),color="black",  nudge_x = 1,
                         force = 3,
                         box.padding = 3,
                         segment.alpha = .5)+
        ggtitle(x$sample_name[1])+ylim(-0.1,55))#+
    #geom_point(data=df_TPP2[df_TPP2$Stroke==1,],
    #pch=21, fill=NA, size=4, colour="black", stroke=1)
    
    
    
    return(df_TPP3)
  }
  
  #for number of curves
  
  df_1<-dplyr::bind_rows(df_TPP1)
  
  
  df_1<-dplyr::bind_rows(df_1) 
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-df_1 %>% dplyr::select(p_dTm,uniqueID,sample_name,dTm) %>% distinct(.)
  df_1<-df_1[!duplicated(df_1),]
  df_1<-df_1[!is.na(df_1$dTm),]
  df_1$p_dTm<-base::round(df_1$p_dTm,3) 
  df_1<-dplyr::bind_rows(df_1) %>% distinct(.)%>% dplyr::group_split(uniqueID,sample_name)
  df_1<-purrr::map(df_1,function(x) x[1,])
  df_1<-dplyr::bind_rows(df_1)
  
  colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_1)$sample_name)))
  colors1$sample_name<-as.character(colors1$sample_name)[1:length(unique(colors1$sample_name))]
  colors1$sample_name<-as.factor(colors1$sample_name)
  colors1$sample_name<-levels(colors1$sample_name)
  colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
  colors1$x<-c(1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
  colors1$y<-c(1,0.42,-0.38,-0.88,-0.88,-0.38,0.42,1)[1:length(unique(colors1$sample_name))]
  df_TPP1<-df_1
  
  if(!nrow(colors1)==8){
    check1<-df_TPP1 %>% count(sample_name) %>%
      #mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
      ggplot() +
      ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name), stat = "pie") +
      ggforce::theme_no_axes()+
      scale_fill_manual(values = colors1$hex,aesthetics="fill")+
      xlim(-1.1,1.45)+
      ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
      ggtitle("Number of fitted curves")
  }else{
    check1<-df_TPP1 %>% count(sample_name) %>%
      mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
      ggplot() +
      ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name,explode=focus), stat = "pie") +
      ggforce::theme_no_axes()+
      scale_fill_manual(values = colors1$hex,aesthetics="fill")+
      xlim(-1.1,1.45)+
      ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
      ggtitle("Number of fitted curves")
    return(check1)
  }
  # 
  # 'Size'=(
  #   intersection_size(
  #     mode=mode,
  #     mapping=aes(fill=exclusive_intersection),
  #     size=0,
  #     text=list(check_overlap=TRUE)
  #   ) + scale_fill_venn_mix(
  #     data=abc_data,
  #     guide=FALSE,
  #     colors=c('A'='red', 'B'='blue', 'C'='green3')
  # 
  df_1<-dplyr::bind_rows(df_TPP)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
    pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
  
  
  
  
  
  #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
  # df_ <- df_ %>%
  #   mutate_if(sapply(df_, is.factor), as.numeric)
  IDs<-df_1$uniqueID
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
  df_1<-cbind(IDs,df_1)
  
  rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                   values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
  df_1<-na.omit(df_1)
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
                     counts=TRUE,
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   #scale_color_brewer(palette="Dark2")+
                   theme(legend.position = 'none',
                         # axis.text=element_text(size=20, face='bold'),
                         # axis.title=element_text(size=20, face='bold')
                   )+
                   ggplot2::ylab("# of fitted curves")
               ),               # base_annotations=list(
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
               
  )+ggtitle(paste0
            ("Top 10 number of fitted curves (TPP) ",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP" | filters=="HQ","filtered","unfiltered")))
  
  #check upset plots for intersections
  check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
  #filter exclusive intersection colors
  colors<-colors1 %>% dplyr::filter(sample_name %in% check_intersections$IDs)
  
  level_data=rev(levels(check$data$intersection))
  #colors$sample_name<-as.character(levels(colors$sample_id))
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
    ))[1:length(unique(colors$sample_name))]
  if(stringr::str_count(check_intersections$IDs[1],"-")<7){
    queries<-queries
  }else{
    queries<-c(queries,list(
      upset_query(
        intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E','C_F_E'),
        color='black',
        fill='black',
        only_components='# of fitted curves'
      )))
  }
  #add query
  check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
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
               height_ratio=0.8,
               queries=queries
               # ,
               # upset_query(
               #   intersect=colors$sample_name[4],
               #   color=colors$hex[4],
               #   fill=colors$hex[4],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect=colors$sample_name[5],
               #   color=colors$hex[5],
               #   fill=colors$hex[5],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect=colors$sample_name[6],
               #   color=colors$hex[6],
               #   fill=colors$hex[6],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect= as.character(colors$sample_name[7]),
               #   color=as.character(colors$hex[7]),
               #   fill=as.character(colors$hex[7]),
               #   only_components='# of fitted curves'
               # )
               # ,
               # upset_query(
               #   intersect=as.character(colors$sample_name[8]),
               #   color=as.character(colors$hex[8]),
               #   fill=as.character(colors$hex[8]),
               #   only_components='# of fitted curves'
               # )
               
               
               
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) (",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP","additional filters)",ifelse(filters=="HQ","high-quality filtered)",ifelse(isFALSE(filter),"unfiltered)","filtered)")))))
  
  return(list(check,check1))
  
  
}