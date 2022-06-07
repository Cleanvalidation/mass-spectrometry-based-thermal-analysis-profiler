volcano_data<-function(f,Trilinear=FALSE,Splines=FALSE,Sigmoidal=TRUE,Peptide=FALSE,benchmark=TRUE,Fhist=TRUE,labels=FALSE,type=NA){
  f<-f %>% purrr::keep(function(x) !class(x)=='try-error')
  f<-data.frame(dplyr::bind_rows(f))
  if(any(names(f)=="sample_id")){
    if(any(names(f)=="sample_id")){
      f<-f %>% dplyr::select(-sample_id)
    }
    f<-f %>% dplyr::rename("sample_id"="sample_id")
  }
  if(any(names(f)=="I3")){
    f<-f %>% dplyr::rename("I"="I3")
  }
  if(isTRUE(benchmark)){
    if(isTRUE(Peptide)){
      
      f<-f %>% 
        dplyr::select(uniqueID,sample_id,sample_name,dTm,p_dTm,AUC,rsq,fStatistic,pValue,pAdj,treatment,Fvals,p_dTm,dRSS,d1,d2,rss1,RSS)
      f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      f<-dplyr::bind_rows(f)
      f$diffexpressed <- "Not Shifted"
      f$diffexpressed[f$dTm > 0 & f$p_dTm<0.05] <- "Stabilized Tm"
      f$diffexpressed[f$dTm < 0& f$p_dTm<0.05] <- "Destabilized Tm"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm > 0 & f$p_dTm<0.05] <- "Stabilized Shift"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm < 0 & f$p_dTm<0.05] <- "Destabilized Shift"
      f$diffexpressed<-as.factor(f$diffexpressed)
      f$targets<-NA
      f$delabel <- NA
      others<-as.character(f$uniqueID[which(f$pAdj < 0.05)])
      f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
      f$delabel[f$uniqueID %in%c("P36507","Q02750",others)] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750",others)])
      flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                          labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
      flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
      
      
      f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f)
      if(isTRUE(benchmark)){
        if(!isTRUE(Fhist)){
          check<-ggplot(data=f,mapping=aes(x=dRSS,y=-log10(pAdj),color=diffexpressed))+geom_point()+
            geom_hline(yintercept=-log10(0.05), col="red")+
            scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
            labs(x=expression(RSS["0"]-RSS["1"]),y=expression(-log["10"]*p["adj"]-value))+
            #geom_label(aes(dTm, pAdj, label = delabel), data = f,color="black")+
            # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
            #                 nudge_x = 1,
            #                 force = 2,
            #                 box.padding = 2,
            #                 segment.alpha = .5)+
            ggtitle(f$sample_name[1])+
            theme(legend.position="bottom", legend.box = "horizontal")+
            xlim(0,max(f$Fvals,na.rm=TRUE))#+
          #ylim(0,4) 
          if(isTRUE(labels)){
            if(type=="targets" | is.na(type)){
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)+
                xlim(-max(sqrt(f$dRSS),na.rm=TRUE),max(sqrt(f$dRSS),na.rm=TRUE))
            }else{
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }
          }
          
        }else{
          
          check<-ggplot(f)+
            geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) +
            theme_bw() +
            coord_cartesian(xlim=c(0,10))+
            ggplot2::xlab("F-values")+
            xlim(0,max(f$Fvals))
        }
        
      }
      
    }else if(!isTRUE(Peptide)){
      
      f<-dplyr::bind_rows(f) %>% 
        dplyr::select(uniqueID,sample_id,sample_name,dTm,AUC,rsq,fStatistic,pValue,pAdj,treatment,Fvals,p_dTm,dRSS,d1,d2,rss1,RSS)
      f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      f<-dplyr::bind_rows(f)
      f$diffexpressed <- "Not Shifted"
      f$diffexpressed[f$dTm > 0] <- "Stabilized Tm"
      f$diffexpressed[f$dTm < 0] <- "Destabilized Tm"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm > 0] <- "Stabilized Shift"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm < 0] <- "Destabilized Shift"
      f$diffexpressed<-as.factor(f$diffexpressed)
      f$targets<-NA
      f$delabel <- NA
      others<-as.character(f$uniqueID[which(f$pAdj < 0.05)])
      f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
      f$delabel[f$uniqueID %in%c("P36507","Q02750",others)] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750",others)])
      flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                          labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
      flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
      
      f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f)
      if(isTRUE(benchmark)){
        if(!isTRUE(Fhist)){
          check<-ggplot(data=f,mapping=aes(x=dRSS,y=-log10(pAdj),color=diffexpressed))+geom_point()+
            geom_hline(yintercept=-log10(0.05), col="red")+
            scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
            labs(x=expression(RSS["0"]-RSS["1"]),y=expression(-log["10"]*p["adj"]-value))+
            #geom_label(aes(dTm, pAdj, label = delabel), data = f,color="black")+
            # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
            #                 nudge_x = 1,
            #                 force = 2,
            #                 box.padding = 2,
            #                 segment.alpha = .5)+
            ggtitle(f$sample_name[1])+
            theme(legend.position="bottom", legend.box = "horizontal")+
            coord_cartesian(xlim=c(0,max(f$dRSS,na.rm=TRUE)),ylim=c(0,max(-log10(f$pAdj),na.rm=TRUE)))
          if(isTRUE(labels)){
            if(type=="targets" | is.na(type)){
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }else{
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }
          }
          
        }else{
          
          check<-ggplot(f)+
            geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) +
            theme_bw() +
            coord_cartesian(xlim=c(0,10))+
            ggplot2::xlab("F-values")+
            xlim(0,max(f$Fvals))
          return(check)
        }
        
      }
      return(check)
      
    }
  }else if(isFALSE(benchmark)){
    if(isTRUE(Peptide)){
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-f %>% 
        dplyr::group_split(uniqueID,treatment,sample_id)
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(Tm=try(stats::approx(.$I,.$C, xout=(min(.$I,na.rm=TRUE)+(0.5*(max(.$I, na.rm=TRUE)-min(.$I, na.rm=TRUE)))),n=100)$y)))
      f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID)
      
      
    }else if(!isTRUE(Peptide)){
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-data.frame(f)
      f$sample_name<-as.factor(f$sample_name)
      f$sample_id<-as.factor(f$sample_id)
      f<-f %>%
        dplyr::group_split(uniqueID,treatment,sample_id)
      f<-f %>% purrr::keep(function(x) nrow(x)>0)
      
      f<-dplyr::bind_rows(f)
      
      
    }else if (isTRUE(Trilinear)){#if this is a trilinear result
      #f<-f %>% dplyr::group_split(uniqueID,treatment)
      # f<-lapply(f,function(x) dplyr::bind_rows(x))
      # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-dplyr::bind_rows(f) %>%
        dplyr::group_by(uniqueID,treatment,sample_id) %>% # group individual curve data to calculate individual Tm
        dplyr::mutate(sample_name=as.factor(sample_name),
                      treatment=as.factor(treatment),
                      Tm=with(.,stats::approx(I,C, xout=min(I,na.rm=TRUE)+(0.5*(max(I, na.rm=TRUE)-min(I, na.rm=TRUE))))$y)) %>% 
        dplyr::select(-rsq,-CI,-data) %>% dplyr::ungroup(.)
      f<-f %>% 
        dplyr::group_split(uniqueID)
      
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(
        # v_Tm=mean(x$Tm[x$treatment == "vehicle"],na.rm=TRUE),
        #             t_Tm=mean(x$Tm[x$treatment == "treated"],na.rm=TRUE),
        dTm=mean(x$Tm[x$treatment == "treated"],na.rm=TRUE)-mean(x$Tm[x$treatment == "vehicle"],na.rm=TRUE),
        #FC=log2(v_Tm/t_Tm),
        variance_equal_vt = var.test(x$Tm ~ x$treatment)$p.value,
        p_dTm= try(ifelse(all(variance_equal_vt < 0.05),
                          t.test(Tm ~ treatment, data = x,
                                 var.equal = ifelse(all(variance_equal_vt<0.05),FALSE,TRUE))$p.value[1],NA)),
        p_dTm = p.adjust(p_dTm,"BH")))
      
      f<-f %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
      f<-dplyr::bind_rows(f) 
      # f<-purrr::map(f,function(x) x %>% dplyr::mutate(p_dTm=calcP(uniqueID,Tm,dTm,300)))
    }
  }else{
    f$sig<-sign(f$dTm)
    fpos<-f %>% dplyr::filter(sig>0)
    fneg<-f %>% dplyr::filter(sig<0)
    f1<-fpos %>%
      dplyr::mutate(pV = (1-stats::pf(sqrt(dRSS),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f1<-f1 %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    
    f2<-fneg %>%
      dplyr::mutate(pV = (1-stats::pf(sqrt(dRSS),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f2<-f2 %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    f<-dplyr::bind_rows(f1,f2)
    f<-f %>%
      dplyr::mutate(pV = (1-stats::pf(log2(Fvals+1),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f<-f %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    
    # ggplot(f)+
    #   geom_density(aes(x=log2(Fvals+1),fill = "steelblue",alpha = 0.5))+ 
    #   geom_line(aes(x=log2(Fvals+1),y= df(log2(Fvals+1),df1=4,df2=8)),color="darkred",size = 1.5)
    
    # 
    # M<-median(f$dRSS,na.rm=TRUE)
    # V<-mad(f$dRSS,na.rm=TRUE)
    # #alternative scaling factor sig0-sq
    # altScale<-0.5*V/M
    # #filter out negative delta rss
    # f<-f %>% dplyr::filter(dRSS>0)
    # #effective degrees of freedom
    # ed1<-MASS::fitdistr(x=f$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # ed2<-MASS::fitdistr(x=f$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # #scale data
    # f <-f %>% 
    #   dplyr::mutate(rssDiff = .$dRSS/altScale,
    #                 rss1 =.$rss1/altScale,
    #                 d1=ed1,
    #                 d2=ed2)
    # #
    # #new F-test
    # f<-f %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
    # 
    # d1<-f$d1
    # d2<-f$d2
    # 
    # #RSS1 numerator
    # f$fNum<-f$dRSS/d1
    # #Rss denominator
    # f$fDen<-f$RSS/d2
    # #     
    # f$fStatistic=f$fNum/f$fDen
    
    f$diffexpressed <- "Not Shifted"
    f$diffexpressed[f$dTm > 2 & f$p_dTm<0.01] <- "Stabilized Tm"
    f$diffexpressed[f$dTm < -2 & f$p_dTm<0.01] <- "Destabilized Tm"
    f$diffexpressed[f$Fvals > f$fStatistic & f$pV<0.05 & f$dTm > 2 & f$p_dTm<0.01] <- "Stabilized Shift"
    f$diffexpressed[f$Fvals > f$fStatistic &f$pV<0.05 & f$dTm < -2 & f$p_dTm<0.01] <- "Destabilized Shift"
    f$diffexpressed<-as.factor(f$diffexpressed)
    f$targets<-NA
    f$delabel <- NA
    others<-as.character(f$uniqueID[which(f$pAdj < 0.01)])
    f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
    f$delabel[f$uniqueID %in%others] <- as.character(f$uniqueID[f$uniqueID %in%others])
    flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                        labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
    flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
    f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
    f<-purrr::map(f,function(x) x[1,])
    f<-dplyr::bind_rows(f)
    if(!isTRUE(Fhist)){
      check<-ggplot(data=f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),color=diffexpressed))+geom_point()+
        #geom_hline(yintercept=-log10(0.01), col="red")+ 
        labs(y=expression(log["2"]*F["vals"]+1),x=expression(sign("k")*sqrt(RSS["0"]-RSS["1"])))+
        scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
        ggtitle(f$sample_name[1])+
        theme(legend.position="bottom", legend.box = "horizontal")+
        xlim(-max(sqrt(f$dRSS),na.rm=TRUE),max(sqrt(f$dRSS),na.rm=TRUE))
      if(isTRUE(labels)){
        if(type=="targets" | is.na(type)){
          check<-check+
            geom_label_repel(f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                             nudge_x = 1,
                             force = 3,
                             box.padding = 3,
                             segment.alpha = .5)
        }else{
          check<-check+
            geom_label_repel(f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                             nudge_x = 1,
                             force = 3,
                             box.padding = 3,
                             segment.alpha = .5)
        }
      }
      
      
    }else{
      check<-ggplot(f)+
        geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) + 
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
    }
    return(check)
  }
  
  if(isTRUE(Trilinear)){
    
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="lm"))
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$treatment=="treated"][1]>x$Tm[x$treatment=="vehicle"][1],"Stabilized","Destabilized"))))
    
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    df_<-dplyr::bind_rows(f) %>% 
      select(uniqueID,FC,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,rsq_greater_than_0.8,model_converged,stabilized) %>% 
      distinct(.)
    df_$diffexpressed <- "Not shifted"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 2 & df_$p_dTm < 0.05] <- "Stabilized Shift"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -2 & df_$p_dTm < 0.05] <- "Destabilized Shift"
    df_$delabel <- NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750","P60033")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750","P60033")])
    
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
      labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
      #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_)+
      geom_text_repel(df_,mapping=aes(Fvals, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")
    
  }else if(isTRUE(Splines)){
    
    f<-dplyr::bind_rows(f) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name)
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(dTm>0 & ptest<0.01,"Stabilized Shift",ifelse(dTm<0 & ptest<0.01,"Destabilized Shift","No")))))
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    
    df_<-f %>% 
      dplyr::select(uniqueID,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,stabilized,ptest,pAdj) %>% 
      distinct(.)
    df_$diffexpressed <- "Not Shifted"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 2 & df_$ptest < 0.01] <- "Stabilized Shift"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -2 & df_$ptest < 0.01] <- "Destabilized Shift"
    
    df_$delabel <- NA
    Names<-as.character(df_$uniqueID[which(df_$pAdj<0.05 & !df_$diffexpressed=="Not Shifted")])
    
    df_$targets<-NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750",Names)] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750",Names)])
    df_$targets[df_$uniqueID %in%c("P36507","Q02750")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750")])
    df_$diffexpressed<-as.factor(df_$diffexpressed)
    flevels<-data.frame(colors=c("blue","black","green"))
    df_<-df_ %>% dplyr::select(-Tm,-rss,-AUC,-rsq,-stabilized)
    #df_<-df_[!duplicated(df_),]
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(ptest),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
      geom_hline(yintercept=-log10(0.01), col="red")+
      scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
      labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
      
      # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
      #                 nudge_x = 1,
      #                 force = 2,
      #                 box.padding = 2,
      #                 segment.alpha = .5)+
      ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")+
      geom_label(aes(x=10,y=max(-log10(ptest),na.rm=TRUE)),label=paste0("S = ",nrow(df_[df_$diffexpressed=="Stabilized Shift",])),show.legend=FALSE)+
      geom_label(aes(x=-10,y=max(-log10(ptest),na.rm=TRUE)),label=paste0("DS = ",nrow(df_[df_$diffexpressed=="Destabilized Shift",])),show.legend=FALSE)+
      xlim(-15,15)
    
    if(isTRUE(labels)){
      if(type=="targets" | is.na(type)){
        check<-check+
          geom_label_repel(df_,mapping=aes(dTm, -log10(ptest),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           nudge_x = 1,
                           force = 3,
                           box.padding = 3,
                           segment.alpha = .5)
      }else{
        check<-check+
          geom_label_repel(df_,mapping=aes(dTm, -log10(ptest),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           nudge_x = 1,
                           force = 3,
                           box.padding = 3,
                           segment.alpha = .5)
      }
    }
    
    
  }
  return(check)
}