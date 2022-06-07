spf<-function(spresults,DFN,filters = TRUE){
  # if(!any("missing" %in% names(spresults))){
  #   spresults<-spresults %>% dplyr::mutate(missing=NA,missing_pct=NA,rank=NA)
  #   DFN<-DFN%>% dplyr::mutate(missing=NA,missing_pct=NA,rank=NA)
  # }
  
  spresults<-spresults %>% 
    dplyr::mutate(replicate=as.numeric(replicate)) %>%
    dplyr::group_split(uniqueID)
  DFN<-DFN %>% dplyr::mutate(replicate=as.numeric(replicate)) 
  if(!isTRUE(filters)){
    sl<-purrr::map(seq_len(length(spresults)),function(x) as.numeric({paste(x)})) 
    sp<-purrr::map2(spresults,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))  
    sp<-dplyr::bind_rows(sp)  
    df1<-data.frame(uniqueID = unique(sp$uniqueID))  
    df2<-dplyr::bind_rows(DFN)  
    df2$C<-as.numeric(as.vector(df2$C)) 
    df2$I<-as.numeric(as.vector(df2$I))  
    
    if(any(names(df2)=="id")){
      df2<-df2 %>% dplyr::select(-id)
    }
    
    names<-intersect(names(df2),names(sp))
    df2<-sp %>% left_join(df2, by = names) 
    
    
    Df1<-spresults
  }else{
    #Apply filters 
    #keep the positive AUC differences
    
    spresults<-spresults %>% keep(function(x) mean(x$AUC[x$treatment=="treated"],na.rm=TRUE)>mean(x$AUC[!x$treatment=="vehicle"],na.rm=TRUE))
    
    spresults<-spresults %>% keep(function(x) max(x$lambda)<1)
    
    if (is.null(nrow(spresults))){
      return(warning("all proteins filtered out by AUC and lambda value"))
    }
    #get Tm and RSS differences
    sp<-purrr::map(spresults, function(x) x %>% dplyr::mutate(Tmd= x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "treated"),'Tm'][[1]] - x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "vehicle"),'Tm'][[1]],
                                                              RSSd = sum(x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']) - sum(x[!stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']),
                                                              AUCd = x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "treated"),'AUC'])[[1]]- x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "vehicle"),'AUC'][[1]])
    #conserve list indexes
    sl<-purrr::map(seq_along(length(sp)),function(x) as.numeric({paste(x)}))
    
    #insert list index column
    sp<-purrr::map2(sp,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))
    sp<-dplyr::bind_rows(sp)
    
    sp<-dplyr::arrange(sp,dplyr::desc(AUCd),dplyr::desc(RSSd),dplyr::desc(Tmd)) %>% dplyr::select(uniqueID,id) %>% unique(.) 
    #arrange results by decreasing AUCd, RSSd and Tmd and standardize the order in spresults
    #Df1 holds the model results and stats for splines 
    
    df1<-data.frame(uniqueID = unique(sp$uniqueID))
    df2<-dplyr::bind_rows(DFN) 
    df2$C<-as.numeric(as.vector(df2$C))
    df2$I<-as.numeric(as.vector(df2$I))
    if(any(names(sp) %in% c("missing.x","LineRegion","N"))){
      if(any(names(df2)=="id")){
        df2<-df2 %>% dplyr::select(-id)
      }
      names<-intersect(names(df2),names(sp))
      df2<-sp %>% left_join(df2, by = names) 
    }else{
      if(any(names(df2)=="id")){
        df2<-df2 %>% dplyr::select(-id)
      }
      names<-intersect(names(df2),names(sp))
      df2<-sp %>% left_join(df2, by = names) 
    }
    Df1<-spresults[sp$id]
  }
  ret<-list()
  ret[[1]]<-df1
  ret[[2]]<-df2
  ret[[3]]<-Df1
  return(ret)
}

spCI<-function(i,df1,df2,Df1,df.temps,overlay=TRUE,alpha,residuals=FALSE,simulations=FALSE,CI=TRUE,Peptide=FALSE,CARRIER=TRUE,Protein=Protein,raw=FALSE){
  
  null<-data.frame()
  i<-i
  if(isTRUE(Peptide)){
    df2<-dplyr::bind_rows(df2)
    df2<-df2[,!stringr::str_detect(names(df2),"File.ID|Channel|RT|Confidence|Protein|p|Percolator|DeltaM|Tm|rsq|CC|k_|AUC")]
    df2<-df2[,!names(df2)=="rss"]
    df2<-df2%>% distinct(.) %>% dplyr::group_split(uniqueID,Annotated_Sequence,treatment,C)
    if (.Platform$OS.type=="windows"){
      df2<-mclapply(df2,function(x){y<-x %>% dplyr::mutate(replicate=as.factor(row.names(x)))
      return(y)})
    }else{
      df2<-mclapply(df2,function(x){y<-x %>% dplyr::mutate(replicate=as.factor(row.names(x)))
      return(y)},
      mc.cores=availableCores())
    }
    
    Df1<-dplyr::bind_rows(Df1) %>% distinct(.) %>% dplyr::group_split(uniqueID,Annotated_Sequence,treatment,C)
    if (.Platform$OS.type=="windows"){
      Df1<-parallel::mclapply(Df1,function(x){
        y<-x %>% dplyr::mutate(replicate=as.factor(row.names(x)))
        return(y)
      })
    }else{
      Df1<-parallel::mclapply(Df1,function(x){
        y<-x %>% dplyr::mutate(replicate=as.factor(row.names(x)))
        return(y)
      },
      mc.cores=availableCores())
    }
  }
  df2<-dplyr::bind_rows(df2) %>% distinct(.)
  Df1<-dplyr::bind_rows(Df1) %>% dplyr::group_split(uniqueID) 
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-as.character(df2$uniqueID)
  
  i<-which(df1$uniqueID %in% Protein)
  #get original data
  ###########################################
  df1<-df1$uniqueID[i]
  DF1<-df2[which(df2$uniqueID %in% df1),]
  Df1<-data.frame(Df1[[i]])
  null<-Df1[which(Df1$uniqueID %in% df1 & Df1$treatment %in% "null"),]
  if(nrow(null)==0){
    null<-Df1 %>% dplyr::mutate(treatment="null")
  }
  ###########################################
  DF_f<-df2 %>% subset(uniqueID %in% df1 & treatment %in% "vehicle")
  vehicle<-Df1 %>% subset(uniqueID %in% df1 & treatment %in% "vehicle")
  ###########################################
  DF_f1<-df2%>% subset(uniqueID %in% df1 & treatment %in% "treated")
  treated<-Df1 %>% subset(uniqueID == df1 & treatment == "treated")
  
  ###########################################
  #get confidence intervals for all conditions
  ###########################################
  
  #return fit and confidence intervals
  
  BSVarN<-NA
  BSVar<-NA
  BSVar1<-NA
  BSvarN<-NA
  BSvar<-NA
  BSvar1<-NA
  df2<-df2 %>% dplyr::filter(uniqueID==as.character(df1))
  BSvarN<-df2 %>% dplyr::mutate(treatment=="null")
  BSvar1 <-df2 %>% dplyr::filter(treatment=="treated")
  BSvar <-df2 %>% dplyr::filter(treatment=="vehicle")
  if(nrow(BSvar)==0){
    return(warning(paste0("No vehicle data found for ",Protein)))
  }
  if(nrow(BSvar1)==0){
    return(warning(paste0("No treated data found for ",Protein)))
  }
  BSvarN$treatment<-as.factor(BSvarN$treatment)
  BSvar$treatment<-as.factor(BSvar$treatment)
  BSvar1$treatment<-as.factor(BSvar1$treatment)
  #####try GAM
  #fit penalized splines
  m <- mgcv::gam(I ~ s(C,k=5), data = BSvar , method = "ML")
  m1<-  mgcv::gam(I ~ s(C,k=5), data =BSvar1, method = "ML")
  mn<-  mgcv::gam(I ~ s(C,k=5), data = BSvarN, method = "ML")
  
  #####try GAM
  
  #get some parmeters
  Vb <- vcov(m)
  newd <- with(BSvar, data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 30)))%>% as.data.frame(.)
  pred <- predict(m, newd, se.fit = TRUE)%>% as.data.frame(.)#get confidence intervals
  se.fit <- pred$se.fit
  #get some parmeters
  Vb1<- vcov(m1) 
  newd1<- with(BSvar1,data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 30)))%>% as.data.frame(.)
  pred1<- predict(m1,newd1,se.fit = TRUE) %>% as.data.frame(.)
  se.fit1<- pred1$se.fit
  #generate seed for randomization
  set.seed(42)
  N <- 1000
  #sample n from mvn dist: generates random multivariate normal deviates
  BUdiff <- mgcv::rmvn(N, mu = rep(0, nrow(Vb)), Vb )
  #sample n from mvn dist generates random multivariate normal deviates
  BUdiff1<-  mgcv::rmvn(N, mu = rep(0, nrow(Vb1)),Vb1)
  #random sampling######################################
  Cg <- predict(m, newd, type = "lpmatrix")
  fits <- Cg %*% t(BUdiff)
  nrnd <- 30 #30 random samples
  rnd <- sample(N, nrnd)
  stackFits <- stack(as.data.frame(fits[, rnd]))
  stackFits <- transform(stackFits, C = rep(newd$C, length(rnd)))
  #simulations for treated
  Cg1 <- predict(m1, newd1, type = "lpmatrix")
  fits1 <- Cg1 %*% t(BUdiff1)
  nrnd1 <- 30 #30 random samples
  rnd1 <- sample(N, nrnd1)
  stackFits1 <- stack(as.data.frame(fits1[, rnd1]))
  stackFits1 <- transform(stackFits1, C = rep(newd1$C, length(rnd1)))
  #calculate deviation
  Cg <- predict(m, newd, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  #calculate deviation
  Cg1<- predict(m1,newd1,type = "lpmatrix")
  simDev1<- Cg1%*% t(BUdiff1)
  #calculate abs deviation
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  #calculate abs deviation
  absDev1<- abs(sweep(simDev1,1, se.fit, FUN = "/"))
  #max abs dev
  masd <- apply(absDev, 2L, max)
  #max abs dev
  masd1<- apply(absDev1,2L, max)
  #95% crit values
  crit <- quantile(masd, prob = alpha/2)
  #95% crit values
  crit1<- quantile(masd1,prob = alpha/2)
  #plot CI
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit))
  pred$treatment<-"vehicle"
  pred$treatment<-as.factor(pred$treatment)
  pred$CI<-"vehicle"
  pred$CI<-as.factor(pred$CI)
  
  plot<-ggplot(pred,mapping= ggplot2::aes(x = C,y=fit,color=treatment))+
    geom_point(BSvar, mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(replicate)))+
    geom_ribbon(aes(ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle("")+ylim(-0.4,1.6)+xlim(37,68)+theme(legend.position="bottom")
  
  pred1<- transform(cbind(data.frame(pred1),newd1),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit))
  pred1$treatment<-"treated"
  pred1$treatment<-as.factor(pred1$treatment)
  pred1$CI<-"treated"
  pred1$CI<-as.factor(pred1$CI)
  pred1$AUC<-pracma::trapz(pred1$fit-pred$fit)
  pred1$AUC<-abs(round(pred1$AUC[1],3))
  
  pred1$RSS<- deviance(m1)+deviance(m)
  
  pred1$RSS<- round(pred1$RSS,3)
  #Residuals
  
  pred1$Tm<-round(treated$Tm[1]-vehicle$Tm[1],1)
  #missing values
  miss_v<-data.frame(NA)
  miss_t<-data.frame(NA)
  #max replicates
  
  miss_v<-DF1%>% dplyr::filter(treatment=="vehicle") %>% unique(.)
  miss_t<-DF1%>% dplyr::filter(treatment=="treated") %>% unique(.)
  
  Pred<-data.frame(m$fitted.values,m$residuals)
  names(Pred)<- c("fit","rn")
  Pred$treatment<-as.factor("vehicle")
  BSvar$treatment<-as.factor("vehicle")
  Pred1<-data.frame(m1$fitted.values,m1$residuals)
  names(Pred1)<-c("fit","rn")
  Pred1$treatment<-as.factor("treated")
  Preds<-rbind(Pred,Pred1)
  BSvar1$treatment<-as.factor("treated")
  #get fitted value data
  fitted.values<-data.frame(C=BSvar$M1[[1]]$model$`x$C`,fit=predict(BSvar$M1[[1]],se.fit=TRUE))
  names(fitted.values)<-c("C","fit","se.fit")
  fitted.values1<-data.frame(C=BSvar1$M1[[1]]$model$`x$C`,fit=predict(BSvar1$M1[[1]],se.fit=TRUE))
  names(fitted.values1)<-c("C","fit","se.fit")
  
  #append Tm values on predicted data
  pred1$Tm<-round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)
  if(isTRUE(residuals)){
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=treatment)) +ggplot2::geom_point()+ 
      ggplot2::ggtitle(paste(Df1$uniqueID[1]," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    print(PLrs)
  }
  if(isTRUE(CI)){
    BSVar <-df2 %>% subset(uniqueID == df1 & treatment== "vehicle") # %>% dplyr::mutate(I=mean(I))
    BSVar1 <-df2 %>% subset(uniqueID == df1 & treatment== "treated")# %>% dplyr::mutate(I=mean(I))
    BSVarN<-df2 
    # BSVar<-BSVar[!is.na(BSVar$I),]
    # BSVar1<-BSVar1[!is.na(BSVar1$I),]
    
    BSVar<-BSVar %>% distinct(.)
    BSVar1<-BSVar1 %>% distinct(.)
    
    
    #fit penalized splines
    m <- mgcv::gam(I ~ s(C,k=5), data = BSVar , method = "ML")
    m1<-  mgcv::gam(I ~ s(C,k=5), data =BSVar1, method = "ML")
    mn<-  mgcv::gam(I ~ s(C,k=5), data = BSVarN, method = "ML")
    
    #####try GAM
    #get some parmeters
    Vb <- vcov(m)
    newd <- with(BSVar, data.frame(C = seq(from=min(BSVar$C,na.rm=TRUE), to=max(BSVar$C,na.rm=TRUE), length.out = 10)))%>% as.data.frame(.)
    BSVar <- BSVar %>% dplyr::mutate(fit=list(predict(m, newd, se.fit = TRUE)))
    
    #get some parmeters
    Vb1<- vcov(m1) 
    newd1<- with(BSVar1,data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 10)))%>% as.data.frame(.)
    BSVar1 <- BSVar1 %>% dplyr::mutate(fit=list(predict(m1, newd1, se.fit = TRUE)))
    # if (any(names(BSVar)=="sample.x")){
    #   BSVar<-BSVar %>% dplyr::rename("sample"="sample.x")
    #   
    # }
    # if (any(names(BSVar1)=="sample.x")){
    #   BSVar1<-BSVar1 %>% dplyr::rename("sample"="sample.x")
    #   
    # }
    # 
    #append missing value data
    if(!isTRUE(Peptide)){
      BSVar<-BSVar %>% dplyr::mutate(missing_v=round(BSVar$missing_pct[!is.na(BSVar$missing_pct)][1],0))
      BSVar<-BSVar %>% dplyr::mutate(missing_t=round(BSVar1$missing_pct[!is.na(BSVar1$missing_pct)][1],0))
    }
    p<-data.frame(BSVar$fit[[1]])
    p1<-data.frame(BSVar1$fit[[1]])
    
    fit_v<-p %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                               uprP=fit+(1.96*se.fit),
                               uprS = fit + (crit * se.fit),
                               lwrS = fit - (crit * se.fit),
                               C= seq(min(BSVar$C,na.rm=TRUE), max(BSVar$C,na.rm=TRUE), length = 10),
                               treatment=BSVar$treatment[1],
                               CI=BSVar$treatment[1])
    fit_t<-p1 %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                                uprP=fit+(1.96*se.fit),
                                uprS = fit + (crit1 * se.fit),
                                lwrS = fit - (crit1 * se.fit),
                                C=seq(min(BSVar1$C,na.rm=TRUE), max(BSVar1$C,na.rm=TRUE), length = 10),
                                treatment=BSVar1$treatment[1],
                                CI=BSVar1$treatment[1])
    # 
    # id<-BSVar %>%dplyr::ungroup(.) %>%dplyr::select(sample,replicate)%>% distinct(.)
    # id1<-BSVar1 %>%dplyr::ungroup(.) %>%dplyr::select(sample,replicate)%>% distinct(.)
    # 
    # BSVar$replicate<-as.factor(BSVar$replicate)
    # BSVar1$replicate<-as.factor(BSVar1$replicate)
    # id$replicate<-as.factor(id$replicate)
    # id1$replicate<-as.factor(id1$replicate)
    # BSVar<-BSVar %>% dplyr::right_join(id,by=c("sample","replicate"))
    # BSVar1<-BSVar1%>% dplyr::right_join(id1,by=c("sample","replicate"))
    
    BSVar<-dplyr::bind_rows(BSVar) %>% distinct(.)
    BSVar1<-dplyr::bind_rows(BSVar1) %>% distinct(.)
    if(isTRUE(Peptide)){
      if(any(names(BSVar)=="replicate")&isTRUE(Peptide)){
        BSVar$Replicate<-as.factor(BSVar$replicate)
        BSVar1$Replicate<-as.factor(BSVar1$replicate)
      }else if(any(names(BSVar)=="Charge")){
        BSVar$PSMs_Charge<-as.factor(BSVar$Charge)
        BSVar1$PSMs_Charge<-as.factor(BSVar1$Charge)
      }
      if(any(names(BSVar)=="rank_l")){
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
      }else{
        BSVar$Stroke<-0
        BSVar1$Stroke<-0
      }
      if(!isTRUE(raw)&!all(BSVar$Stroke)==0){#if I isnt raw and there's no rank column for the data
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, fill=NA, size=4, colour="black", stroke=1)+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        
        
        plot<-plot1+
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, fill=NA, size=4)+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.60,2)+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }else if(!isTRUE(raw)&BSVar$Stroke==0){#if I is raw and there's no rank column for the data
        BSVar$Num_PSMs<-as.factor(BSVar$Num_PSMs)#the number of PSMs found for this particular annotated_Sequence
        BSVar1$Num_PSMs<-as.factor(BSVar1$Num_PSMs)
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point()+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        
        
        plot<-plot1+
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point()+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.60,2)+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }else{#if there is a rank column present, use the stroke to label low-ranking Peptides
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, size=4, colour="black", stroke=1)+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        
        
        plot<-plot1+
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, size=4, color=treatment)+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(min(c(min(BSVar$I,na.rm=TRUE),min(BSVar1$I,na.rm=TRUE))),0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }
      
    }else{#if this is a protein file
      if(!any(stringr::str_detect(names(BSVar),"replicate"))){
        BSVar$Replicate<-as.factor(BSVar$replicate)
        BSvar1$Replicate<-as.factor(BSVar1$replicate)
      }
      if(!any(stringr::str_detect(names(BSVar),"Replicate"))){#if there is no replicate column
        BSVar<-BSVar %>% dplyr::group_by(C,treatment) %>% distinct(.) %>%dplyr::group_split()
        BSVar<-lapply(BSVar,function(x) x %>% dplyr::mutate(replicate=row.names(x)))
        BSVar<-dplyr::bind_rows(BSVar) %>% dplyr::ungroup()
        
        BSvar1<-BSvar1  %>% dplyr::group_by(C,treatment) %>% distinct(.) %>%dplyr::group_split()
        BSvar1<-lapply(BSvar1,function(x) x %>% dplyr::mutate(replicate=row.names(x)))
        BSvar1<-dplyr::bind_rows(BSvar1) %>% dplyr::ungroup()
        
        BSVar$Replicate<-as.factor(BSVar$replicate)
        BSvar1$Replicate<-as.factor(BSvar1$replicate)
      }
      if(any(names(BSVar)=="rank_l")){#if this protein file has a rank column
        
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
        
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, size=4, color=treatment, stroke=1)+
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",signif(pred1$Tm[1],3),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        
        
        plot<-plot1+
          ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, size=4, color=treatment, stroke=1)+
          ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred1),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = 2+round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5,
                   yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0,
                   yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.75,0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+xlim(37,68)+
          theme(legend.position="bottom")
        return(plot)
      }else{#if there is no rank column for this protein file
        
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x=C,y=I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          ggplot2::geom_point()+
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")
        
        
        plot<-plot1+
          ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=Replicate))+
          ggplot2::geom_point()+
          ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred1),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values1$C), xend = signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = BSVar1$Tm[1],
                   xend = BSVar1$Tm[1], y = 0, yend = 0.5,
                   colour = "red",linetype=2)+
          ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.75,0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+
          xlim(37,68)+
          theme(legend.position="bottom")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",signif(pred1$Tm[1],3),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        return(plot)
      }
    }
  }
}

spSim<-function(df2,species=9606){#df2 would be the data 
  
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-df2$Accession
  df1<-df2
  df1$uniqueID<-df1$Accession
  df2$uniqueID<-as.character(df2$uniqueID)
  
  ###########################################
  df1<-as.character(df1$uniqueID)
  
  ###########################################
  #get confidence intervals for all conditions
  ###########################################
  
  #return fit and confidence intervals
  
  BSVarN<-NA
  BSVar<-NA
  BSVar1<-NA
  BSvarN<-NA
  BSvar<-NA
  BSvar1<-NA
  
  BSVarN<-df2[df2$uniqueID %in% df1,]
  BSVar1 <-df2[df2$uniqueID %in% df1 & df2$treatment== "treated",]
  BSVar <-df2[df2$uniqueID %in% df1 & df2$treatment== "vehicle",]
  
  BSVarN<-BSVarN %>% dplyr::mutate(treatment="null") 
  BSVar<-BSVar %>% dplyr::mutate(treatment="vehicle")
  BSVar1<-BSVar1 %>% dplyr::mutate(treatment="treated")
  BSVarN$treatment<-as.factor(BSVarN$treatment)
  BSVar$treatment<-as.factor(BSVar$treatment)
  BSVar1$treatment<-as.factor(BSVar1$treatment)
  
  BSVar<-BSVar[!is.na(BSVar$I),]
  BSVar1<-BSVar1[!is.na(BSVar1$I),]
  BSVarN<-BSVarN[!is.na(BSVarN$I),]
  
  
  #load stringdb data
  if(species==9606){
    string_db <- STRINGdb$new( version="10", species=9606,score_threshold=200, input_directory="")
    Known_proteins<-string_db$mp(c("MAP2K1","MAP2K2"))
    
  }else if(species==7955){
    string_db <- STRINGdb$new( version="10", species=7955,score_threshold=900, input_directory="")
    
    Known_proteins<-string_db$mp(c("Q7ZUY3",#histone
                                   "Q6NV46",# stat3, 
                                   "Q68SP3",# stat5a, 
                                   "A4QNT9",# stat5b, 
                                   "Q90Y03",# aldh1a2,
                                   "Q0H2G3",# aldh1a3,
                                   "F1QZU7",# aldh1a2.2, 
                                   "A2BGR9",# aldh1a2.1, 
                                   "F1QLV5",#nqo1,
                                   "F1Q7F3"))
  }
  #CD81<-string_db$mp("CD81")
  #string_db$get_neighbors( c(MEK1, MEK2,CD81) ) [1:10]
  proteins<-string_db$get_proteins()
  graph<-string_db$get_graph()
  
  #get interactors and ensembl ids
  Interactors<-string_db$get_interactions( Known_proteins )
  
  if(length(Interactors)>5000){
    TP<-c(Known_proteins,Interactors,string_db$get_neighbors( Known_proteins ))[1:5000]
  }else{
    TP<-c(Known_proteins,Interactors,string_db$get_neighbors( Known_proteins ))
  }
  TP<-unlist(TP)
  #get protein data for interactors
  xx<-proteins[proteins$protein_external_id %in% TP,]
  symbols <- xx$preferred_name
  ext_id<-xx$protein_external_id
  #get uniprotID
  if(species=="9606"){
    ss<-mapIds(org.Hs.eg.db,symbols,"UNIPROT",'SYMBOL')
  }else if (species ==7955){
    ss <- as.character(mapIds(org.Dr.eg.db, symbols, "UNIPROT", 'SYMBOL'))
  }
  
  ss<-ss[!is.na(ss)]#remove missing values
  #plot ppi network
  example1_mapped <- string_db$map( data.frame(gene=ss), "gene", removeUnmappedRows = TRUE )
  #remove unmapped identifiers
  
  #hits
  if(!nrow(example1_mapped)>400){
    hits <- example1_mapped$STRING_id[1:400]
  }else{
    hits <- example1_mapped$STRING_id
  }
  hits<-hits[!is.na(hits)]
  enrich<-STRINGdb::ppi_enrichment_full(hits,graph)
  
  # # see how many proteins do you have    
  # vcount(graph)
  # 
  # # find top 200 proteins with the highest degree
  # top.degree.verticies <- names(tail(sort(degree(graph)), 200))
  # 
  # # extract the relevant subgraph
  # top.subgraph <- induced_subgraph(graph, top.degree.verticies)
  # 
  # # count the number of proteins in it
  # vcount(top.subgraph)
  
  BSVar_<-BSVar %>% dplyr::filter(uniqueID %in% ss) %>%  dplyr::ungroup(.) %>% dplyr::group_split(uniqueID,sample_name)
  BSVar1_<-BSVar1 %>% dplyr::filter(uniqueID %in% ss) %>% dplyr::ungroup(.) %>% dplyr::group_split(uniqueID,sample_name)
  BSVarN_<-BSVarN %>% dplyr::filter(uniqueID %in% ss) %>% dplyr::ungroup(.) %>% dplyr::group_split(uniqueID,sample_name)
  #add some noise to the original data
  set.seed(1)
  y_data <- purrr::map(BSVar_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$C), 0, 0.05)))
  y_data1 <- purrr::map(BSVar1_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$C), 0, 0.05)))
  y_dataN<-purrr::map(BSVarN_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$C), 0, 0.05)))
  #show simulation results with TP
  fT<-TRUE
  show_results<-TRUE
  #show results with some noise
  spresults<-spstat(y_dataN,y_data,y_data1,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  TP<-spresults %>% dplyr::mutate(outcome="TP")
  #set false positives by overlaying vehicle curves
  spresults<-spstat(y_dataN,y_data,y_data,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  FP<-spresults %>% dplyr::mutate(outcome="FP")
  test<-rbind(TP,FP)
  
  if(any(names(test)=="dRSS")){
    test<-test %>% dplyr::rename("rssDiff"='dRSS')
  }
  # 
  # roc4 <- roc(test$outcome,
  #             test$Fvals, percent=TRUE,
  #             # arguments for ci
  #             ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
  #             # arguments for plot
  #             plot=TRUE,
  #             print.auc=TRUE, show.thres=TRUE)+title(test$sample_name[1])
  # 
  rocs <- pROC::roc(outcome ~Fvals + rssDiff + dTm+AUC,data = test)
  rocs$Fvals$auc[1]<-round(rocs$Fvals$auc[1],3)
  rocs$rssDiff$auc[1]<-round(rocs$rssDiff$auc[1],3)
  rocs$dTm$auc[1]<-round(rocs$dTm$auc[1],3)
  rocs$AUC$auc[1]<-round(rocs$AUC$auc[1],3)
  ROC<- ggroc(rocs)+ggtitle(paste0(str_replace(test$sample_name[1],"S",paste0("\u03A6"))," simulated data"))+
    scale_color_colorblind("Parameters",labels=c("F-stat",
                                                 expression(Delta*RSS),
                                                 expression(Delta*Tm),
                                                 expression(Delta*AUC)))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    ggplot2::annotate("text", x=0.25, y=0.25, label= paste("F-stat AUC = ", rocs$Fvals$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.15, label= paste("RSS AUC = ", rocs$rssDiff$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.05, label= paste("Tm AUC = ", rocs$dTm$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=-0.05, label= paste("AUC = ", rocs$AUC$auc[1]),size=3.5)
  
  
  
  
  #show results without
  spresults<-spstat(BSVarN_,BSVar_,BSVar1_,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  TP<-spresults %>% dplyr::mutate(outcome="TP")
  #set false positives by overlaying vehicle curves
  spresults<-spstat(BSVarN_,BSVar_,BSVar_,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  FP<-spresults %>% dplyr::mutate(outcome="FP")
  test<-rbind(TP,FP)
  if(any(names(test)=="dRSS")){
    test<-test %>% dplyr::rename("rssDiff"='dRSS')
  }
  
  rocs1<- roc(outcome ~Fvals + rssDiff + dTm+AUC,data = test)
  
  rocs1$Fvals$auc[1]<-round(rocs1$Fvals$auc[1],3)
  rocs1$rssDiff$auc[1]<-round(rocs1$rssDiff$auc[1],3)
  rocs1$dTm$auc[1]<-round(rocs1$dTm$auc[1],3)
  rocs1$AUC$auc[1]<-round(rocs1$AUC$auc[1],3)
  ROC1<- ggroc(rocs1)+ggtitle(paste0(str_replace(test$sample_name[1],"S",paste0("\u03A6"))," real data"))+
    scale_color_colorblind("Parameters",labels=c("F-stat",
                                                 expression(Delta*RSS),
                                                 expression(Delta*Tm),
                                                 expression(Delta*AUC)))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    ggplot2::annotate("text", x=0.25, y=0.25, label= paste("F-stat AUC = ", rocs1$Fvals$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.15, label= paste("RSS AUC = ", rocs1$rssDiff$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.05, label= paste("Tm AUC = ", rocs1$dTm$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=-0.05, label= paste("AUC = ", rocs1$AUC$auc[1]),size=3.5)
  
  ROC_list<-list(ROC,ROC1)
  
}

