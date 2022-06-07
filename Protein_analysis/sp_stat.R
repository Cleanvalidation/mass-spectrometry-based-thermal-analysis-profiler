spstat<-function(DF,df,df1,Ftest=TRUE,show_results=TRUE,filters=TRUE,scaled_dof=FALSE,Peptide=FALSE){
  if(!any(names(df)=="C")&!any(names(df)=="temperature")){
    DF<-DF %>% dplyr::rename("C"="temperature")
    df<-df %>% dplyr::rename("C"="temperature")
    df1<-df1 %>% dplyr::rename("C"="temperature")
  }
  
  if(any(class(df)=="list")){
    df<-df %>% purrr::keep(function(x) is.data.frame(x))
    df1<-df1 %>% purrr::keep(function(x) is.data.frame(x))
    DF<-DF %>% purrr::keep(function(x) is.data.frame(x))
    
    df<-dplyr::bind_rows(df)
    df1<-dplyr::bind_rows(df1)
    DF<-dplyr::bind_rows(DF)
  }
  #plot spline results
  df1$treatment<-as.factor(df1$treatment)
  df$treatment<-as.factor(df$treatment)
  DF$treatment<-as.factor(DF$treatment)
  
  df1$I<-as.numeric(as.character(df1$I))
  df$I<-as.numeric(as.character(df$I))
  DF$I<-as.numeric(as.character(DF$I))
  
  df1$C<-as.numeric(as.character(df1$C))
  df$C<-as.numeric(as.character(df$C))
  DF$C<-as.numeric(as.character(DF$C))
  #mutate to get CV values
  DF<-DF %>% dplyr::group_split(C,uniqueID) 
  DF<- purrr::map(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df<-df %>% dplyr::group_split(C,uniqueID,treatment) 
  df<- purrr::map(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df1<-df1 %>% dplyr::group_split(C,uniqueID,treatment) 
  df1<- purrr::map(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  
  #convert to data frame
  
  df<-dplyr::bind_rows(df)
  df1<-dplyr::bind_rows(df1)
  DF<-dplyr::bind_rows(DF)
  
  #convert back to list
  if(any(names(DF)=="replicate")&any(names(DF)=="Fraction")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    
  }else if(any(names(DF)=="Fraction")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,Fraction)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,Fraction)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,Fraction)
  }else if (any(names(DF)=="replicate")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,replicate)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,replicate)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,replicate)
    
  }else{
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID)
  }
  
  #alternative spline fit method : Generalized Additive Models
  #fit penalized splines
  # if(any(names(df[[1]])=="time_point")){
  #   df<-purrr::map(df,function(x) dplyr::bind_rows(x) %>% dplyr::group_split(uniqueID,time_point))
  # }

  if (.Platform$OS.type=="windows"){
    m <- parallel::mclapply(df,fit_gam)
  }else{
    m <- parallel::mclapply(df,fit_gam,mc.cores = future::availableCores())
  }
  
  #remove fitted functions that failed
  m<-m %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(m)==0){
    warning("No fits were possible for the vehicle treatment")
  }
  #parallelize summary of results from gam fit 
  if (.Platform$OS.type=="windows"){
    m <- parallel::mclapply(m,populate_fit)
  }else{
    m <- parallel::mclapply(m,populate_fit,mc.cores = future::availableCores())
  }
  m<-dplyr::bind_rows(m) %>% dplyr::group_split(uniqueID,treatment,sample_id)
  
  
  msub<-purrr::map(m,function(x) x %>% dplyr::select(M1) %>% dplyr::first(.))
  msub<-purrr::map(msub,function(x) data.frame(fitted_values=x[[1]]$fitted.values %>% unique(.)))
  m<-purrr::map2(m,msub,function(x,y)x %>%  
                   dplyr::mutate(fitted_values=ifelse(length(x$M1[[1]]$fitted.values)==length(x$C),
                                                      list(x$M1[[1]]$fitted.values),list(y$fitted_values))))
  if (.Platform$OS.type=="windows"){
    m1 <- parallel::mclapply(df1,fit_gam)
  }else{
    m1 <- parallel::mclapply(df1,fit_gam,mc.cores = future::availableCores())
  }
  m1<-m1 %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(m1)==0){
    warning("No fits were possible for the treated treatment")
  }
  
  #parallelize summary of results from gam fit 
  if (.Platform$OS.type=="windows"){
    m1 <-parallel::mclapply(m1,populate_fit)
  }else{
    m1 <-parallel::mclapply(m1,populate_fit,mc.cores = future::availableCores())
  }
  m1<-dplyr::bind_rows(m1) %>% dplyr::group_split(uniqueID,treatment,sample_id)
  
  
  msub1<-purrr::map(m1,function(x) x %>% dplyr::select(M1) %>% dplyr::first(.))
  msub1<-purrr::map(msub1,function(x) data.frame(fitted_values=x[[1]]$fitted.values %>% unique(.)))
  m1<-purrr::map2(m1,msub1,function(x,y)x %>%  
                    dplyr::mutate(fitted_values=ifelse(length(x$M1[[1]]$fitted.values)==length(x$C),
                                                       list(x$M1[[1]]$fitted.values),list(y$fitted_values))))
  
  
  if (.Platform$OS.type=="windows"){
    mn <- parallel::mclapply(DF,fit_gam)
  }else{
    mn <- parallel::mclapply(DF,fit_gam,mc.cores = future::availableCores())
  }
  mn<-mn %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(mn)==0){
    warning("No fits were possible for the null treatment")
  }
  #free up memory
  DFN<-NA
  df<-NA
  df1<-NA
  #parallelize summary of results from gam fit 
  mn <-furrr::future_map(mn,function(x) populate_fit(x))
  
  mn<-dplyr::bind_rows(mn) %>% dplyr::group_split(uniqueID,treatment,sample_id)
  
  mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M2))=="gam"))
  msubn<-purrr::map(mn,function(x) x %>% dplyr::select(M1) %>% dplyr::first(.))
  msubn<-purrr::map(msubn,function(x) data.frame(fitted_values=x[[1]]$fitted.values %>% unique(.)))
  mn<-purrr::map2(mn,msubn,function(x,y)x %>%  
                    dplyr::mutate(fitted_values=ifelse(length(x$M1[[1]]$fitted.values)==length(x$C),
                                                       list(x$M1[[1]]$fitted.values),list(y$fitted_values))))
  
  
  lm<-length(m)
  lm1<-length(m1)
  
  if(lm<lm1){
    CID<-unique(dplyr::bind_rows(m)$uniqueID)
  }else{
    CID<-unique(dplyr::bind_rows(m1)$uniqueID)
  }
  
  m<-dplyr::bind_rows(m)
  m1<-dplyr::bind_rows(m1)
  mn<-dplyr::bind_rows(mn)
  
  #filter
  m <-m %>% dplyr::filter(uniqueID %in% CID)
  m1<-m1%>% dplyr::filter(uniqueID %in% CID)
  mn<-mn%>% dplyr::filter(uniqueID %in% CID)
  
  
  CID<-dplyr::intersect(unique(m$uniqueID),unique(m1$uniqueID))
  m <-m %>% dplyr::filter(uniqueID %in% CID)
  m1<-m1%>% dplyr::filter(uniqueID %in% CID)
  mn<-mn%>% dplyr::filter(uniqueID %in% CID)
  #split 
  m<-m %>% tibble::as_tibble()%>% dplyr::group_split(uniqueID,sample_id)
  m1<-m1%>% tibble::as_tibble()%>% dplyr::group_split(uniqueID,sample_id)
  mn<-mn%>% tibble::as_tibble()%>% dplyr::group_split(uniqueID,sample_id)
  
  if(length(m[[1]])<length(m1[[1]])){
    m1<-purrr::map(m1,function(x) x[1:length(m[[1]])])
  }else if (length(m[[1]])>length(m1[[1]])){
    m<-purrr::map(m,function(x) x[1:length(m1[[1]])])
  }
  
  m<-dplyr::bind_rows(m) %>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  m1<-dplyr::bind_rows(m1)%>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  mn<-dplyr::bind_rows(mn)%>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  #calculate RSS
  m<-purrr::map2(m,m1,function(x,y) x %>% dplyr::mutate(RSSv=as.numeric(x$rss[x$treatment=="vehicle"][1]),
                                                        RSSt=as.numeric(y$rss[y$treatment=="treated"][1]),
                                                        RSS1=as.numeric(x$rss[x$treatment=="vehicle"][1])+as.numeric(y$rss[y$treatment=="treated"][1]),
                                                        treatment=x$treatment[1],
                                                        Tm=x$Tm[1]
  ))
  
  m1<-purrr::map2(m,m1,function(x,y) y %>% dplyr::mutate(RSSv=x$RSSv[1],
                                                         RSSt=x$RSSt[1],
                                                         RSS1=x$RSS1[1],
                                                         Tm=y$Tm[1]
  ))
  mn<-purrr::map(mn,function(x) x %>% dplyr::mutate(RSSv=NA,
                                                    RSSt=NA,
                                                    RSS0=x$rss,
                                                    Tm=x$Tm[1]
  ))
  if(any(names(m)=="Fraction")){
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,Fraction)
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,Fraction)
    mn<-dplyr::bind_rows(mn) %>%
      dplyr::group_split(uniqueID,Fraction)
  }else if(any(names(m)=="replicate")){
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,replicate)
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,replicate)
    mn<-dplyr::bind_rows(mn) %>%
      dplyr::group_split(uniqueID,replicate)
  }
  m<-purrr::map2(m,mn,function(x,y) x %>%
                   dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  m1<-purrr::map2(m1,mn,function(x,y)x %>%
                    dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  mn<-purrr::map2(m1,mn,function(x,y)y %>%
                    dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  #4031 proteins have rssDiff> 0 
  
  #convert to df and split by uniqueID 
  mean1<-dplyr::bind_rows(m)
  mean1_1<-dplyr::bind_rows(m1)
  mean3<-dplyr::bind_rows(mn)
  #equal column names
  mean1<-mean1 %>% dplyr::select(-RSSv,-RSSt,-RSS1)%>% distinct(.)
  if(any(names(mean1_1)=="RSS1")){
    mean1_1<-mean1_1  %>% dplyr::select(-RSSv,-RSSt,-RSS1)%>% distinct(.)
  }
  mean3<-mean3  %>% dplyr::select(-RSS0,-RSSv,-RSSt) %>% distinct(.)
  
  #filter out proteins with neg RSSdiff
  if(isTRUE(filters)){
    check<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
    
    mean1<-mean1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    mean3<-mean3 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    
    
    mean1<-mean1 %>% dplyr::filter(uniqueID %in% check)
    mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check)
    mean3<-mean3 %>% dplyr::filter(uniqueID %in% check)
    
  }
  
  #split into lists by uniqueID
  mean1<-mean1%>% dplyr::group_split(uniqueID)
  mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
  mean3<-mean3 %>% dplyr::group_split(uniqueID)
  
  
  #Cliff
  results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID,replicate)
  results<-purrr::map(results,function(x) x %>% dplyr::mutate(dTm=x[which(x$treatment=="treated"),]$Tm[1]-x[which(x$treatment=="vehicle"),]$Tm[1]))
  results<-dplyr::bind_rows(results)
  if(!isTRUE(Peptide)){
    results_<-results %>% dplyr::select(uniqueID,sample_name,Tm,sample_id,treatment,M1,M2,missing_pct,C,I) %>%
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name,treatment)
    results_<-purrr::map(results_,function(x)x %>% dplyr::mutate(replicate=row.names(.)))
    results_<-dplyr::bind_rows(results_) %>% dplyr::group_split(uniqueID,sample_name,replicate)
    results_<-results_ %>% purrr::keep(function(x) nrow(x)>1)
    results_<-purrr::map(results_,function(x) x %>%dplyr::group_by(replicate) %>% dplyr::mutate(dTm=(.$Tm[which(.$treatment=="treated")]-.$Tm[which(.$treatment=="vehicle")])))
  }else{#if this is a peptide file
    results_<-dplyr::bind_rows(results) %>%
      dplyr::select(uniqueID,Annotated_Sequence,sample_name,-Tm,sample_id,treatment,C,I,replicate) %>%
      distinct(.) %>% dplyr::group_split(uniqueID,Annotated_Sequence,sample_name,sample_id)
    data<-purrr::map(results_,function(x) x[1,] %>% dplyr::select(-C,-I))
    
    data<-dplyr::bind_rows(data) %>%  dplyr::group_split(uniqueID)
    results_<-dplyr::bind_rows(results_) %>%  dplyr::group_split(uniqueID)
    name<-dplyr::intersect(names(data[[1]]),names(results_[[1]]))
    #join
    results_<-purrr::map2(results_,data,function(x,y)x %>% dplyr::right_join(y,by=name))
    results_1<-dplyr::bind_rows(results) %>% 
      dplyr::mutate(replicate=as.character(replicate)) %>% 
      dplyr::group_split(uniqueID)
    name<-dplyr::intersect(names(results_1[[1]]),names(results_[[1]]))
    results_<-purrr::map2(results_,results_1,function(x,y)x %>% dplyr::right_join(y,by=name))
    # if(isTRUE(Frac)){
    #   results_<-purrr::map(results_,function(x) x %>%dplyr::group_by(replicate,Fraction) %>% dplyr::mutate(dTm=(.$Tm[.$treatment=="treated"]-.$Tm[.$treatment=="vehicle"])))
    # }else{
    results_<-dplyr::bind_rows(results_) %>%dplyr::group_split(uniqueID,replicate)
    results1<-purrr::map(results_,function(x) x %>% dplyr::mutate(dTm=(.$Tm[.$treatment=="treated"][1]-.$Tm[.$treatment=="vehicle"][1])))
    # }
  }
  #only keep data with vehicle and treated values
  
  results_<-dplyr::bind_rows(results_) %>%
    dplyr::group_by(uniqueID) %>%
    dplyr::group_split(.)
  #make sure we keep proteins with both treated and vehicle data
  results_1<-results_ %>%
    purrr::keep(function(x) length(unique(x$treatment))>1)
  results_2<-dplyr::bind_rows(results_1) %>%
    dplyr::group_split(uniqueID,treatment)
  #nest the columns that dont involve Tm t-test calculation 
  nesting<-names(results_2[[1]])[!names(results_2[[1]]) %in% c("uniqueID","sample_id","treatment","Tm","replicate","sample_name","dTm","missing_pct")]
  
  results_2<-purrr::map(results_2,function(x) x %>% tidyr::nest(data=nesting))
  #only keep replicates with Tm values
  results_2<-dplyr::bind_rows(results_2) %>%
    dplyr::group_split(uniqueID)
  results_2<-purrr::map(results_2,function(x) x %>%
                          dplyr::filter(replicate %in% unique(x$replicate[duplicated(x$replicate)])))
  
  results_3<-dplyr::bind_rows(results_2) %>%  dplyr::group_split(uniqueID,replicate)
  
  results_<-results_3 %>% purrr::keep(function(x) any(!is.na(x$dTm)))
  if(length(results_)==0){
    warning(paste0(results_2[[1]]$sample_name[1], " has at least one Tm value missing"))
  }
  results_3<-purrr::map(results_3,function(x) x %>%
                          distinct(.) %>%
                          dplyr::mutate(
                            dTm=x$Tm[x$treatment %in% "treated"][1]-x$Tm[x$treatment %in% "vehicle"][1])%>% dplyr::ungroup(.)) 
  results_<-results_3 %>% purrr::keep(function(x) !is.na(x$dTm[1]))
  results_<-dplyr::bind_rows(results_) %>% dplyr::group_split(sample_name)
  results_<-results_ %>% purrr::keep(.,function(x) length(unique(x$treatment))>1)
  
  results_<-purrr::map(results_,function(x) 
    y<-x %>% 
      dplyr::mutate(hypothesis=ifelse(treatment=="null",as.factor("null"),as.factor("alternative")),
                    p_dTm= try(p.adjust(t.test(Tm ~ treatment, data = .,
                                               var.equal = FALSE,conf.level=0.975)$p.value[1],"BH"))))
  
  
  
  results1<-results_ %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
  
  results<-dplyr::bind_rows(results1) %>% tidyr::unnest(cols=data)
  #results<-results %>%  dplyr::mutate(p_dTm=calcP(uniqueID,Tm,dTm,30000))
  results2<-dplyr::bind_rows(results)$uniqueID
  if(!length(mean3)==length(mean1)){
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    
    CID<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-dplyr::intersect(CID,mean3$uniqueID)
    CID<-dplyr::intersect(CID,results2)
    
    #filter
    mean1 <-mean1 %>% dplyr::filter(uniqueID %in% CID)
    mean1_1<-mean1_1%>% dplyr::filter(uniqueID %in% CID)
    mean3<-mean3%>% dplyr::filter(uniqueID %in% CID)
    
    #split 
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1%>% dplyr::group_split(uniqueID)
    mean3<-mean3%>% dplyr::group_split(uniqueID)
    
    if(!nrow(mean3[[1]])==nrow(mean1[[1]])){
      mean1<-dplyr::bind_rows(mean1) %>% distinct(.)
      mean1_1<-dplyr::bind_rows(mean1_1)%>% distinct(.)
      mean3<-dplyr::bind_rows(mean3)%>% distinct(.)
      
      CID<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-dplyr::intersect(CID,mean3$uniqueID)
      #filter
      mean1 <-mean1 %>% dplyr::filter(uniqueID %in% CID)
      mean1_1<-mean1_1%>% dplyr::filter(uniqueID %in% CID)
      mean3<-mean3%>% dplyr::filter(uniqueID %in% CID)
      
      #split 
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1%>% dplyr::group_split(uniqueID)
      mean3<-mean3%>% dplyr::group_split(uniqueID)
    }
    results<-dplyr::bind_rows(mean1,mean1_1,mean3)
  }
  
  if(isTRUE(Ftest)){
    stats_summary_null<-function(x){
      y =data.frame(uniqueID=x$uniqueID[1],
                    sample_id=x$sample_id[1],
                    sample_name=x$sample_name[1],
                    AUC=x$AUC[1],
                    rsq=x$rsq[1],
                    RSS0=x$rss,
                    se0=summary.gam(x$M1[[1]])$se[1],#standard error
                    edf0=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                    rsq0=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
                    np0=summary.gam(x$M1[[1]])$np[1],#number of parameters
                    rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m=summary.gam(x$M1[[1]])$m[1]) %>% distinct(.)
    }
    #Calculate stats summary
    rss0<-purrr::map(mean3,function(x)stats_summary_null(x))
    #this is the alternative hypothesis (one curve per treatment)
    stats_summary_alt<-function(x,y){
      z =data.frame(uniqueID=x$uniqueID[1],
                    sample_id=x$sample_id[1],
                    sample_name=x$sample_name[1],
                    AUC=x$AUC[1],
                    dTm=y$Tm[1]-x$Tm[1],
                    RSS=x$rss,
                    se=summary.gam(x$M1[[1]])$se[1],#standard error
                    pN1=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                    rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m=summary.gam(x$M1[[1]])$m[1],
                    se1=summary.gam(y$M1[[1]])$se[1],#standard error
                    pN2=summary.gam(y$M1[[1]])$edf[1],#effective degrees of freedom treated
                    rsq=mean(x$M1[[1]]$r.sq[1],y$M1[[1]]$r.sq[1],na.rm=TRUE),#r-squared
                    rdf1=summary.gam(y$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m1=summary.gam(y$M1[[1]])$m[1],
                    npa=sum(summary.gam(x$M1[[1]])$np[1],summary.gam(y$M1[[1]])$np[1],na.rm=TRUE),
                    pNA=sum(summary.gam(x$M1[[1]])$edf[1],summary.gam(y$M1[[1]])$edf[1],na.rm=TRUE)) %>% distinct(.)
      return(z)
    }
    rss1<-purrr::map2(mean1,mean1_1,function(x,y)stats_summary_alt(x,y))
    
    
    
    
    # rss1<-purrr::map2(rss1,mean1,function(x,y) x %>% dplyr::mutate(fit_v=list(ifelse(class(try(mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",
    #                                                                                                              se.fit=TRUE,newdata.guaranteed = TRUE)))=='try-error',NA,
    #                                                                                  mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",se.fit=TRUE,newdata.guaranteed = TRUE)))))
    int<-dplyr::intersect(names(mean3[[1]]),names(rss0[[1]]))
    int1<-dplyr::intersect(names(mean1[[1]]),names(rss1[[1]]))
    
    
    # mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
    # mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
    # 
    #bind predicted values to main data frame
    mean3<-purrr::map2(mean3,rss0,function(x,y)x %>% dplyr::right_join(y))
    mean1<-purrr::map2(mean1,rss1,function(x,y)x %>% dplyr::right_join(y))
    
    
    #params for null and alternative models
    f0<-lapply(mean3,function(x) data.frame(nfv=length(x$M1[[1]]$fitted.values)))#number of measurements
    if(length(mean1[[1]]$M1[[1]]$fitted.values)<length(mean1_1[[1]]$M1[[1]]$fitted.values)){
      f1<-purrr::map2(mean1,mean1_1,function(x,y) sum(length(x$M1[[1]]$fitted.values),
                                                      length(y$M1[[1]]$fitted.values)))#number of measurements alternative
    }else{
      f1<-purrr::map2(mean1,mean1_1,function(x,y) sum(length(x$M1[[1]]$fitted.values),
                                                      length(y$M1[[1]]$fitted.values)))#number of measurements alternative
    }
    f1<-purrr::map(f1,function(x) data.frame(nfvA=sum(x)))
    int<-dplyr::intersect(names(mean3[[1]]),names(f0[[1]]))
    int1<-dplyr::intersect(names(mean1[[1]]),names(f1[[1]]))
    
    # 
    # mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
    # mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
    
    #bind data frames
    mean3<-purrr::map2(mean3,f0,function(x,y) merge(x,y,all=TRUE))
    mean1<-purrr::map2(mean1,f1,function(x,y) merge(x,y,all=TRUE))
    mean1_1<-purrr::map2(mean1_1,f1,function(x,y) merge(x,y,all=TRUE))
    #
    pN<-NA
    pA<-NA
    
    d1<-NA
    d2<-NA
    
    #calculate parameters 
    pN<-purrr::map(mean3,function(x) x[1,])
    pA<-purrr::map(mean1,function(x)x[1,])
    
    
    #degrees of freedom before
    d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pNA-y$np0))
    d2<-purrr::map(pA,function(x)data.frame(d2=x$npa-x$pNA))
    
    #delta RSS
    rssDiff<-purrr::map2(rss0,rss1,function(x,y) data.frame(dRSS=x$RSS0-y$RSS,#rss null minus rss alt
                                                            dTm=y$dTm[1],
                                                            rss1=y$RSS[1]))
    
    d2<-lapply(d2,function(x) x$d2[1])#edf2
    d1<-lapply(d1,function(x) x$d1[1])#edf1
    
    #RSS1 numerator
    Fvals<-purrr::map2(rssDiff,d1,function(x,y) data.frame(fNum=x$dRSS/y[1]))
    #Rss denominator
    Fd<-purrr::map2(rss1,d2,function(x,y) data.frame(fDen=x$RSS/y[1]))
    
    Fvals<-purrr::map2(Fvals,Fd,function(x,y) data.frame(fStatistic=x$fNum/y$fDen))
    Fvals<-purrr::map2(Fvals,d1,function(x,y) x %>% dplyr::mutate(df1=y[1]))
    Fvals<-purrr::map2(Fvals,d2,function(x,y) x %>% dplyr::mutate(df2=y[1]))
    
    Fvals<-purrr::map2(Fvals,rssDiff,function(x,y) cbind(x,y))
    Fvals<-purrr::map(Fvals,function(x) x %>% dplyr::mutate(Fvals=(x$dRSS/x$rss1)*(x$df2/x$df1)))
    #calculate p and p-adj vals
    
    Fvals<-purrr::map(Fvals,tryCatch({function(x) x %>% dplyr::mutate(pValue = 1 - pf(fStatistic, df1 = x$df1, df2 = x$df2),
                                                                      pAdj = p.adjust(pValue,"BH"))
    },error = function(cond){
      message(cond)
    }))
    Fvals<-purrr::map2(Fvals,mean1,function(x,y) x %>% dplyr::mutate(uniqueID=y$uniqueID[1]))
    
    
    mean1<-purrr::map(mean1,function(x) data.frame(x)%>% dplyr::select(-M1) )
    
    mean1<-purrr::map(mean1,function(x) x %>% distinct(.) )
    names1<-dplyr::intersect(names(mean1),names(Fvals))
    #convert results to list
    testResults<-purrr::map2(mean1,Fvals,function(x,y)x%>% dplyr::right_join(y,by=names1))
    testResults<-purrr::map(testResults,function(x) x[1,])
    testResults<-dplyr::bind_rows(testResults)
    mean1<-dplyr::bind_rows(mean1)
    Unscaled<-ggplot(testResults)+
      geom_density(aes(x=fStatistic),fill = "steelblue",alpha = 0.5) +
      geom_line(aes(x=fStatistic,y= df(fStatistic,df1=df1,df2=df2)),color="darkred",size = 1.5) +
      theme_bw() +
      ggplot2::xlab("F-values")+
      ggplot2::xlim(0,0.05)
    # print(Unscaled)
    #scale variables
    M<-median(testResults$dRSS,na.rm=TRUE)
    V<-mad(testResults$dRSS,na.rm=TRUE)
    #alternative scaling factor sig0-sq
    altScale<-0.5*V/M
    #filter out negative delta rss
    testResults<-testResults %>% dplyr::filter(dRSS>0)
    #effective degrees of freedom
    ed1<-tryCatch({MASS::fitdistr(x=testResults$dRSS, densfun = "chi-squared", start = list(df=2))[["estimate"]]},
                  error= function (cond){message(cond)})
    ed2<-tryCatch({MASS::fitdistr(x=testResults$rss, densfun = "chi-squared", start = list(df=2))[["estimate"]]},
                  error = function(cond){message(cond)})
    #scale data
    testScaled <-testResults %>%
      dplyr::mutate(rssDiff = .$dRSS/altScale,
                    rss1 =.$RSS/altScale,
                    d1=ed1,
                    d2=ed2)
    #
    #new F-test
    if(!class(ed1)=="NULL"&!class(ed2)=="NULL"){
      testResults<-testScaled %>% dplyr::mutate(Fvals=(dRSS/rss1)*(d2/d1))
      Fvals<-testResults$Fvals
      d1<-testResults$d1
      d2<-testResults$d2
      
      #scaled values
      TestScaled<-ggplot(testResults)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() + 
        ggplot2::xlab("F-values")+
        ggplot2::xlim(0,0.05)
      # print(TestScaled)
      #Define checked as filtered protein IDs
      check<-testResults$uniqueID
      test<-testResults
      test$d1<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test$d2<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test<-test %>% dplyr::filter(test$pAdj<0.01)
      testS<-ggplot(test)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
      # print(testS)
      
      
      mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
      mean3<-dplyr::bind_rows(mean3)
      mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) 
      if(!isTRUE(Peptide)){
        results1<-dplyr::bind_rows(results) %>%
          dplyr::select(uniqueID,treatment,sample_id,p_dTm)
        nam<-dplyr::intersect(names(results),names(results1))
        results<-results %>% dplyr::right_join(results1,by=nam)
        
        nam<-dplyr::intersect(names(testResults),names(results1))
        testResults<-testResults %>% dplyr::right_join(results1,by=nam)
        
        names<-dplyr::intersect(names(test),names(results1))
        test<-test %>% dplyr::right_join(results1,by=names)
      }else{
        results1<-dplyr::bind_rows(results1) %>% dplyr::select(uniqueID,treatment,sample_id,p_dTm)
        names<-dplyr::intersect(names(results),names(results1))
        results<-results %>% dplyr::right_join(results1,by=names)
        
        names<-dplyr::intersect(names(testResults),names(results1))
        testResults<-testResults %>% dplyr::right_join(results1,by=names)
        
        names<-dplyr::intersect(names(test),names(results1))
        test<-test %>% dplyr::right_join(results1,by=names)
      }
    }
    
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) 
    # results<-results %>% dplyr::group_by(treatment) %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(performance_k5 = list(performance::model_performance(unlist(.$M1))),
    #                 performance_k6 = list(performance::model_performance(unlist(.$M2))))
    
  }
  if(isTRUE(show_results)&exists("testResults")){
    
    if(isTRUE(scaled_dof)){
      return(test)
    }
    
    return(testResults)
  }else{
    
    return(results)
  }
  
}

fit_gam<-function(x){
  y = x %>% dplyr::filter(!is.infinite(I)) %>% dplyr::mutate(M1 = list(tryCatch(mgcv::gam(x$I ~ s(x$C,by = treatment,k=5), data = x , method = "REML"),
                                                                                error = function(cond) {
                                                                                  message("Here's the original error message:")
                                                                                  message(cond)
                                                                                  # Choose a return value in case of error
                                                                                  return(NA)})),
                                                             M2 = list(tryCatch(mgcv::gam(x$I ~ s(x$C,by = treatment,k=6), data = x , method = "REML"),
                                                                                error = function(cond) {
                                                                                  message("Here's the original error message:")
                                                                                  message(cond)
                                                                                  # Choose a return value in case of error
                                                                                  return(NA)}))
  )
  return(y)
}

populate_fit<-function(x) {
  if(any(stringr::str_detect(names(x),"File.ID"))&!any(names(x)=="sample_id")){
    x<-x %>% dplyr::rename("sample_id"="File.ID")
  }
  x<-x %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
  y<-purrr::map(x,function(x) x %>% dplyr::mutate(pr=list(predict(.$M1[[1]]))))
  y<-purrr::map(y,function(x) x %>% dplyr::mutate(Tm = stats::approx(y[[1]]$pr[[1]],y[[1]]$M1[[1]]$model$`x$C`, xout=min(y[[1]]$pr[[1]],na.rm=TRUE)+(0.5*(max(y[[1]]$pr[[1]], na.rm=TRUE)-min(y[[1]]$pr[[1]], na.rm=TRUE))))$y))
  y<-purrr::map(y,function(x)x %>% dplyr::mutate(uniqueID=.$uniqueID[1],
                                                 k_ = .$M1[[1]]$rank,
                                                 treatment=.$treatment[1],
                                                 rss=deviance(.$M1[[1]]),
                                                 CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
                                                 AUC = pracma::trapz(.$M1[[1]]$fitted.values[(which(abs(.$M1[[1]]$fitted.values-0.5)==min(abs(.$M1[[1]]$fitted.values-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fitted.values-0.5)))+1)]),
                                                 rsq=summary(x$M1[[1]])$r.sq,
                                                 n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
                                                 sample_name=.$sample_name[1],
                                                 missing_pct=ifelse(any(names(x)=="missing_pct"),.$missing_pct[1],NA),
                                                 replicate=replicate) %>% 
                  dplyr::ungroup(.))
  
  return(y)
}