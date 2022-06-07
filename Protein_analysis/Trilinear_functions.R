#Trilinear functions
DLR<-function(d){
  #preallocate final result as a list
  d<-lapply(d,function(x) x %>% dplyr::mutate(C=as.numeric(C)))
  
  df_n<-vector(mode = "list", length(d))
  df_n[[1]]<-data.frame()
  df1<-df_n
  df2<-df1
  df3<-df1
  df_1<-df_n
  df0<-df_n
  df_0<-df_n
  
  #d<-lapply(d,function(x) x %>% dplyr::mutate(I=as.numeric(I)))
  
  df_1<-purrr::map(d, function(x) { 
    x %>%
      dplyr::group_by(C) %>%
      dplyr::mutate(I=mean(I,na.rm=TRUE)) %>% 
      dplyr::ungroup(.) %>% 
      distinct(.)
    
  })
  #rank intensity values using 3 regions,  rename column as LineRegion
  LR<-purrr::map(df_1, function(x) { 
    dplyr::ntile(dplyr::desc(x$I),3)%>%
      as.data.frame(.) %>% dplyr::rename("LineRegion"=".")})
  df_1<-purrr::map(df_1,function(x){x %>% dplyr::select(-LineRegion)})#remove Line Region column from one treatment before merging
  
  #Add LR to the list
  df_1<-purrr::map2(df_1,LR, function(x,y) {c(x,y) %>% as.data.frame(.)})
  df_1 <-purrr::map(df_1,function(x){x %>% dplyr::mutate(C = C,I=I,CC=as.factor(CC))})
  
  #separate by Line Regions
  df1<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==1) %>% as.data.frame(.)})
  df2<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==2) %>% as.data.frame(.)})
  df3<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==3) %>% as.data.frame(.)})
  
  #preallocate model data per line region
  LM1<-list(NA)
  LM2<-list(NA)
  LM3<-list(NA)
  df1<-purrr::map(df1,function(x) x[order(x$C),])
  df2<-purrr::map(df2,function(x) x[order(x$C),])
  df3<-purrr::map(df3,function(x) x[order(x$C),])
  # #Flag NA values
  df1<-purrr::map(df1,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df2<-purrr::map(df2,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df3<-purrr::map(df3,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  # #remove NA values
  df1<-purrr::map(df1,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df2<-purrr::map(df2,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df3<-purrr::map(df3,function(x) x %>% dplyr::filter(!is.na(x$I)))
  # #remove empty rows for proteins
  df1<-df1 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  df2<-df2 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  df3<-df3 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  #get common uniqueIDs
  d1<-dplyr::intersect(dplyr::bind_rows(df1)$uniqueID,dplyr::bind_rows(df2)$uniqueID)
  CID<-dplyr::intersect(d1,dplyr::bind_rows(df3)$uniqueID)
  #keep common uniqueIDs
  
  df1<-df1 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  df2<-df2 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  df3<-df3 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  #find fitted curves L3<-purrr::map(df3,function(x) tryCatch(lm(formula = I~C,data = x ,na.action=na.omit), error = function(e){NA}))
  L1<-purrr::map(df1,function(x){ 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM1<-purrr::map2(df1,L1,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L2<-purrr::map(df2,function(x){ 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM2<-purrr::map2(df2,L2,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L3<-purrr::map2(df3,seq(df3),function(x,y) { 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM3<-purrr::map2(df3,L3,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  #linear fit per line region
  LM1<-purrr::map2(df1,L1,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM2<-purrr::map2(df2,L2,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM3<-purrr::map2(df3,L3,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  
  
  #fitted curves
  x1<-purrr::map(LM1, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x2<-purrr::map(LM2, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x3<-purrr::map(LM3, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  
  #fit per line region with confidence intervals
  fit1<-purrr::map(LM1,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  fit2<-purrr::map(LM2,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  fit3<-purrr::map(LM3,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  
  #keep last value for CI 
  fit1 <- purrr::map(fit1,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))
  fit2 <- purrr::map(fit2,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))
  fit3 <- purrr::map(fit3,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))
  
  #append # of fitted curves to original data (columns must have the same rows for map2)
  df1<-purrr::map2(df1,x1,function(x,y)x %>% dplyr::mutate(fitn=y))
  df2<-purrr::map2(df2,x2,function(x,y)x %>% dplyr::mutate(fitn=y))
  df3<-purrr::map2(df3,x3,function(x,y)x %>% dplyr::mutate(fitn=y))
  
  #append # of fitted curves to original data (columns must have the same rows for map2)
  df1<-purrr::map2(df1,fit1,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  df2<-purrr::map2(df2,fit2,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  df3<-purrr::map2(df3,fit3,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  
  #Reassign line Regions if intensity falls within previous Line Region's CI
  
  df2<-purrr::map2(df1,df2,function(x,y)y %>%
                     dplyr::mutate(LineRegion=ifelse(any(y$I<tail(x$I-x$CI,1)),2,1))) 
  
  df3<-purrr::map2(df2,df3,function(x,y)y %>% 
                     dplyr::mutate(LineRegion=ifelse(any(y$I<tail(x$I-x$CI,1)),3,2))) 
  
  df1<-df1 %>% dplyr::bind_rows(.)
  df2<-df2 %>% dplyr::bind_rows(.)
  df3<-df3 %>% dplyr::bind_rows(.)
  #merge all prepared lists to one data frame
  df_0<-rbind(df1,df2,df3) 
  
  #define line Region as a factor
  df_0$LineRegion<-as.factor(df_0$LineRegion)
  df_0<-df_0 %>% dplyr::group_split(uniqueID)
  return(df_0)
}
#this function takes original data with replicates as an input
CP<-function(df_0,d,PSM=FALSE,CARRIER=TRUE){#df_0 is the result data frame and d is the orginal data with replicates
  
  df_0<-df_0 %>% unique(.)
  d<-d %>% unique(.)
  
  #remove duplicate values
  
  
  df_0<-df_0 %>% dplyr::select(-missing,-CI)
  
  #keep the IDs in df_0 which are present in d
  df_0<-df_0 %>% dplyr::filter(uniqueID %in% d$uniqueID)
  
  #Split into lists once again
  d<-d %>% dplyr::group_split(uniqueID) 
  
  #remove IDs that are not common in both datasets
  d<-d %>% purrr::keep(function(x) x$uniqueID[1] %in% df_0$uniqueID)
  #split into list
  df_0<-df_0 %>% dplyr::group_split(uniqueID)
  #For the original data (unlabeled LR) define LR with intensities
  df_0<-suppressWarnings(purrr::map2(d,df_0,function(x,y) x %>% #if the intensity in DF is greater than the max(LR2) label 1 else if the intensity is less than min (LR=2)label 3
                                       dplyr::mutate(LineRegion=as.numeric(ifelse(x$I>=min(y$I[y$LineRegion==1]),1,ifelse(x$I<min(y$I[y$LineRegion==2]),3,2))))))
  
  df_n<-vector(mode = "list", length(df_0))
  ctest<-df_n
  dap<-data.frame()
  Split<-df_n
  #This function is to verify consistent line Region assignments for C (temperature) across replicates
  df_0<-purrr::map(df_0,function(x)x %>% dplyr::arrange(C) %>% dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup())
  
  #subset of the data with shared line region values, using purrr::map to keep size constant
  Split<-purrr::map(df_0,function(x)x %>%subset(n<max(n)) %>% data.frame(.) %>% dplyr::mutate(LineRegion=as.numeric((.$LineRegion))))
  #remove NA values
  Split<-dplyr::bind_rows(Split)
  df_0<-dplyr::bind_rows(df_0)
  names<-dplyr::intersect(names(df_0),names(Split))
  names<-names[which(!names == "LineRegion")]
  if(isTRUE(PSM)){
    
    if(!any(names(df_0)=="missing_pct")){
      if(isTRUE(CARRIER)){
        df.temps<-length(unique(df.temps$temperature))-1
      }else{
        df.temps<-length(unique(df.temps$temperature))
      }
      
      #missing values
      miss_v<-data.frame(NA)
      #max replicates
      roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
        if(length(x) != 1) stop("'x' must be of length 1")
        10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
      }
      
      num<-roundUpNice(length(unique(df_0$C)))
      #append missing value data
      if(length(unique(miss_v$C))==df.temps & length(miss_v$C)>df.temps){
        
      }else{
        df_0$missing_pct<-(100*(num-df.temps)/num)
      }
      
    }
    #Join df_0 with the subset of values
    if(any(names(df_0) =="XCor_l")){
      
      dap<-df_0 %>% dplyr::left_join(Split,by=names)
    }else{
      
      dap<-df_0 %>% dplyr::left_join(Split,by=names)
    }
    dap<-dap %>% dplyr::group_split(uniqueID)
    dap<-purrr::map(dap,function(x) x %>% dplyr::mutate(LineRegion=as.numeric(as.character(x$LineRegion.x)),
                                                        C=as.numeric(as.character(x$C))))
    # dap<-purrr::map(dap,function(x) x %>% dplyr::mutate(C=as.numeric(as.character(x$C))))
    dap<-purrr::map(dap,function(x)x %>% dplyr::mutate(LineRegion=ifelse(x$C<=max(x$C[x$LineRegion==1],na.rm=TRUE),1,ifelse(x$C>=min(x$C[x$LineRegion==3],na.rm=TRUE),3,2))) %>% dplyr::select(-LineRegion.y,-LineRegion.x))
    
    return(dap)
  }else{
    
    #Join df_0 with the subset of values
    dap<-df_0 %>% dplyr::left_join(Split,by=names)
    dap<-dap %>% dplyr::group_split(uniqueID)
    dap<-purrr::map(dap,function(x) x %>% dplyr::mutate(LineRegion=as.numeric(as.character(x$LineRegion.x)),
                                                        C=as.numeric(as.character(x$C))))
    
    dap<-purrr::map(dap,function(x)x %>% dplyr::mutate(LineRegion=ifelse(x$C<=max(x$C[x$LineRegion==1],na.rm=TRUE),1,ifelse(x$C>=min(x$C[x$LineRegion==3],na.rm=TRUE),3,2))) %>% dplyr::select(-LineRegion.y,-LineRegion.x))
    
    
    
    
    return(dap)
  }
  
  
  #Trilinear functions
  tlstat<-function(DF,df,df1,norm=FALSE,Filters=FALSE,Ftest=FALSE,show_results=FALSE){
    i<-1
    #convert to df for numeric variables C and I 
    DF<-dplyr::bind_rows(DF)
    df1<-dplyr::bind_rows(df1)
    df<-dplyr::bind_rows(df)
    #convert factor to numeric columns
    df1$C<-as.numeric(as.vector(df1$C))
    df$C<-as.numeric(as.vector(df$C))
    DF$C<-as.numeric(as.vector(DF$C))
    
    df1$I<-as.numeric(as.vector(df1$I))
    df$I<-as.numeric(as.vector(df$I))
    DF$I<-as.numeric(as.vector(DF$I))
    
    df1$uniqueID<-as.character(df1$uniqueID)
    df$uniqueID<-as.character(df$uniqueID)
    DF$uniqueID<-as.character(DF$uniqueID)
    #convert back to list
    df1<- df1 %>% dplyr::group_split(uniqueID)
    df<-df %>% dplyr::group_split(uniqueID)
    DF<-DF %>% dplyr::group_split(uniqueID)
    
    if(!isTRUE(norm)){
      mean1<-list()
      mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
      mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                           dplyr::group_nest(LineRegion,uniqueID) %>%
                           dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=purrr::map(M1,function(x){deviance(x)}),
                                         Rsq=purrr::map(M1,function(x){summary(x)$r.squared}), 
                                         treatment="vehicle",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
      
      mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      
      
      #define linear models with outputs
      
      mean1_1<-list()
      mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
      
      mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                             dplyr::group_nest(LineRegion,uniqueID) %>% 
                             dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                           CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                           Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                           slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                           intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                           rss=map(M1,function(x){deviance(x)}),
                                           Rsq=map(M1,function(x){summary(x)$r.squared}), 
                                           treatment="treated",
                                           uniqueID=x$uniqueID[1],
                                           n=ifelse(class(M1)=="lm",1,0)))
      
      
      
      mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      # null hypothesis
      #null
      mean3<-list()
      mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
      
      
      mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx( x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}),
                                         treatment="null",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
      
      
      mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      if(isTRUE(show_results)){
        results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(sample_name)
        return(results)
      }
      if (isTRUE(Filters)){
        #Apply lax Rsq and negative slope filter to remove flat melt curves
        mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        
        mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
        mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
        mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      }
      #convert to df and split by uniqueID 
      mean1<-dplyr::bind_rows(mean1)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean3<-dplyr::bind_rows(mean3)
      
      #obtain common uniqueIDs
      CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-intersect(CID,mean3$uniqueID)
      #subset common uniqueIDs
      mean1<-mean1 %>% subset(uniqueID %in% CID)
      mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
      mean3<-mean3  %>% subset(uniqueID %in% CID)
      #split into lists by uniqueID
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
      mean3<-mean3 %>% dplyr::group_split(uniqueID)
      if(isTRUE(Ftest)){
        #Calculate rss0 and rss1 null vs alt
        rss0<-purrr::map(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
        rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
                                                                Tm = y$Tm[[1]]-x$Tm[[1]]))
        #params for null and alternative models
        pN<-purrr::map(mean3,function(x)x %>% dplyr::summarise(pN = 4))
        pA<-purrr::map(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
        
        #sum residuals
        n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = as.numeric(nrow(dplyr::bind_rows(x$data))) + as.numeric(nrow(dplyr::bind_rows(y$data)))))
        #degrees of freedom before
        d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
        d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA)) 
        #delta RSS
        rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS-y$RSS %>% as.data.frame(.))
        #bind rows
        rssDiff<-dplyr::bind_rows(rssDiff)$.
        rss0<-dplyr::bind_rows(rss0)$RSS
        rss1<-dplyr::bind_rows(rss1)
        d2<-dplyr::bind_rows(d2)$d2
        d1<-dplyr::bind_rows(d1)$d1
        #F-test
        Fvals<-(rssDiff/rss1$RSS)*(d2/d1)
        #append results to data
        ResF<-purrr::map2(mean1,Fvals,function(x,y) x %>% dplyr::mutate(Fvals=y))
        ResF<-purrr::map2(ResF,rss0,function(x,y) x %>% dplyr::mutate(rss0=y))
        ResF<-purrr::map2(ResF,rss1,function(x,y) x %>% dplyr::mutate(rss1=y$RSS,Tm=y$Tm))
        ResF<-purrr::map2(ResF,rssDiff,function(x,y) x %>% dplyr::mutate(rssDiff=y))
        ResF<-purrr::map2(ResF,d1,function(x,y) x %>% dplyr::mutate(d1=y))
        ResF<-purrr::map2(ResF,d2,function(x,y) x %>% dplyr::mutate(d2=y))
        
        #convert to df
        mean1<-dplyr::bind_rows(mean1)
        ResF<-dplyr::bind_rows(ResF)
        
        #convert results to list
        testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-treatment)
        testResults<-testResults%>% dplyr::left_join(ResF,by="uniqueID")
        
        #p-val
        testResults<-testResults %>%
          dplyr::mutate(pV = 1-pf(testResults$Fvals,df1=testResults$d1,df2=testResults$d2))
        testResults<-testResults %>% dplyr::mutate(pAdj = p.adjust(.$pV,method="BH"))
        if(isTRUE(Ftest)){
          return(testResults)
        }
        #V is zero, so it would not work as a scaling factor
        ggplot(testResults)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
          geom_line(aes(x=Fvals,y= df(Fvals,df1=4,df2=8)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,10))+
          ggplot2::xlab("F-values")
        #testResults<-testResults %>% dplyr::filter(pAdj<0.05)
        #scale variables
        M<-median(testResults$rssDiff,na.rm=TRUE)
        V<-mad(testResults$rssDiff,na.rm=TRUE)
        #alternative scaling factor sig0-sq
        altScale<-0.5*V/M
        #filter out negative delta rss
        testResults<-testResults %>% dplyr::filter(rssDiff>0)
        #effective degrees of freedom
        ed1<-MASS::fitdistr(x=testResults$rssDiff, densfun = "chi-squared", start = list(df=1))[["estimate"]]
        ed2<-MASS::fitdistr(x=testResults$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
        #scale data
        testScaled <-testResults %>% 
          dplyr::mutate(rssDiff = .$rssDiff/altScale,
                        rss1 =.$rss1/altScale,
                        d1=ed1,
                        d2=ed2)
        #
        #new F-test
        testScaled<-testScaled %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
        Fvals<-testScaled$Fvals
        d1<-testScaled$d1
        d2<-testScaled$d2
        
        #scaled values 
        ggplot(testScaled)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
          geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,10))+
          ggplot2::xlab("F-values")
        #Define checked as filtered protein IDs
        check<-testScaled$uniqueID
        test<-testScaled %>% dplyr::filter(.$pAdj<0.05)
        ggplot(test)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
          geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,10))+
          ggplot2::xlab("F-values")
        
        mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
        mean1_1<-dplyr::bind_rows(mean1_1)
        mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
        mean3<-dplyr::bind_rows(mean3)
        mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      }
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      return(results)
    }else if (isTRUE(norm)){
      mean1<-list()
      mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="vehicle",uniqueID=df[[i]]$uniqueID[1],Tm=rep(0,1))
      mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                           dplyr::group_nest(LineRegion,uniqueID) %>%
                           dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}),
                                         treatment="vehicle",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
      
      
      mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      
      #define linear models with outputs
      
      mean1_1<-list()
      mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
      
      mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                             dplyr::group_nest(LineRegion,uniqueID) %>% 
                             dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                           CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                           Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                           slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                           intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                           rss=map(M1,function(x){deviance(x)}),
                                           Rsq=map(M1,function(x){summary(x)$r.squared}),
                                           treatment="treated",
                                           uniqueID=x$uniqueID[1],
                                           n=ifelse(class(M1)=="lm",1,0)))
      
      
      mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      
      
      
      # null hypothesis
      #null
      mean3<-list()
      mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),treatment="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
      
      
      mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}),
                                         treatment="null",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
      
      
      mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
      #convert to df and split by uniqueID 
      mean1<-dplyr::bind_rows(mean1)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean3<-dplyr::bind_rows(mean3)
      
      #obtain common uniqueIDs
      CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-intersect(CID,mean3$uniqueID)
      #subset common uniqueIDs
      mean1<-mean1 %>% subset(uniqueID %in% CID)
      mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
      mean3<-mean3  %>% subset(uniqueID %in% CID)
      #split into lists by uniqueID
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
      mean3<-mean3 %>% dplyr::group_split(uniqueID)
      
      
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      
      if(isTRUE(Filters)){
        
        #Apply lax Rsq and negative slope filter to remove flat melt curves
        mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
        
        mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
        mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
        mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      }
      #convert to df and split by uniqueID 
      mean1<-dplyr::bind_rows(mean1)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean3<-dplyr::bind_rows(mean3)
      
      #obtain common uniqueIDs
      CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-intersect(CID,mean3$uniqueID)
      #subset common uniqueIDs
      mean1<-mean1 %>% subset(uniqueID %in% CID)
      mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
      mean3<-mean3  %>% subset(uniqueID %in% CID)
      #split into lists by uniqueID
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
      mean3<-mean3 %>% dplyr::group_split(uniqueID)
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      if(isTRUE(Ftest)){
        #Calculate rss0 and rss1 null vs alt
        rss0<-purrr::map(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
        rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
                                                                Tm = y$Tm[[1]]-x$Tm[[1]]))
        #params for null and alternative models
        pN<-purrr::map(mean3,function(x)x %>% dplyr::summarise(pN = 4))
        pA<-purrr::map(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
        #sum residuals
        n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = as.numeric(nrow(dplyr::bind_rows(x$data))) + as.numeric(nrow(dplyr::bind_rows(y$data)))))
        #degrees of freedom before
        d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
        d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA))
        #delta RSS
        rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS-y$RSS %>% data.frame(.))
        #bind rows
        rssDiff<-dplyr::bind_rows(rssDiff)$.
        rss0<-dplyr::bind_rows(rss0)$RSS
        rss1<-dplyr::bind_rows(rss1)
        d2<-dplyr::bind_rows(d2)$d2
        d1<-dplyr::bind_rows(d1)$d1
        #F-test
        Fvals<-(rssDiff/rss1$RSS)*(d2/d1)
        #append results to data
        #append results to data
        ResF<-purrr::map2(mean1,Fvals,function(x,y) x %>% dplyr::mutate(Fvals=y))
        ResF<-purrr::map2(ResF,rss0,function(x,y) x %>% dplyr::mutate(rss0=y))
        ResF<-purrr::map2(ResF,rss1,function(x,y) x %>% dplyr::mutate(rss1=y,Tm=y$Tm))
        ResF<-purrr::map2(ResF,rssDiff,function(x,y) x %>% dplyr::mutate(rssDiff=y))
        ResF<-purrr::map2(ResF,d1,function(x,y) x %>% dplyr::mutate(d1=y))
        ResF<-purrr::map2(ResF,d2,function(x,y) x %>% dplyr::mutate(d2=y))
        
        #convert to df
        mean1<-dplyr::bind_rows(mean1)
        ResF<-dplyr::bind_rows(ResF)
        
        #convert results to list
        testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-treatment)
        testResults<-testResults%>% dplyr::left_join(ResF,by="uniqueID")
        
        
        
        #p-val
        testResults<-testResults %>%
          dplyr::mutate(pV = 1-pf(testResults$Fvals,df1=testResults$d1,df2=testResults$d2))
        testResults<-testResults %>% dplyr::mutate(pAdj = p.adjust(.$pV,method="BH"))
        
        #V is zero, so it would not work as a scaling factor
        ggplot(testResults)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
          geom_line(aes(x=Fvals,y= df(Fvals,df1=4,df2=8)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,100))+
          ggplot2::xlab("F-values")
        #scale variables
        M<-median(testResults$rssDiff,na.rm=TRUE)
        V<-mad(testResults$rssDiff,na.rm=TRUE)
        #alternative scaling factor sig0-sq
        altScale<-0.5*V/M
        #filter out negative delta rss
        testResults<-testResults %>% dplyr::filter(rssDiff>0)
        #effective degrees of freedom
        ed1<-MASS::fitdistr(x=testResults$rssDiff, densfun = "chi-squared", start = list(df=1))[["estimate"]]
        ed2<-MASS::fitdistr(x=testResults$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
        #scale data
        testScaled <-testResults %>%
          dplyr::mutate(rssDiff = .$rssDiff/altScale,
                        rss1 =.$rss1/altScale,
                        d1=ed1,
                        d2=ed2)
        #
        #new F-test
        testScaled<-testScaled %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
        Fvals<-testScaled$Fvals
        d1<-testScaled$d1
        d2<-testScaled$d2
        
        #scaled values
        ggplot(testScaled)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
          geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,10))+
          ggplot2::xlab("F-values")
        #Define checked as filtered protein IDs
        check<-testScaled$uniqueID
        test<-testScaled %>% dplyr::filter(.$pAdj<0.01)
        ggplot(test)+
          geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
          geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
          theme_bw() +
          coord_cartesian(xlim=c(0,10))+
          ggplot2::xlab("F-values")
        
        mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
        mean1_1<-dplyr::bind_rows(mean1_1)
        mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
        mean3<-dplyr::bind_rows(mean3)
        mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      }
      if(isTRUE(Ftest)){
        return(testResults)
      }
      #convert to df and split by uniqueID 
      mean1<-dplyr::bind_rows(mean1)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean3<-dplyr::bind_rows(mean3)
      
      #obtain common uniqueIDs
      CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-intersect(CID,mean3$uniqueID)
      #subset common uniqueIDs
      mean1<-mean1 %>% subset(uniqueID %in% CID)
      mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
      mean3<-mean3  %>% subset(uniqueID %in% CID)
      
      #split into lists by uniqueID
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
      mean3<-mean3 %>% dplyr::group_split(uniqueID)
      #apply Tm data from LR 2
      mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
      mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
      mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
      
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      return(results)
    }
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
    mean3<-mean3  %>% subset(uniqueID %in% CID)
    #split into lists by uniqueID
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    #apply Tm data from LR 2
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    
    return(results)
  }
  tlf<-function(tlresults,DFN,APfilt=TRUE,PF=TRUE){
    ##Apply Filters
    #####################
    if(isTRUE(APfilt)){
      tlresults1<-tlresults#save unfiltered data
      #apply filters prior to hypothesis testing
      tlresults<-tlresults %>% keep(function(x) min(as.numeric(x$Rsq),na.rm=TRUE) >= 0.40)
      tlresults<-tlresults %>% keep(function(x) mean(as.numeric(x$slope),na.rm=TRUE) <= -0.02)
      #tlresults<-tlresults %>% keep(function(x)  sum(data.frame(x)[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss'],na.rm=TRUE) <10)#move data with extremely large RSS values 
      # tlresults<-tlresults %>% keep(function(x) sum(data.frame(x)[!stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss'],na.rm=TRUE) <1.5)
      tlresults<-tlresults %>% keep(function(x) sum(unlist(x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']),na.rm=TRUE) > sum(unlist(x[!stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']),na.rm=TRUE))#remove data with extremely large RSS values 
      tlresults<-tlresults %>% keep(function(x) mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "vehicle"),'Tm']),na.rm=TRUE) < mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "treated"),'Tm']),na.rm=TRUE))
      #tlresults<-tlresults %>% keep(function(x) max(data.frame(x)$slope[x$LineRegion==2],na.rm=TRUE) < -0.03)#the linear region have the largest slope < 0.03
      #tlresults<-tlresults %>% keep(function(x) length(x$slope)>8)#remove list values with less than 5 rows
      #tlresults<-tlresults %>% keep(function(x) abs(max(x$slope[!x$LineRegion==2] ,na.rm=TRUE)) < 0.1)#eeps plateau values where the min abs(slope) < 0.06
      #steepest slope in vehicle and treatment has to be less than 0.06C
    }
    Nsum<-list()
    Nsum[[1]]<-data.frame(RSS=0,Tm=0)
    
    tlresults<-purrr::map(tlresults,function(x) x %>% dplyr::mutate(sample_name=data[[1]]$sample_name[1]))
    if(any(class(tlresults)=="data.frame")){
      tlresults<-tlresults %>% dplyr::group_split(sample_name)
    }
    #get the summed rss values for null
    Nsum<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(treatment), pattern = "null")) %>% 
                       dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,treatment,uniqueID)%>% head(.,1))
    
    #get the summed rss values for vehicle
    Rssv<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(treatment), pattern = "vehicle")) %>% 
                       dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,treatment,uniqueID)%>%head(.,1))
    #get the summed rss values for treated
    Rsst<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(treatment), pattern = "treated")) %>% 
                       dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,treatment,uniqueID)%>% head(.,1))
    #find the rss difference between treated and vehicle 
    
    Rssv<-purrr::map(Rssv,function(x)na.omit(x))
    Rsst<-purrr::map(Rsst,function(x)na.omit(x))
    #find common IDs
    CID<-intersect(dplyr::bind_rows(Rsst)$uniqueID,dplyr::bind_rows(Rssv)$uniqueID)
    #keep common IDs
    Rssv<-Rssv %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID)) 
    Rsst<-Rsst %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))
    Nsum<-Nsum %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))                           
    K1<-data.frame(dplyr::bind_rows(purrr::map2(Rsst,Rssv,function(x,y) data.frame(RSSd = x$RSS-y$RSS, Tma = x$Tm[1] - y$Tm[1])))) 
    K2<-data.frame(uniqueID = dplyr::bind_rows(Rssv)$uniqueID)
    Dsum<-data.frame(K1,K2)
    
    Dsum$rank<- dplyr::ntile(Dsum$Tma,7)
    #keep data where the difference in RSS is less than the null
    #nsum converted to data frame
    Nsum<-data.frame(RSSn=dplyr::bind_rows(Nsum))
    names(Nsum)<-c("RSSn","Tmn","treatment","uniqueID")
    Nsum<-Nsum %>% dplyr::filter(uniqueID %in% CID)
    Nsum<-Nsum %>% dplyr::mutate(id=rownames(Nsum))
    
    Nsum$treatment<-as.factor(Nsum$treatment)
    #mutate data frame
    #join two data frames by uniqueID
    Dsum1<-Dsum %>% dplyr::left_join(Nsum,by = c("uniqueID"="uniqueID"))
    #Childs
    Dsum2<-Dsum %>% dplyr::right_join(Nsum,by = c("uniqueID"="uniqueID"))
    
    Dsum<-Dsum2
    Dsum$RSSd<-Dsum1$RSSd
    Dsum$Tma<-Dsum1$Tma
    Dsum<-Dsum %>% dplyr::mutate(rank = dplyr::ntile(Dsum$Tma,7))
    if (isTRUE(PF)){
      #rank the data by Tm change
      
      #arrange data from greater Tm and RSS difference to lowest
      Dsum<-dplyr::arrange(Dsum, dplyr::desc(Tma), dplyr::desc(RSSd))  %>% dplyr::filter(RSSd>0) 
      
      test<-data.frame()
      test<-Dsum[which(Dsum$RSSn>Dsum$RSSd),] %>% data.frame()#get the stable proteins (+ = Rsstreated-Rssvehicle)
      rssdec<-data.frame()
      rssdec<-data.frame(data.table::fsort(test$RSSd,decreasing=TRUE))#decreasing Rss differences
      names(rssdec)<-"Rssd"
      
      tmdec<-data.frame()
      tmdec<-data.table::fsort(test$Tma,decreasing=TRUE) %>% data.frame()
      names(tmdec)<-"Tm"
      
      test<-tmdec %>% inner_join(test,by=c("Tm"="Tma"))#orders data by decreasing Tm
      orows<-data.frame()
      orows <- test
      orows$id<-sapply(orows$id, function(x) as.numeric(as.character(x)))
      Df1<-tlresults[orows$id] #divide 1=highly destabilized,4=noeffect,7=highly stabilized
    }else{
      tlresults<-dplyr::bind_rows(tlresults)
      
      #order by RSS differences while keeping original rownames for index
      #create an external data frame for stabilized proteins
      Df1<-Dsum %>% dplyr::left_join(tlresults,by=c("uniqueID")) %>% as.data.frame(.) %>% dplyr::rename("treatment"="treatment.y") %>% 
        dplyr::select(-treatment.x)
      Df1<-Df1 %>% dplyr::group_split(uniqueID)
      
    }
    
    df1<-list()
    #get uniqueID and treatment for stable proteins with decreasing RSS differences
    df1<-purrr::map(Df1,function(x) x %>% dplyr::select(uniqueID,treatment) %>% head(.,1))
    df1<-data.frame(dplyr::bind_rows(df1))
    
    #unlist to data.frame
    #order the original data by RSS differences
    #
    DFN<- dplyr::bind_rows(DFN)
    DFN$uniqueID<-as.vector(DFN$uniqueID)
    df1$uniqueID<-as.vector(df1$uniqueID)
    
    
    df2<-df1 %>% dplyr::right_join(DFN,by=c("uniqueID")) %>% dplyr::rename("treatment"="treatment.y") %>% 
      dplyr::select(-treatment.x)
    
    
    return(list(df1,df2,Df1))
  }
  tlCI<-function(i,df1,df2,Df1,overlay=TRUE,residuals=FALSE,df.temps,PSMs,CARRIER=TRUE){
    roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
      if(length(x) != 1) stop("'x' must be of length 1")
      10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
    }
    null<-data.frame()
    i<-i
    df1<-df1
    df2<-df2[!is.na(df2$I),]
    Df1<-Df1[[i]]
    DF1<-data.frame(NA)
    DF1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) 
    DF1<-DF1[!is.na(DF1$I),]
    
    null<-Df1 %>% subset(treatment == "null")
    
    pred1<-predict(null$M1[[1]], interval="confidence") %>% as.data.frame(.)
    
    if(nrow(null)==2|nrow(null)==3){
      pred2<-predict(null$M1[[2]], interval="confidence")%>% as.data.frame(.)
      pred2<-na.omit(pred2)
      
    }else{
      pred2<-data.frame()
    }
    if(nrow(null)==3){
      pred3<-predict(null$M1[[3]], interval="confidence")%>% as.data.frame(.)
      pred3<-na.omit(pred3)
      
    }else{
      pred3<-data.frame()
    }
    Pred1<-NA
    pred1<-na.omit(pred1)
    
    
    
    FIT<- NA
    LOW<-NA
    HI<-NA
    if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
      Pred<-data.frame(rbind(pred1,pred2,pred3))
    } else if (nrow(pred2)>0 & nrow(pred3)>0){
      Pred<-data.frame(rbind(pred2,pred3))
    } else if (nrow(pred1)>0 & nrow(pred2)>0){
      Pred<-data.frame(rbind(pred1,pred2))  
    }else if (nrow(pred1)>0 & nrow(pred3)>0){
      Pred<-data.frame(rbind(pred1,pred3))
    }else if(nrow(pred1)>0){
      Pred<-data.frame(pred1)
    }
    rownames(Pred)<-as.vector(1:nrow(Pred))
    
    #Pred<-Pred[1:length(DF1$C),]##############
    Pred<-cbind(Pred,DF1$C[1:nrow(Pred)],DF1$I[1:nrow(Pred)])################
    names(Pred)<-c("fit","lower","upper","C","I")
    
    Pred$Treatment<-null$treatment[1]##################
    Pred<-na.omit(Pred)
    Pred$C<-as.numeric(as.vector(Pred$C))
    Pred$I<-as.numeric(as.vector(Pred$I))
    PLN<-ggplot2::ggplot(Pred, ggplot2::aes(x = C,y = I,color=Treatment)) +
      ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+
      ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ 
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
      annotate("text", x=-0.15, y=min(Pred$I),label=paste("RSS= ",round(sum(unlist(null$rss)),3)))+
      annotate("text",
               x = Pred[which(round(Pred$fit,1)==0.5)[1],]$C,
               y = -0.15,
               label=paste(Pred[which(round(Pred$fit,1)==0.5)[1],]$C),
               colour="red"
      )+theme(legend.position="bottom")
    
    
    DF_f<-df2 %>%subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(treatment=ifelse(CC==0,'vehicle','treated')) %>% subset(treatment=="vehicle")
    
    vehicle<-Df1 %>% subset(treatment == "vehicle")
    
    pred1<-predict(vehicle$M1[[1]], interval="confidence")%>% as.data.frame(.)
    
    if(nrow(vehicle)==2|nrow(vehicle)==3){
      pred2<-predict(vehicle$M1[[2]], interval="confidence")%>% as.data.frame(.)
      pred2<-na.omit(pred2)
      
    }else{
      pred2<-data.frame()
    }
    if(nrow(vehicle)==3){
      pred3<-predict(vehicle$M1[[3]], interval="confidence")%>% as.data.frame(.)
      pred3<-na.omit(pred3)
      
    }else{
      pred3<-data.frame()
    }
    Pred1<-NA
    pred1<-na.omit(pred1)
    
    
    
    FIT<- NA
    LOW<-NA
    HI<-NA
    if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
      Pred1<-data.frame(rbind(pred1,pred2,pred3))
    } else if (nrow(pred2)>0 & nrow(pred3)>0){
      Pred1<-data.frame(rbind(pred2,pred3))
    } else if (nrow(pred1)>0 & nrow(pred2)>0){
      Pred1<-data.frame(rbind(pred1,pred2))  
    }else if (nrow(pred1)>0 & nrow(pred3)>0){
      Pred1<-data.frame(rbind(pred1,pred3))
    }else if(nrow(pred1)>0){
      Pred1<-data.frame(pred1)
    }
    
    #Pred<-Pred[1:length(DF1$C),]##############
    Pred1<-data.frame(Pred1,DF_f$C[1:nrow(Pred1)],DF_f$I[1:nrow(Pred1)])################
    names(Pred1)<-c("fit","lower","upper","C","I")
    
    Pred1$Treatment<-vehicle$treatment[1]##################
    Pred1<-na.omit(Pred1)
    rownames(Pred1)<-1:nrow(Pred1)
    Pred1$C<-as.numeric(as.vector(Pred1$C))
    Pred1$I<-as.numeric(as.vector(Pred1$I))
    
    
    DF_f1<-data.frame()
    DF_f1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(treatment=ifelse(CC==0,'vehicle','treated'))
    if(length(unique(DF_f1$treatment))==1){
      DF_f1<-DF_f1 
    }else{
      DF_f1<-DF_f1 %>% subset(treatment =="treated")
    }
    
    treated<-data.frame()
    treated<-Df1 %>% subset(treatment == "treated")
    
    pred1<-predict(treated$M1[[1]], interval="confidence")
    if(nrow(treated)==2|nrow(treated)==3){
      pred2<-predict(treated$M1[[2]], interval="confidence")%>% as.data.frame(.)
      pred2<-na.omit(pred2)
    }else{
      pred2<-data.frame()
    }
    if(nrow(treated)==3){
      pred3<-predict(treated$M1[[3]], interval="confidence")%>% as.data.frame(.)
      pred3<-na.omit(pred3)
    }else{
      pred3<-data.frame()
    }
    
    pred1<-na.omit(pred1)
    
    
    Pred2<-NA
    FIT<- NA
    LOW<-NA
    HI<-NA
    if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
      Pred2<-data.frame(rbind(pred1,pred2,pred3))
    } else if (nrow(pred2)>0 & nrow(pred3)>0){
      Pred2<-data.frame(rbind(pred2,pred3))
    } else if (nrow(pred1)>0 & nrow(pred2)>0){
      Pred2<-data.frame(rbind(pred1,pred2))  
    }else if (nrow(pred1)>0 & nrow(pred3)>0){
      Pred2<-data.frame(rbind(pred1,pred3))
    }else if(nrow(pred1)>0){
      Pred2<-data.frame(pred1)
    }
    rownames(Pred2)<-as.vector(1:nrow(Pred2))
    
    #Pred<-Pred[1:length(DF1$C),]##############
    Pred2<-data.frame(Pred2,DF_f1$C[1:nrow(Pred2)],DF_f1$I[1:nrow(Pred2)])################
    names(Pred2)<-c("fit","lower","upper","C","I")
    
    Pred2$Treatment<-treated$treatment[1]##################
    Pred2<-na.omit(Pred2)
    rownames(Pred2)<-as.vector(1:nrow(Pred2))
    #Area under the curve using trapezoid rule
    
    P1_AUC <- pracma::trapz(Pred1$I)
    P2_AUC <- pracma::trapz(Pred2$I)
    #Residuals
    rn<-data.frame(residuals(null$M1[[1]]))
    
    if(nrow(null)>2){
      rn<-data.frame(residuals=c(residuals(null$M1[[1]]),residuals(null$M1[[2]]),residuals(null$M1[[3]])))
    }else if (nrow(null)>1){
      rn<-data.frame(residuals=c(residuals(null$M1[[1]]),residuals(null$M1[[2]])))
    }else if (nrow(null)==1){
      rn<-data.frame(residuals=c(residuals(null$M1[[1]])))
    }
    Pred<-cbind(Pred,rn[1:nrow(Pred),])
    names(Pred)<-c("fit","lower","upper","C","I","Treatment",'residuals')
    rn<-data.frame(residuals(vehicle$M1[[1]]))
    if(nrow(vehicle)==3){
      rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]]),residuals(vehicle$M1[[3]])))
    }else if (nrow(vehicle)>1){
      rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]])))
    }else if (nrow(vehicle)==1){
      rn<-data.frame(residuals=c(residuals(vehicle$M1[[1]])))
    }
    Pred1<-cbind(Pred1,rn[1:nrow(Pred1),])
    names(Pred1)<- c("fit","lower","upper","C","I","Treatment",'residuals')
    
    Pred1$uniqueID<-vehicle$uniqueID[1]
    rn<-data.frame(residuals(treated$M1[[1]]))
    if(nrow(treated)==3){
      rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]]),residuals(treated$M1[[3]])))
    }else if (nrow(treated)>2){
      rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]])))
    }else if (nrow(treated)==1){
      rn<-data.frame(residuals=c(residuals(treated$M1[[1]])))
    }
    
    Pred2<-cbind(Pred2,rn[1:nrow(Pred2),])
    names(Pred2)<-c("fit","lower","upper","C","I","Treatment",'residuals')
    
    Pred2$uniqueID<-treated$uniqueID[1]
    Preds<-rbind(Pred1,Pred2)
    Preds$C<-as.numeric(as.vector(Preds$C))
    Preds$I<-as.numeric(as.vector(Preds$I))
    Pred2$C<-as.numeric(as.vector(Pred2$C))
    Pred2$I<-as.numeric(as.vector(Pred2$I))
    DF1$treatment<-as.factor(DF1$treatment)
    
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = residuals,color=Treatment)) +
      ggplot2::geom_point()+ ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+
      ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    if(isTRUE(residuals)){
      print(PLrs)
    }
    Tm_d<-round(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),1)
    
    #If dealing with peptides are true
    if(isTRUE(PSMs)){
      Pred1<-dplyr::bind_rows(Pred1) %>% unique(.)
      Pred2<-dplyr::bind_rows(Pred2) %>% unique(.)
      
      Pred1<-Pred1 %>% dplyr::group_split(uniqueID,C)
      Pred1<-purrr::map2(Pred1,seq(Pred1),function(x,y) x %>% dplyr::mutate(replicate=y))
      Pred2<-Pred2 %>% dplyr::group_split(uniqueID,C)
      Pred2<-purrr::map2(Pred2,seq(Pred2),function(x,y) x %>% dplyr::mutate(replicate=y))
      
      Pred1<-dplyr::bind_rows(Pred1)
      Pred2<-dplyr::bind_rows(Pred2)
      PSMs_v<-max(as.numeric(Pred1$replicate),na.rm=TRUE)
      PSMs_t<-max(as.numeric(Pred2$replicate),na.rm=TRUE)
      
      
      PLR_P1<-ggplot2::ggplot(Pred1, ggplot2::aes(x = C,y = fit,color=Treatment))+
        ggplot2::geom_point(Pred1, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
        ggplot2::geom_ribbon(data=Pred1,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
        annotate("segment", x = min(Pred1$C), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "blue",linetype=2)+
        annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "blue",linetype=2)+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+theme(legend.position="bottom")+
        annotate("text",
                 x = 43,
                 y = 1.2,
                 label=paste0("Peptides = ", PSMs_v),
                 colour="#00BFC4",
                 size=3.5
        )+ggplot2::annotate("text", x=43, y=-0.35, label= paste("\u03A3","RSS = ", round(sum(unlist(treated$rss),unlist(vehicle$rss)),3)),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.45, label=  paste("\u0394", "AUC = ",abs(round(P2_AUC-P1_AUC,3))),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.55, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"),size=3.5)+
        annotate("text",
                 x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="blue",
                 size=3.5)
      
      
      
      PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
        ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ylim(-0.6,1.6)+xlim(37,68)+theme(legend.position="bottom")+
        annotate("text",
                 x = 43,
                 y = 1.1,
                 label=paste0("Peptides = ", PSMs_t),
                 colour="#F8766D",
                 size=3.5
        )+
        annotate("text",
                 x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="red",
                 size=3.5
        )
      return(PLR_P2)
      if(overlay=="TRUE"){
        AUCd<-abs(round(P2_AUC-P1_AUC,2))
        p<-expression(paste(Delta, "AUCdiff"))
        AUCd<-as.numeric(AUCd)
        miss_v<-data.frame(NA)
        miss_t<-data.frame(NA)
        miss_v<-DF1%>% dplyr::filter(treatment=="vehicle")
        miss_t<-DF1 %>% dplyr::filter(treatment=="treated")
        #get unique TMT channels
        if(isTRUE(CARRIER)){
          df.temps<-length(unique(df.temps$temperature))-1
        }else{
          df.temps<-length(unique(df.temps$temperature))
        }
        num<-max(roundUpNice(length(unique(miss_v$C))),roundUpNice(length(unique(miss_t$C))))
        
        if(length(unique(miss_v$C))==df.temps & length(miss_v$C)>df.temps){
          getmode <- function(v) {
            uniqv <- unique(v)
            uniqv[which.max(tabulate(match(v, uniqv)))]
          }
          Pred1$missing_v<-rep(getmode(miss_v$missing_pct),nrow(Pred1))
        }else{
          
          Pred1$missing_v<-rep((100*(num-df.temps)/num),nrow(Pred1))
        }
        if(length(unique(miss_t$C))==df.temps & length(miss_t$C)>df.temps){
          getmode <- function(v) {
            uniqv <- unique(v)
            uniqv[which.max(tabulate(match(v, uniqv)))]
          }
          Pred2$missing_t<-rep(getmode(miss_t$missing_pct),nrow(Pred2))
        }else{
          Pred2$missing_t<-rep((100*(num-df.temps)/num),nrow(Pred2))
        }
        
        Pred1$missing_v<-round(Pred1$missing_v,0)
        Pred2$missing_t<-round(Pred2$missing_t,0)
        
        PSMs_v<-max(as.numeric(Pred1$replicate),na.rm=TRUE)
        PSMs_t<-max(as.numeric(Pred2$replicate),na.rm=TRUE)
        
        
        PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
          ggplot2::geom_line(data=Pred2,ggplot2::aes(color=Annotated_Sequence))+
          ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
          ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+
          ggplot2::annotate("text", x=43, y= -0.15, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$treatment), pattern = "vehicle"),'rss']))+
                                                                                          sum(unlist(Df1[stringr::str_detect(tolower(Df1$treatment), pattern = "treated"),'rss'])),3)),size=3.5)+
          ggplot2::annotate("text", x=43, y= -0.25, label=  paste("\u0394", "AUC = ",AUCd),size=3.5)+ 
          ggplot2::annotate("text", x=43, y= -0.35, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"),size=3.5)+ 
          ggplot2::annotate("text", x=43, y= -0.45, label= paste("missing  ",Pred1$missing_v[1],"%"),colour="#00BFC4",size=3.5)+ 
          ggplot2::annotate("text", x=43, y= -0.55, label= paste("missing  ",Pred2$missing_t[1],"%"),colour="#F8766D",size=3.5)+
          annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ylim(-0.8,1.6)+xlim(37,68)+theme(legend.position="bottom")+
          annotate("text",
                   x = max(50,na.rm=TRUE),
                   y = max(1.2,na.rm=TRUE),
                   label=paste0("PSMs = ", PSMs_t),
                   colour="#F8766D",
                   size=3.5
          )
        
        par(mfrow=c(2,2))
        return(PLR_P2)
      }
    }else{
      Pred1<-Pred1 %>% dplyr::group_split(uniqueID,C)
      Pred1<-lapply(Pred1,function(x) x %>% dplyr::mutate(replicate=row.names(.)))
      Pred2<-Pred2 %>% dplyr::group_split(uniqueID,C)
      Pred2<-lapply(Pred2,function(x) x %>% dplyr::mutate(replicate=row.names(.)))
      
      Pred1<-dplyr::bind_rows(Pred1)
      Pred2<-dplyr::bind_rows(Pred2)
      PLR_P1<-ggplot2::ggplot(Pred1, ggplot2::aes(x = C,y = fit,color=Treatment))+
        ggplot2::geom_point(Pred1, mapping=ggplot2::aes(x = C,y = I,color=Treatment,shape=replicate)) +
        ggplot2::geom_ribbon(data=Pred1,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
        annotate("text",
                 x = 2+round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="blue",
                 size=3.5
        )+
        annotate("segment", x = min(Pred1$C), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "blue",linetype=2)+
        annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "blue",linetype=2)+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+ylim(-0.6,1.6)+xlim(37,68)+theme(legend.position="bottom")
      
      
      PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment,shape=replicate)) +
        ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        annotate("text",
                 x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="red",
                 size=3.5
        )+
        annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ylim(-0.6,1.6)+xlim(37,68)+theme(legend.position="bottom")
      if(overlay=="TRUE"){
        AUCd<-abs(round(P2_AUC-P1_AUC,2))
        p<-expression(paste(Delta, "AUCdiff"))
        
        AUCd<-as.numeric(AUCd)
        miss_v<-data.frame(NA)
        miss_t<-data.frame(NA)
        miss_v<-DF1%>% dplyr::filter(treatment=="vehicle")
        miss_t<-DF1 %>% dplyr::filter(treatment=="treated")
        #get unique TMT channels
        df.temps<-length(unique(df.temps$temp_ref))
        num<-max(roundUpNice(length(miss_v$C)),roundUpNice(length(miss_t$C)))
        
        if(length(unique(miss_v$C))==df.temps & length(miss_v$C)>df.temps){
          Pred1$missing_v<-rep(max(miss_v$missing_pct,na.rm=TRUE),nrow(Pred1))
        }else{
          Pred1$missing_v<-rep((100*(num-df.temps)/num),nrow(Pred1))
        }
        if(length(unique(miss_t$C))==df.temps & length(miss_t$C)>df.temps){
          Pred2$missing_t<-rep(max(miss_t$missing_pct,na.rm=TRUE),nrow(Pred2))
        }else{
          Pred2$missing_t<-rep((100*(num-df.temps)/num),nrow(Pred2))
        }
        
        Pred1$missing_v<-round(Pred1$missing_v,0)
        Pred2$missing_t<-round(Pred2$missing_t,0)
        #add shape parameters
        #group_split to generate replicates
        Pred1<-Pred1 %>% dplyr::group_split(C)
        Pred2<-Pred2 %>% dplyr::group_split(C)
        #mutate
        Pred1 <-lapply(Pred1,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        Pred2 <-lapply(Pred2,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        #bind rows
        Pred1<-dplyr::bind_rows(Pred1)
        Pred2<-dplyr::bind_rows(Pred2)
        
        
        PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment,shape=replicate)) +
          ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
          ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+
          ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.15, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$treatment), pattern = "vehicle"),'rss']))+
                                                                                                      sum(unlist(Df1[stringr::str_detect(tolower(Df1$treatment), pattern = "treated"),'rss'])),3)),size=3.5)+
          ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.25, label=  paste("\u0394", "AUC = ",AUCd),size=3.5)+ 
          ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.35, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"),size=3.5)+ 
          ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.45, label= paste("missing  ",Pred1$missing_v[1],"%"),colour="#00BFC4",size=3.5)+ 
          ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.55, label= paste("missing  ",Pred2$missing_t[1],"%"),colour="#F8766D",size=3.5)+
          annotate("text",
                   x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.15,
                   label=paste0(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ylim(-0.8,1.6)+xlim(37,68)+theme(legend.position="bottom")
        
        par(mfrow=c(2,2))
        
        if(isTRUE(PSMs)){
          if(all(Pred1$I==Pred2$I)){
            PLR_P1<-PLR_P1+ guides(shape=guide_legend(title="PSMs"))
            return(PLR_P1)
          }else{
            PLR_P2<-PLR_P2+ guides(shape=guide_legend(title="PSMs"))
            return(PLR_P2)
          }
        }else{
          return(PLR_P2)
        }
      }else if(overlay=="FALSE"){
        miss_v<-data.frame(NA)
        miss_t<-data.frame(NA)
        miss_v<-DF1%>% dplyr::filter(treatment=="vehicle")
        miss_t<-DF1 %>% dplyr::filter(treatment=="treated")
        #get unique TMT channels
        if(isTRUE(CARRIER)){
          df.temps<-length(unique(df.temps$temperature))-1
        }else{
          df.temps<-length(unique(df.temps$temperature))
        }
        num<-max(roundUpNice(length(unique(miss_v$C))),roundUpNice(length(unique(miss_t$C))))
        
        if(length(unique(miss_v$C))==df.temps & length(miss_v$C)>df.temps){
          Pred1$missing_v<-rep(max(miss_v$missing_pct,na.rm=TRUE),nrow(Pred1))
        }else{
          Pred1$missing_v<-rep((100*(num-df.temps)/num),nrow(Pred1))
        }
        if(length(unique(miss_t$C))==df.temps & length(unique(miss_t$C))>df.temps){
          Pred2$missing_t<-rep(max(miss_t$missing_pct,na.rm=TRUE),nrow(Pred2))
        }else{
          Pred2$missing_t<-rep((100*(num-df.temps)/num),nrow(Pred2))
        }
        
        Pred1$missing_v<-round(Pred1$missing_v,0)
        Pred2$missing_t<-round(Pred2$missing_t,0)
        #add shape parameters
        #group_split to generate replicates
        Pred1<-Pred1 %>% dplyr::group_split(C)
        Pred2<-Pred2 %>% dplyr::group_split(C)
        #mutate
        Pred1 <-lapply(Pred1,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        Pred2 <-lapply(Pred2,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        #bind rows
        Pred1<-dplyr::bind_rows(Pred1)
        Pred2<-dplyr::bind_rows(Pred2)
        #add shape parameters
        #group_split to generate replicates
        Pred1<-Pred1 %>% dplyr::group_split(C)
        Pred2<-Pred2 %>% dplyr::group_split(C)
        #mutate
        Pred1 <-lapply(Pred1,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        Pred2 <-lapply(Pred2,function(x) x %>% dplyr::mutate(replicate=as.factor(row.names(.))))
        #bind rows
        Pred1<-dplyr::bind_rows(Pred1)
        Pred2<-dplyr::bind_rows(Pred2)
        PLR<-PLR_P2+
          facet_wrap("Treatment") + 
          ggplot2::annotate("text", x=43, y=-0.55, label= paste("missing % v",Pred2$missing_v[1]))+ 
          ggplot2::annotate("text", x=43, y=-0.65, label= paste("missing % t",Pred2$missing_t[1]))+
          annotate("text",
                   x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.15,
                   label=paste(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red"
          )+
          annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)
        
        if(isTRUE(PSMs)){
          if(all(Pred1$I==Pred2$I)){
            PLR_P1<-PLR_P1+ guides(shape=guide_legend(title="PSMs"))
            print(PLR_P1)
          }else{
            PLR_P2<-PLR_P2+ guides(shape=guide_legend(title="PSMs"))
            print(PLR_P2)
          }
        }else{
          print(PLR_P2)
        }
        
        if(bootstrap==TRUE){
          set.seed(233)
          n<-length(Pred$I)
          mean_orig_v<-mean(Pred$I)
          mean_orig_t<-mean(Pred2$I)
          N<-1000
          for (i in 1:N){
            bs_vehicle[i]<-sample(Pred$I,n,replace=TRUE)
            mean_v[i]<-mean(bs_vehicle[i])
            
            bs_treated[i]<-sample(Pred2$I,n,replace=TRUE)
            mean_t[i]<-mean(bs_treated[i])
          }
          mean_bs_v<-mean(mean_v)
          mean_bs_t<-mean(mean_t)
          bias_v<-mean_orig_v-mean_bs_v
          bias_t<-mean_orig_t-mean_bs_t
          
          summary<-data.frame(uniqueID=rep(Pred$uniqueID[1],2),
                              treatment=c(Pred$treatment[1],Pred2$treatment[1]),
                              Orig_mean=c(mean_orig_v,mean_orig_t),
                              BS_mean=c(mean_bs_v,mean_bs_t),
                              bias_bs=c(bias_v,bias_t),
                              sd_bs_v=sd(mean_v),
                              sd_bs_t=sd(mean_t),
                              CI_2.5_v=quantile(mean_v,0.025),
                              CI_97.5_v=quantile(mean_v,0.975),
                              CI_2.5_t=quantile(mean_t,0.025),
                              CI_97.5_t=quantile(mean_t,0.975))
          
          
        }
      }
    }
  }
  
  
}