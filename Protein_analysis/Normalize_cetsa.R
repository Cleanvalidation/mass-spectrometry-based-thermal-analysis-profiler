#' normalize CETSA data
#'
#' Normalize data according to Pelago paper.
#' Normalizaion works by establishing a subset of curves which pass a set of criteria
#' (see below) and then calculating the median fold change through this set.  A curve is
#' then fit through these points and used to generate correction factors.
#'
#' Curve selection criteria:
#'   \itemize{
#'   \item All 10 temperature points present
#'   \item T7/T1 between 0.4 and 0.6
#'   \item T9/T1 < 0.3
#'   \item T10/T1 < 0.2
#' }
#'
#' @param df. A cleaned CETSA data frame
#' @param temperatures. A list of experimental temperatures
#'
#' @return normalized CETSA data frame
#'
#' @import dplyr
#' @import tidyr
#' @import nls2
#' @export
normalize_cetsa <- function(df, temperatures,Peptide=FALSE,filters=FALSE,CARRIER=TRUE,baseline="min") {
  if(isTRUE(CARRIER)){
    df<-df %>% dplyr::filter(!temperature=="68")
    temperatures<-temperatures[temperatures$temperature<68,]
  }
  if(any(names(df)=="value")&!any(names(df)=="I")){
    df<-df %>% dplyr::mutate(I=value)
  }
  if(!any(names(df)=="sample_id")&any(names(df)=="sample_id")){
    df<-df %>% dplyr::mutate(sample_id=sample_id)
  }
  if(isTRUE(Peptide)){
    if(any(names(df)=="uniqueID")){
      df<-df %>% dplyr::rename("Accession"="uniqueID","value"="I")
    }
    
    df$Accession<-as.factor(df$Accession)
    df<-data.frame(df)
    
    #rank by intensity at the lowest temperature channel
    df_filt<-df %>% dplyr::filter(temp_ref==temp_ref[temperature=min(df$temperature,na.rm=TRUE)]) %>%
      dplyr::mutate(rank=dplyr::ntile(dplyr::desc(I),3)) %>%
      dplyr::select(Accession,Annotated_Sequence,I,sample_id,treatment,sample_name,rank)
    
    df_filt$rank<-as.factor(df_filt$rank)
    
    #select top6000 top12000 and all peptides according to rank
    df_filt3 <- df_filt %>%
      dplyr::filter(!is.na(I)) %>%  
      dplyr::arrange(dplyr::desc(I)) %>% 
      group_by(rank) %>% 
      dplyr::slice(1:6000)
    df_filt5 <- df_filt %>%
      dplyr::filter(!is.na(I)) %>%  
      arrange(I) %>% 
      group_by(rank) %>%
      dplyr::slice(1:12000)
    df_filt10 <- df_filt %>%
      dplyr::filter(!is.na(I)) %>%  
      dplyr::arrange(I) %>% 
      dplyr::group_by(rank) 
    
    df_filt3<-df_filt3 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
    df_filt5<-df_filt5 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
    df_filt10<-df_filt10 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
    #preserve original treatment
    df1<-df
    df<-df %>% dplyr::group_by(sample_id)
    
    name<-dplyr::intersect(names(df1),names(df_filt3))
    #Only keep curves with the topN values
    df3<-df1 %>% dplyr::right_join(df_filt3,by=name)
    df5<-df1 %>% dplyr::right_join(df_filt5,by=name)
    df10<-df1 %>% dplyr::right_join(df_filt10,by=name)
    #remove missing values 
    df3<-df3[!is.na(df3$I),]
    df5<-df5[!is.na(df5$I),]
    df10<-df10[!is.na(df10$I),]
    
    
    df3<-dplyr::bind_rows(df3) %>%
      dplyr::group_split(Accession,treatment,sample_id,sample_name)
    df5<-dplyr::bind_rows(df5) %>%
      dplyr::group_split(Accession,treatment,sample_id,sample_name)
    df10<-dplyr::bind_rows(df10) %>%
      dplyr::group_split(Accession,treatment,sample_id,sample_name)
    
    #Calculate fold changes acros temperatures 7, 9 and 10
    if(order(df3[[1]]$temperature)[1]==length(df3[[1]]$temperature)){
      
      #if the temperatures are reversed, relabel 7,9th and 10th temperature channels
      df.jointP3 <- suppressWarnings(purrr::map(df3,function(x) FC_to_ref(x,baseline)))
      
      df.jointP5 <- suppressWarnings(purrr::map(df5,function(x) FC_to_ref(x,baseline)))
      
      df.jointP10 <- suppressWarnings(purrr::map(df10,function(x) function(x) FC_to_ref(x,baseline)))
    }else{
      df.jointP3 <- suppressWarnings(purrr::map(df3,function(x) FC_to_ref(x,baseline)))
      
      df.jointP5 <- suppressWarnings(purrr::map(df5,function(x) function(x) FC_to_ref(x,baseline)))
      
      df.jointP10 <- suppressWarnings(purrr::map(df10,function(x) function(x) FC_to_ref(x,baseline)))
      
      
    }
    
    
    df.jointP3<- dplyr::bind_rows(df.jointP3)%>%
      dplyr::group_split(treatment,sample_name)
    
    df.jointP5<- dplyr::bind_rows(df.jointP5)%>%
      dplyr::group_split(treatment,sample_name)
    
    df.jointP10<- dplyr::bind_rows(df.jointP10)%>%
      dplyr::group_split(treatment,sample_name)
    
    #this would implement fold-change filters upon request acccording to Franken et al.
    if(isTRUE(filters)){
      #top3
      
      df.jointP3<-FC_filter(df.jointP3)
      #top 5
      df.jointP5<-FC_filter(df.jointP5)
      #top 10
      df.jointP10<-FC_filter(df.jointP10)
      #also filter original data by FC filters
      df1<-dplyr::bind_rows(df1)
      df1 <- suppressWarnings(df1 %>% dplyr::group_split(Accession,sample_id,sample_name) %>% 
                                purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n())))
      
      df1<-purrr::map(df1,function(x) FC_to_ref(x,baseline))
      
      if(isTRUE(filters)){
        df1<-purrr::map(df2,function(x) FC_filter(x))
        
      }
      df1<-dplyr::bind_rows(df1)
    }
    #check that all the data wasn't filtered out
    if(nrow(df.jointP3)==0){
      return(warning("Please disable filters, all data was filtered out for 6k peptides."))
    }
    if(nrow(df.jointP5)==0){
      return(warning("Please disable filters, all data was filtered out for top 12k peptides."))
    }
    if(nrow(df.jointP10)==0){
      return(warning("Please disable filters, all data was filtered out for all peptides."))
    }
    
    ## split data by Replicate 
    l.bytype3 <- split.data.frame(df.jointP3, df.jointP3$sample_id)
    l.bytype5 <- split.data.frame(df.jointP5, df.jointP5$sample_id)
    l.bytype10 <- split.data.frame(df.jointP10, df.jointP10$sample_id)
    
    ## determine which Replicate (F1 through FN in PD) contains the greatest number of PSM curves and use this for normalization
    n.filter3 <- lapply(l.bytype3, nrow)
    n.filter5 <- lapply(l.bytype5, nrow)
    n.filter10 <- lapply(l.bytype10, nrow)
    #subset the replicates by selecting the one with the greatest # of PSMs
    df.normP3 <- l.bytype3[[which.max(n.filter3)]]
    df.normP5 <- l.bytype5[[which.max(n.filter5)]]
    df.normP10<- l.bytype10[[which.max(n.filter10)]]
    #choose the accessions
    norm.accessions3 <- df.normP3$Accession
    norm.accessions5 <- df.normP5$Accession
    norm.accessions10 <- df.normP10$Accession
    ## turn lists back to data frames
    
    df.jointP3<-dplyr::bind_rows(df.jointP3)
    df.jointP5<-dplyr::bind_rows(df.jointP5)
    df.jointP10<-dplyr::bind_rows(df.jointP10)
    
    ## calculate median for each selected replicate
    df.median3 <- df.jointP3 %>%#top 6k psms
      dplyr::group_by(temperature) %>%
      dplyr::mutate(value = median(I,na.rm=TRUE))
    
    df.median5 <- df.jointP5 %>%#top 12k PSMs
      dplyr::group_by(temperature) %>%
      dplyr::mutate(value = median(I,na.rm=TRUE))
    
    df.median10 <- df.jointP10 %>%#all PSMs
      dplyr::group_by(temperature) %>%
      dplyr::mutate(value = median(I,na.rm=TRUE))
    
    
    #fit curves to median 
    ## fit curves to the median data for each replicate (F1 through FN)
    df.fit3 <- df.median3 %>%
      dplyr::group_by(sample_id) %>% 
      dplyr::mutate(fit = list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                       start = c(Pl=0, a = 550, b = 10),
                                       data = list(x=temperature,y=value),
                                       na.action = na.exclude,
                                       algorithm = "port",
                                       lower = c(0.0,1e-5,1e-5),
                                       upper = c(1.5,15000,300),
                                       control = nls.control(maxiter = 50)
      )),silent=TRUE))
    df.fit3<-df.fit3 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
    df.fit3<-df.fit3 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
    
    df.fit3<-dplyr::bind_rows(df.fit3)
    df.fit3<-df.fit3%>% #calculate fitted values and record as scaling factors for top 6k PSMs
      dplyr::mutate(fitted_values3 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
      dplyr::select(sample_id,fitted_values3,temperature,treatment) %>% dplyr::ungroup(.)
    
    df.fit5 <- df.median5 %>%
      dplyr::group_by(sample_id) %>% 
      dplyr::mutate(fit = list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                       start = c(Pl=0, a = 550, b = 10),
                                       data = list(x=temperature,y=value),
                                       na.action = na.exclude,
                                       algorithm = "port",
                                       lower = c(0.0,1e-5,1e-5),
                                       upper = c(1.5,15000,300),
                                       control = nls.control(maxiter = 50)
      )),silent=TRUE))
    df.fit5<-df.fit5 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
    df.fit5<-df.fit5 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
    df.fit5<-dplyr::bind_rows(df.fit5)
    df.fit5<- df.fit5%>% #calculate fitted values and record as scaling factors for top 12k PSMs
      dplyr::mutate(fitted_values5 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
      dplyr::select(sample_id,fitted_values5,temperature,treatment) %>% dplyr::ungroup(.)
    
    df.fit10 <- df.median10 %>%
      dplyr::group_by(sample_id) %>% 
      dplyr::mutate(fit = list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                       start = c(Pl=0, a = 550, b = 10),
                                       data = list(x=temperature,y=value),
                                       na.action = na.exclude,
                                       algorithm = "port",
                                       lower = c(0.0,1e-5,1e-5),
                                       upper = c(1.5,15000,300),
                                       control = nls.control(maxiter = 50)
      )),silent=TRUE))
    df.fit10<-df.fit10 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
    df.fit10<-df.fit10 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
    df.fit10<-dplyr::bind_rows(df.fit10)
    
    df.fit10<- df.fit10%>% #calculate fitted values and record as scaling factors for all PSMs
      dplyr::mutate(fitted_values10 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
      dplyr::select(sample_id,fitted_values10,temperature,treatment) %>% dplyr::ungroup(.)
    
    ## split data by replicate, treatment and temperature
    d3<-df.fit3 %>% dplyr::group_split(sample_id,treatment,temperature)
    d5<-df.fit5 %>% dplyr::group_split(sample_id,treatment,temperature)
    d10<-df.fit10 %>% dplyr::group_split(sample_id,treatment,temperature)
    
    
    #we only need the first row out of each group since the fitted values are nested
    check3 <-purrr::map(d3,function(x) x [1,])
    check5 <-purrr::map(d5,function(x) x [1,])
    check10 <-purrr::map(d10,function(x) x [1,])
    #unnest fitted values
    check3<-purrr::map(check3,function(x) x %>% tidyr::unnest(cols=fitted_values3) %>% unique(.))
    check5<-purrr::map(check5,function(x) x %>% tidyr::unnest(cols=fitted_values5) %>% unique(.)) 
    check10<-purrr::map(check10,function(x) x %>% tidyr::unnest(cols=fitted_values10) %>% unique(.)) #%>% dplyr::mutate(temperature=temperatures))
    # 
    #bind_rows
    check3<-dplyr::bind_rows(check3)
    check5<-dplyr::bind_rows(check5)
    check10<-dplyr::bind_rows(check10)
    #check column names in common between original data and PSM fitted_values
    name3<-dplyr::intersect(names(check3),names(df1))
    #perform a right join to add fitted_value column (this will be a scaling factor)
    test3<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check3,name3)
    test5<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check5,name3)
    test10<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check10,name3)
    
    ## calculate ratios between the fitted curves and the median values (this will be the correction factor)
    df.out3 <- test3 %>%
      dplyr::mutate(correction3 = ifelse(is.na(fitted_values /I),NA,fitted_values/I)) %>%
      dplyr::select('sample_id','temperature','I','fitted_values','correction3')
    df.out3<-df.out3 %>% dplyr::select(-fitted_values,-I)
    
    df.out5 <- test5 %>%
      dplyr::mutate(correction5 = ifelse(is.na(fitted_values / I),NA,fitted_values / I)) %>%
      dplyr::select('sample_id','temperature','I','fitted_values','correction5')
    df.out5<-df.out5 %>% dplyr::select(-fitted_values,-I)
    
    df.out10 <- test10 %>%
      dplyr::mutate(correction10 = ifelse(is.na(fitted_values /I),NA,fitted_values /I)) %>%
      dplyr::select('sample_id','temperature','I','fitted_values','correction10')
    df.out10<-df.out10 %>% dplyr::select(-fitted_values,-I)
    ## join correction factor to data
    df1$temperature<-as.factor(df1$temperature)
    df1<-df1 %>% dplyr::group_split(temperature)
    #group split data by temperature
    df.out3$temperature<-as.factor(df.out3$temperature)
    df.out3<-df.out3 %>% dplyr::group_by(temperature) %>% dplyr::group_split(.)
    
    df.out5$temperature<-as.factor(df.out5$temperature)
    df.out5<-df.out5 %>% dplyr::group_by(temperature) %>% dplyr::group_split(.)
    
    df.out10$temperature<-as.factor(df.out10$temperature)
    df.out10<-df.out10 %>% dplyr::group_by(temperature) %>% dplyr::group_split(.)
    #since the correction factor is global and is different for each temperature, split the data by temperature
    df1<-dplyr::bind_rows(df1)
    df1<-df1 %>% dplyr::group_split(temperature)
    #apply correction factor by temperature to original data
    df3<-purrr::map2(df1,df.out3,function(x,y)tryCatch(x %>% dplyr::mutate(correction3=y$correction3[1])))
    df5<-purrr::map2(df3,df.out5,function(x,y)tryCatch(x %>% dplyr::mutate(correction5=y$correction5[1])))
    df<-purrr::map2(df5,df.out10,function(x,y)tryCatch(x %>% dplyr::mutate(correction10=y$correction10[1])))
    #bind_rows
    df<-dplyr::bind_rows(df)
    
    df3 <- df%>% 
      dplyr::mutate(norm_value3= ifelse(is.na(correction3),I,I* correction3)) %>% dplyr::ungroup(.)
    df5 <- df3%>% 
      dplyr::mutate(norm_value5= ifelse(is.na(correction5),I,I* correction5)) %>% dplyr::ungroup(.)
    df<- df5%>% 
      dplyr::mutate(norm_value10= ifelse(is.na(correction10),I,I* correction10)) %>% dplyr::ungroup(.)
    
    df<-df %>% dplyr::mutate(norm_value3=ifelse(is.na(norm_value3),NA,norm_value3),
                             norm_value5=ifelse(is.na(norm_value5),NA,norm_value5),
                             norm_value10=ifelse(is.na(norm_value10),NA,norm_value10)
    )
    
    df <- df %>% 
      dplyr::rename("uniqueID"="Accession", "C"="temperature","I3"="norm_value3","I5"="norm_value5","I10"="norm_value10")
    df<-df %>% dplyr::ungroup(.) %>% dplyr::select(-I,-correction3,-correction5,-correction10)
    
    if(isTRUE(filters)){
      if (any(names(df)==T7)){
        df<-df %>% dplyr::select(-T10,-T7,-T9)
      }
    }
    return(df)
  }else{#if this is a protein file
    
    if(any(names(df)=="uniqueID")){
      df<-df %>% dplyr::rename("Accession"="uniqueID")
    }
    df$Accession<-as.factor(df$Accession)
    if(any(names(df)=="replicate")){
      df<-df %>%
        dplyr::group_split(Accession,sample_id,replicate)
    }else if(any(names(df)=="Fraction")){
      df<-df %>%
        dplyr::group_split(Accession,sample_id,Fraction)
    }
    
    
    #check baseline temperature
    baseline<-purrr::map(df,function(x) check_baseline(x))
    
    
    if(any(!isTRUE(order(unique(df[[1]]$temperature))==order(unique(temperatures$temperature))))){
      #if the temperatures are reversed, relabel 7,9th and 10th temperature channels
      df.jointP <- suppressWarnings(purrr::map2(df,baseline,function(x,y) FC_calc(x,y)))
      
    }else{
      df.jointP <- suppressWarnings(purrr::map2(df,baseline,function(x,y) FC_calc(x,y)))
    }
    
    df.jointP<- dplyr::bind_rows(df.jointP)
    
    if(isTRUE(filters)){
      df.jointP<-FC_filter(df.jointP)
    }
    if(nrow(df[[1]])==0){
      return(warning("Please disable filters, all data was filtered out from original treatment."))
    }
    
    ## split the data by replicate 
    l.bytype <- split.data.frame(df.jointP, df.jointP$sample_id)
    
    ## determine which replicate (F1 through FN) contains the greatest number of curves and use this for normalization
    n.filter <- lapply(l.bytype, nrow)
    df.normP <- l.bytype[[which.max(n.filter)]]
    norm.accessions <- df.normP$Accession
    
    ## calculate median for each sample_id group
    df<-dplyr::bind_rows(df)
    df.mynormset <- df %>% base::subset(Accession %in% norm.accessions)
    
    df.median <- df %>%
      dplyr::group_by(temperature) %>%
      dplyr::mutate(value = median(I,na.rm=TRUE)) %>% dplyr::ungroup(.)
    #
    # nls3 = purrr::quietly(.f = nls)
    # qtwolevel_fun = function(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl
    #                          start = c(Pl=0, a = 550, b = 10)
    #                          data = list(x=temperature,y=value)
    #                          na.action = na.exclude
    #                          algorithm = "port"
    #                          lower = c(0.0,1e-5,1e-5)
    #                          upper = c(1.5,15000,300)
    #                          control = nls.control(maxiter = 50))
    
    
    ## fit curves to the median data for each sample_id (F1 through FN)
    df.fit <- df.median %>%
      dplyr::group_by(sample_id) %>% 
      dplyr::mutate(fit=list(try(cetsa_fit(.,norm=FALSE))))
    # dplyr::mutate(fit = list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
    #                                  start = c(Pl=0, a = 550, b = 10),
    #                                  data = list(x=temperature,y=value),
    #                                  na.action = na.exclude,
    #                                  algorithm = "port",
    #                                  lower = c(0.0,1e-5,1e-5),
    #                                  upper = c(1.5,15000,300),
    #                                  control = nls.control(maxiter = 50))
    #)))
    df.fit<-df.fit %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
    df.fit<-df.fit %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
    df.fit<-dplyr::bind_rows(df.fit)
    df.fit<-df.fit%>% dplyr::group_by(sample_id) %>% 
      dplyr::mutate(fitted_values = ifelse(!is.logical(fit[[1]]),list(data.frame(fitted_values=stats::predict(fit[[1]]))),NA))  %>% 
      dplyr::select(sample_id,fitted_values,temperature) %>% dplyr::ungroup(.)
    
    d<-df.fit %>% dplyr::group_split(sample_id)
    check <-data.frame(temperature=unique(d[[1]]$temperature),
                       fitted_values=unique(d[[1]]$fitted_values[[1]]))
    #check<-purrr::map(check,function(x) x %>% unnest(fitted_values) %>% unique(.))#check<-purrr::map(d,function(x) x %>% unnest(c(fitted_values)) %>% unique(.) %>% dplyr::mutate(temperature=temperatures))
    check<-dplyr::bind_rows(check) %>% unique(.) 
    
    test<-df %>% dplyr::right_join(check,c('temperature'))
    
    
    ## calculate ratios between the fitted curves and the median values
    df.out <- test %>% dplyr::group_by(sample_id,temperature) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(correction = ifelse(is.na(fitted_values/I),NA,fitted_values/I)) %>%
      dplyr::select('temperature','correction','sample_id') %>%
      unique(.) %>%
      dplyr::ungroup(.) %>% 
      dplyr::group_by(sample_id) %>% 
      dplyr::group_split()
    
    df.jointP<-dplyr::bind_rows(df.jointP) %>% unique(.)
    df.jointP$temperature<-as.factor(df.jointP$temperature)
    df.jointP<-df.jointP %>% dplyr::group_split(sample_id) 
    
    #apply correction factor by temperature to original data
    df3<-purrr::map2(df.jointP,df.out,function(x,y)tryCatch(x %>% dplyr::mutate(correction=y$correction[1])))
    
    
    df3 <- dplyr::bind_rows(df3)%>% 
      dplyr::mutate(norm_value= ifelse(is.na(.$correction),.$I,.$I*correction)) %>% 
      dplyr::ungroup(.) 
    df3<-df3 %>% dplyr::mutate(norm_value=ifelse(is.na(.$norm_value),"NA",I)) %>% 
      dplyr::select(-I)
    df3 <- df3 %>% 
      dplyr::rename("uniqueID"="Accession", "C"="temperature","I"="norm_value")
    df3<-df3 %>% dplyr::ungroup(.) %>% dplyr::select(-correction)
    
    if(isTRUE(filters)){
      if (any(names(df3)==T7)){
        df3<-df3 %>% dplyr::select(-T10,-T7,-T9)
      }
    }
    return(df3)
  }
}