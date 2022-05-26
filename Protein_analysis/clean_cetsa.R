clean_cetsa <- function(df, temperatures = NULL,samples = NA,Peptide=FALSE,solvent,CFS=TRUE,CARRIER=TRUE,baseline="min"){
  df<-df %>% as.data.frame()
  if(isTRUE(CARRIER)){#if the carrier channel is present Default 131C
    df<-df %>% dplyr::filter(!temp_ref=="131C")
    temperatures<-temperatures %>% dplyr::filter(!temp_ref=="131C")
  }
  if(!any(names(df)=="sample")&any(names(df)=="sample_id")){
    df %>% dplyr::mutate(sample=sample_id)
  }
  if(!any(names(df)=="rank")){#if there's no rank column, add NA values to the df
    df$rank<-NA
  }
  if(any(names(df)=="I")&!any(names(df)=="value")){#if there's no rank column, add NA values to the df
    df$value<-df$I
  }
  samples<-samples
  if(any(stringr::str_detect(names(temperatures),"sample_name"))&any(stringr::str_detect(names(df),"replicate"))){#if this is a timepoint study
    if(!nchar(temperatures$sample_name[1])==nchar(df$sample_name[1])){
      hi<-df %>% dplyr::select(sample_name) %>% dplyr::mutate(sample=stringr::str_extract(sample_name,"_S[[:digit:]]+_"))
      hi<-hi %>% dplyr::mutate(sample=stringr::str_remove_all(sample,"_")) %>% distinct(.)
      df<-df %>% dplyr::right_join(hi,by="sample_name")
      temperatures<-temperatures %>% dplyr::select(-temp_ref) %>% 
        dplyr::mutate(dataset=ifelse(stringr::str_detect(.$sample,"no"),"vehicle","treated")) %>% dplyr::select(-sample)
      temperatures<-temperatures %>% dplyr::mutate(sample=.$sample_name,
                                                   replicate=as.character(.$Replicate)) %>% dplyr::select(-sample_name,-Replicate)
    }else if(stringr::str_detect(unique(df$sample_name),temperatures$sample_name[1])){
      hi<-df %>%
        dplyr::select(sample_name) %>%
        dplyr::mutate(sample=stringr::str_extract(sample_name,"_S[[:digit:]]+_"))
      hi<-hi %>% dplyr::mutate(sample=stringr::str_remove_all(sample,"_")) %>% distinct(.)
      
      df<-df %>% dplyr::right_join(hi,by="sample")
      temperatures<-temperatures %>% dplyr::select(-temp_ref) %>% 
        dplyr::mutate(dataset=ifelse(stringr::str_detect(.$sample,"no"),"vehicle","treated")) %>% dplyr::select(-sample)
      temperatures<-temperatures %>% dplyr::mutate(sample=.$sample_name,
                                                   replicate=as.character(.$Replicate)) %>% dplyr::select(-sample_name,-Replicate)
    }
  }
  
  if (is.null(temperatures)) {return(warning('No temperature data'))}
  #Rename column headers if necessary
  if(any(names(df)=="I10")){#if any Peptide groups were subset previously
    df<-df %>% mutate(value=I10)
  }
  if(any(names(df)=="uniqueID")&!any(names(df)=="Accession")){
    df<-df %>% dplyr::rename("Accession"="uniqueID")
  }
  
  df <- df %>% dplyr::ungroup(.)
  if(ncol(temperatures)>2){#If dealing with a time_point study 
    name<-dplyr::intersect(names(df),names(temperatures))
    df_ <- data.frame(merge(df,temperatures, all.y=TRUE))
    if(!any(names(df)=="sample_name")&!any(names(df)=="dataset")){
      df$sample_name<-df$Spectrum.File
      df_<-df
      df_<-df_ %>% dplyr::mutate(dataset=ifelse(stringr::str_detect(df_$sample_name,solvent),"vehicle","treated"))
    }else if(!is.na(samples)){
      df_ <- df_ %>%
        dplyr::right_join(samples, by = intersect(names(df_),names(samples)))
    }
    
  }else if (!any(stringr::str_detect(names(df),"temperature"))& is.na(samples)){#if samples data is missing
    hi<-dplyr::intersect(names(df),names(temperatures))
    df <- df %>%
      dplyr::right_join(temperatures, by = hi) 
  }
  
  df<-dplyr::bind_rows(df)
  if(baseline=="min"){
    bl <-min(df$temperature,na.rm=TRUE)
  }else{
    bl <-max(df$temperature,na.rm=TRUE)
  }
  
  if(any(names(df)=="Fraction")){
    if(any(stringr::str_detect(df$Fraction,"."))){
      df<-df %>% dplyr::mutate(Fraction = stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."))
    }
    df <- df %>%
      dplyr::filter(!is.na(.$Accession),
                    !is.na(.$temperature),
                    !is.na(.$value)) %>%
      dplyr::group_by(Accession,Annotated_Sequence,dataset,sample_name,Fraction) %>%
      dplyr::group_split()
  }else{
    df <- df %>%
      dplyr::filter(!is.na(.$Accession),
                    !is.na(.$temperature),
                    !is.na(.$value)) %>%
      dplyr::group_by(Accession,Annotated_Sequence,dataset,sample_name) %>%
      dplyr::distinct() %>% 
      dplyr::group_split()
  }
  
  
  if (.Platform$OS.type=="windows"){
    df<-parallel::mclapply(df,tryCatch({function(x){x<-x %>% 
      dplyr::mutate(value = x$value/mean(x$value[x$temperature == bl],na.rm=TRUE)) %>%
      distinct(.) %>% ungroup(.)
    return(x)
    }},
    error = function(cond){return(NA)}))
  }else{
    df<-parallel::mclapply(df,tryCatch({function(x){x<-x %>% 
      dplyr::mutate(value = x$value/mean(x$value[x$temperature == bl],na.rm=TRUE)) %>%
      distinct(.) %>% ungroup(.)
    return(x)
    }},
    error = function(cond){return(NA)}),mc.cores=availableCores())
  }
  df<-dplyr::bind_rows(df)
  if(!any(names(df)=="missing")){
    df<-df %>% dplyr::mutate(missing=is.na(value))
  }
  
  if(any(names(df)=="Fraction")){
    if(any(stringr::str_detect(df$Fraction,"."))){
      df<-df %>% dplyr::mutate(Fraction = stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."))
    }
  }
  
  df1<-df %>%
    dplyr::filter(temperature==min(temperature,na.rm=TRUE)) %>% 
    dplyr::mutate(rank=dplyr::ntile(value,3)) %>% 
    dplyr::select(Accession,dataset,rank,sample_name) %>%
    distinct(.)
  df1<-df1 %>% dplyr::group_by(Accession,sample_name) %>% 
    dplyr::group_split()
  df1<-furrr::future_map(df1,function(x) x %>% 
                           dplyr::filter(rank==min(rank,na.rm=TRUE)) %>% 
                           distinct())
  df1<-dplyr::bind_rows(df1)
  name<-dplyr::intersect(names(df),names(df1))
  df<-df %>% dplyr::left_join(df1,by=name)
  df<-dplyr::bind_rows(df)
  df$sample_id<-as.factor(df$sample_id)
  
  df<-dplyr::bind_rows(df) 
  return(df)
}
FC_to_ref<-function(x,baseline){ 
  if(baseline=="min"){
    y<-x %>%  
      dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==unique(x$temperature)[length(unique(x$temperature))-2,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)])),
                    T9 = try(mean(x[which(x$temperature==unique(x$temperature)[length(unique(x$temperature))-1,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)])),
                    T10 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)))
    return(y)
  }else{
    y<-x %>%  
      dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==unique(x$temperature)[min(unique(x$temperature),na.rm=TRUE)+2,"I"]$I/x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I)])),
                    T9 = try(mean(x[which(x$temperature==unique(x$temperature)[min(unique(x$temperature),na.rm=TRUE)+1,"I"]$I/x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I)])),
                    T10 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I/x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I)))
    return(y)
  }
}
FC_calc<-function(x,y) {
  y<-x %>%  
    dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==unique(x$temperature)[order(unique(x$temperature))=="7"]),]$I/x[which(x$temperature==y),"I"]$I)),
                  T9 = try(mean(x[which(x$temperature==unique(x$temperature)[order(unique(x$temperature))=="9"]),]$I/x[which(x$temperature==y),"I"]$I)),
                  T10 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I/x[which(x$temperature==y),"I"]$I)))
  return(y)
}
FC_filter<-function(x){
  y<-x%>% dplyr::filter(T7 >= 0.4, T7 <= 0.6,T9 < 0.3) %>% 
    dplyr::select(-T7,-T9,-n)
  if(any(names(x)=="T10" & all(!is.na(x$T10)))){
    y<- y %>% subset(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
  }
  return(y)
}
check_baseline<-function(x) {
  if(min(temperatures$temperature,na.rm=TRUE)==min(x$temperature,na.rm=TRUE)){
    baseline<-min(x$temperature,na.rm=TRUE)
  }else{
    warning("minimum temperature reporter ion values not found for this replicate")
    baseline<-min(x$temperature,na.rm=TRUE)
  }
  return(baseline)
}

