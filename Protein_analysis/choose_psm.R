choose_PSM<-function(x,Frac=Frac,NORM=NORM,CARRIER=CARRIER,subset=subset,baseline=baseline,Prot_Pattern=Prot_Pattern,Peptide=Peptide){
  if(isTRUE(CARRIER)){
    x<-x[,!stringr::str_detect(colnames(x),"131C")]
  }
  
  df2<- dplyr::bind_rows(x)
  DF2<-dplyr::bind_rows(x)
  #collect a sub-subset of PSMs present for each protein,psm_sequence, replicate/condition and show all replicates where PSM sequences are present
  united<-dplyr::bind_rows(df2) %>%
    dplyr::select(Accession,Annotated_Sequence,sample_id,sample_name) %>% 
    tidyr::pivot_wider(.,
                       names_from="sample_id",
                       values_from="sample_id",
                       values_fill=NA,
                       values_fn=unique) %>% 
    dplyr::filter(!is.na(Accession))#remove data with missing protein accession values
  
  # #determine the threshold for minimum number of replicates present for normalization > 40% default
  if(ncol(united)==4){
    numsamples<-0.50*(ncol(united)-2)
    united1<-united[rowSums(!is.na(united))-2>=numsamples,]#remove PSMs with 40% or less replicates
    print(paste0("Removed ",round(100*(nrow(united)-nrow(united1))/nrow(united),2),
                 "% of PSMs where each sequence has less than 50% of the total bioreplicates present"))
  }else{
    numsamples<-0.80*(ncol(united)-2)
    united1<-united[rowSums(!is.na(united))-2>=numsamples,]#remove PSMs with 40% or less replicates
    print(paste0("Removed ",round(100*(nrow(united)-nrow(united1))/nrow(united),2),
                 "% of PSMs where each sequence has less than 80% of the total bioreplicates present"))
  }
  
  #from the PSM subset, keep those accessions where the conditional min # of replicates are met
  if(any(names(dplyr::bind_rows(df2)=="Fraction"))){
    df2<-dplyr::bind_rows(df2) %>%
      dplyr::filter(Accession %in%united1$Accession & Annotated_Sequence %in% united1$Annotated_Sequence) %>% 
      dplyr::group_by(Accession,Annotated_Sequence,sample_name,dataset,Fraction) %>% 
      dplyr::group_split()
  }else{
    df2<-dplyr::bind_rows(df2) %>%
      dplyr::filter(Accession %in%united1$Accession & Annotated_Sequence %in% united1$Annotated_Sequence) %>% 
      dplyr::group_by(Accession,Annotated_Sequence,sample_name,dataset) %>% 
      dplyr::group_split()
  }
  mat_norm<-which(stringr::str_detect(names(df2[[1]]),"Abundance"))
  
  #transform abundance to log2 abundance
  df2_log<-furrr::future_map(df2,function(x) log2(x[,mat_norm]))
  #apply selected normalization method
  if(is.na(NORM)){
    check<-(df2_log)
  }else if(NORM=="QUANTILE"){
    check<-furrr::future_map(df2_log,function(x) {
      tryCatch(limma::normalizeBetweenArrays(as.matrix(x),ties=FALSE,method="quantile"),
               error = function(e) {print(e)
                 return(NA)},
               finally = {print("Quantile normalization finished.")})
    }
    )
  }else if(NORM=="EQ_Median"){
    check<-furrr::future_map(tryCatch(DEqMS::equalMedianNormalization(as.matrix(x)),
                                      error = function(e) {print(e)
                                        return(NA)},
                                      finally = {print("Equal median normalization finished")})
    )
  }
  #convert log2check are now abundance values 
  df2_<-purrr::map2(df2,check,tryCatch({function(x,y)
    if(is.na(y[[1]])){
      return(x)
    }else{
      x<-x[,!stringr::str_detect(colnames(x),"Abundance")]
      x<-cbind(x,y)
      return(x)
    }
  },
  error = function(e) {message("Error adding normalized data to the original data frame")
    return(NA)},
  finally = {print("Adding normalized data to original data frame finished.")}
  ))
  
  check<-NA
  
  #group data by sample name
  df2<-dplyr::bind_rows(df2_) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::group_split()
  
  #if the column names are in PD format, change them
  if(any(names(df2[[1]])=="Average.Reporter.S.N")){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("Average_Reporter_S/N"="Average.Reporter.S.N")%>% 
                      dplyr::mutate(sample_id=File.ID))
  }
  if(any(names(df2[[1]])=="Isolation.Interference.")){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("Isolation_Interference_[%]"="Isolation.Interference.") %>% 
                      dplyr::mutate(sample_id=File.ID))
  }
  if(any(names(df2[[1]])=="Ion.Inject.Time.ms.")){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("Ion_Inject_Time_[ms]"="Ion.Inject.Time.ms."))
    
  }
  if(any(names(df2[[1]])=="Percolator.PEP.")){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("Percolator_PEP"="Percolator.PEP") %>% 
                      dplyr::mutate(sample_id=File.ID))
    
  }
  
  #change PSMs from wide to long and create TMT channel "temp_ref" and abundance "value" columns
  if (.Platform$OS.type=="windows"){
    df2<-parallel::mclapply(df2,function(x){
      x1 <- x %>% 
        tidylog::pivot_longer(cols=colnames(x)[stringr::str_detect(colnames(x),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(temp_ref = unlist(stringr::str_remove(id,'Abundance.')),
                      value = as.numeric(value))
      return(x1)
    })
  }else{
    df2<-parallel::mclapply(df2,function(x){
      x1 <- x %>% 
        tidylog::pivot_longer(cols=colnames(x)[stringr::str_detect(colnames(x),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(temp_ref = unlist(stringr::str_remove(id,'Abundance.')),
                      value = as.numeric(value))
      return(x1)
    },mc.cores = availableCores())
  }
  # if(any(stringr::str_detect(names(df2[[1]]),"File.ID"))){
  #   df2<-purrr::map(df2,function(x) x %>% dplyr::rename("sample_id"="File.ID"))
  # }
  # #uncomment to Diagnose the plot after normalization uncomment for debugging
  # df2<-dplyr::bind_rows(df2) %>%
  #   dplyr::group_split(sample_name)
  # ggplot(df2,mapping=aes(x=temp_ref,y=value))+geom_boxplot()
  # #
  df2<-dplyr::bind_rows(df2) %>% 
    dplyr::group_by(Accession,Annotated_Sequence,sample_name,sample_id,dataset) %>%
    dplyr::distinct() %>% 
    dplyr::group_split()
  #if the sample_name isn't shortened, shorten it
  if(any(stringr::str_detect(dplyr::bind_rows(df2)$sample_name,"_[[:digit:]]+_")) & !isTRUE(CARRIER)){
    df2<-dplyr::bind_rows(df2) %>%#add replicate value from the name assigned in PD if it exists
      dplyr::mutate(replicate=stringr::str_remove_all(stringr::str_extract(sample_name,"_[[:digit:]]+_"),"_"))
    df2$sample_name<-sub("_[[:digit:]]+_", "", dplyr::bind_rows(df2)$sample_name)
    df2$sample_name<-sub("[[:digit:]]+", "", dplyr::bind_rows(df2)$sample_name)
    
  }
  df2<-dplyr::bind_rows(df2) %>% 
    dplyr::group_by(Accession,sample_name,sample_id,dataset)
  #data is grouped by  protein, bioreplicate and condition to aggregate PSMs to proteins
  united<-df2 %>% 
    dplyr::select(-id) %>% 
    dplyr::filter(!is.na(temp_ref)) %>%
    dplyr::distinct(.) %>% 
    dplyr::group_by(Accession,Annotated_Sequence,dataset,sample_name,sample_id) %>% 
    dplyr::group_split()
  
  
  #for each protein, bioreplicate and condition, nest PSM values and TMT channels
  united<-purrr::map(united,function(x) x %>% 
                       dplyr::mutate(temp_ref = as.character(temp_ref),
                                     value = as.double(value)) %>%
                       tidyr::nest(data = c(value,temp_ref)))
  
  united<-dplyr::bind_rows(united)
  #Create a new variable to pivot from long to wide before aggregating PSM values
  United<- united %>% 
    dplyr::mutate(data=furrr::future_map(data,function(y)y %>%  
                                           tidyr::pivot_wider(.,#need to transform to wide for PSM aggregation
                                                              names_from = "temp_ref",
                                                              values_from ="value",
                                                              names_repair =unique,
                                                              values_fn = unique,
                                                              values_fill = NA
                                           )
    )
    )
  United<-United %>%
    dplyr::group_by(Accession,sample_name,dataset) %>%
    dplyr::group_split()#Split the data to get multiple peptides of the same master protein
  
  #Merge abundances contribution from each sequence back to the data frame
  United<-purrr::map(United,function(x){
    x<-x %>% dplyr::mutate(data = list(dplyr::bind_rows(x$data)))
    return(x)
  }
  )
  #apply median polish to perform local normalization from the PSM to the protein level
  United<-purrr::map(United,function(x){
    x<-x %>%
      dplyr::mutate(data=list(
        medianPolish(as.matrix(x$data[[1]]),ncol(x$data[[1]]))))
    return(x)
  })
  #reverse log2 transform then transpose and change matrix to data frame
  # United<-purrr::map(United,function(x){
  #   x<-x %>% dplyr::mutate(data=list(2^x$data[[1]]))
  #   return(x)
  # })
  #each accession has redundant list of data frames after aggregation, therefore we can select the aggregated abundance
  United1<-purrr::map(United,tryCatch({function(x){
    x<- x %>% dplyr::mutate(data=list(data.frame(t(x$data[[1]]))))
    y<-magrittr::set_names(x$data[[1]],unique(united$data[[1]]$temp_ref))
    return(y)
  }
  },
  error = function(e){print(e)
    return(NA)
  },
  finally = {'Finalizing PSM data cleanup.'})
  )
  
  #append aggregated PSM values to the data
  United<-purrr::map2(United,United1,function(x,y)
    x %>%
      dplyr::select(-data) %>%
      dplyr::bind_cols(y))
  #Join all the PSM data together with aggregated abundances and split by method or sample name
  United<-dplyr::bind_rows(United) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::group_split()
  #pivot data to long format after aggregating PSM abundances to protein level
  df2<-furrr::future_map(United,function(x) x %>%
                           tidyr::pivot_longer(cols=colnames(x)[stringr::str_detect(colnames(x),c('[:digit:][:digit:][:digit:][N|C]|126|131'))],
                                               names_to="temp_ref",
                                               values_to="value"))
  #multiply value by the scaling factor
  df2<-furrr::future_map(df2,function(x)x %>% 
                           dplyr::mutate(value = ifelse(!is.na(x$value)&(!is.na(x$scal_fac[[1]])),x$value*x$scal_fac[[1]],NA)))
  
  #Keep proteins where accession is not NA and nrow of peptides has information available
  df2<-purrr::keep(df2,function(x) !is.logical(x$Accession)&nrow(x)>0)
  #group data and split by sample_name
  df2<-dplyr::bind_rows(df2) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::group_split(.)
  df2<-purrr::map(df2,function(x) clean_cetsa(x,#normalize data from 0 to 1
                                              temperatures = df.temps,
                                              Peptide=Peptide,
                                              solvent=solvent,
                                              CFS=CFS,
                                              CARRIER=CARRIER,
                                              baseline=baseline))
  
  df2<-purrr::keep(df2,function(x)
    !is.logical(x$Accession) & nrow(x)>0
  )
  
  # #uncomment to Diagnose the plot after normalization uncomment for debugging
  # df2<-dplyr::bind_rows(df2) %>%
  #   dplyr::group_split(sample_name)
  # ggplot(df2[[1]],mapping=aes(x=temp_ref,y=value))+geom_boxplot()+ylim(0,5)
  # #
  
  if(any(stringr::str_detect(names(df2[[1]]),"Accession"))){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("uniqueID"="Accession"))
  }
  return(df2)
}

