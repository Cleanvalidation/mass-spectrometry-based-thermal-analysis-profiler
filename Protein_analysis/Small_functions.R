mutate_missing<-function(x){
  if(any(names(x)=="replicate")){
    x<-x %>%  dplyr::group_by(uniqueID,Annotated_Sequence,dataset,sample_id,replicate) %>%
      dplyr::mutate(missing_pct=ifelse(all(is.na(.$value)),100,(100*(sum(as.numeric(is.na(.$value))))/nrow(.))))%>% ungroup(.)
    
  }else{
    x<-x %>%  dplyr::group_by(uniqueID,Annotated_Sequence,dataset,sample_id,Fraction) %>%
      dplyr::mutate(missing_pct=ifelse(all(is.na(.$value)),100,(100*(sum(as.numeric(is.na(.$value))))/nrow(.))))%>% ungroup(.)
  }
}


missing_label<-function(x) {
  x<-x %>% dplyr::group_by(Accession,sample_id,sample_name,dataset) %>%
    distinct(.) %>% dplyr::mutate(missing=is.na(value),
                                  missing_pct=sum(is.na(value))/nrow(.))
  return(x)
}




medianPolish <- function(intensities, num_channels){
  wide <- matrix(intensities, byrow = TRUE, ncol = num_channels)
  tmp_fit <- stats::medpolish(wide, na.rm = TRUE, trace.iter = FALSE)
  tmp_fit$overall + tmp_fit$col
}

rank_label<-function(x,baseline,Frac){
  if(isTRUE(Frac))
  {#if this is a fractionated dataset      
    if(baseline=="min"){
      x<-x%>% dplyr::filter(temperature==min(temperature,na.rm=TRUE)) %>%
        dplyr::mutate(rank=dplyr::ntile(value,3))%>%
        dplyr::select(sample_id,sample_name,Accession,rank,id,Fraction,Spectrum.File)%>%
        dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>% dplyr::select(-id) %>% distinct(.)
      
      
    }else{#if baseline is max and this is a fractionated dataset   
      x<-x%>% dplyr::filter(temp_ref==unique(x$temperature)[length(unique(x$temperature))]) %>%
        dplyr::mutate(rank=dplyr::ntile(value,3))%>%
        dplyr::select(sample_id,sample_name,Accession,rank,id,Fraction,Spectrum.File)%>%
        dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>% dplyr::select(-id) %>% distinct(.)
      
    }
  }
else{#if this is unfractionated
  if(baseline=="min"){#if baseline is min and this is an unfractionated dataset
    
    x<-x%>% dplyr::filter(temperature==min(temperature,na.rm=TRUE)) %>%
      dplyr::mutate(rank=dplyr::ntile(value,3))%>%
      dplyr::select(sample_id,sample_name,Accession,rank,id,Spectrum.File)%>%
      dplyr::filter(!is.na(rank),!is.na(id)) %>% dplyr::select(-id) %>% distinct(.)
    
    
  }else{#if baseline is max and this is an unfractionated dataset
    
    x<-x%>% dplyr::filter(temp_ref==unique(x$temperature)[length(unique(x$temperature))]) %>%
      dplyr::mutate(rank=dplyr::ntile(value,3))%>%
      dplyr::select(sample_id,sample_name,Accession,rank,id,Spectrum.File)%>%
      dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>%
      dplyr::select(-id) %>% distinct(.)
    
  }
}
return(x)
}




