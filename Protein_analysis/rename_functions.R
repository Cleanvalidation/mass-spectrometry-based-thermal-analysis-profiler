#Renaming the protein, peptide, psms


read_xl<-function(x,solvent) {

  
  df2<-readxl::read_xlsx(x,trim_ws=TRUE,.name_repair = function(x) gsub("[[:punct:][:blank:]]+",".",x))
  df2<-df2 %>% dplyr::select(-tidyr::starts_with("Found"))
  
  if(any(names(df2)=="Percolator.PEP.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"=names(df2)[stringr::str_detect(names(df2),"PEP")])
  }else if(any(stringr::str_detect(names(df2),"PEP"))){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"=names(df2)[stringr::str_detect(names(df2),"PEP")])
  }
  if(any(names(df2)=="Charge.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Charge"=names(df2)[stringr::str_detect(names(df2),"Charge")])
  }else if (any(stringr::str_detect(names(df2),"Charge"))){
     df2<-df2 %>% dplyr::rename("Charge"=names(df2)[stringr::str_detect(names(df2),"Charge")])
    }
  
  if(any(names(df2)=="DeltaM.ppm.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.by.Search.Engine.Sequest.HT")
    
  }else if(any(stringr::str_detect(names(df2),"DeltaM"))){
    df2<-df2 %>% dplyr::rename("DeltaM"=names(df2)[stringr::str_detect(names(df2),"DeltaM")])
  }
  if(any(names(df2)=="XCorr.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("XCorr"="XCorr.by.Search.Engine.Sequest.HT")
  }else if(any(stringr::str_detect(names(df2),"XCorr"))){
    df2<-df2 %>% dplyr::rename("XCorr"=names(df2)[stringr::str_detect(names(df2),"XCorr")])
}
  if(any(names(df2)==".Missed.Cleavages")){
    df2<-df2 %>% dplyr::rename("MissedCleavages"=".Missed.Cleavages")
  }else if(any(stringr::str_detect(names(df2),"Missed"))){
    df2<-df2 %>% dplyr::rename("XCorr"=names(df2)[stringr::str_detect(names(df2),"XCorr")])
}
  if(any(names(df2)=="Annotated.Sequence")){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"="Annotated.Sequence") %>% 
      dplyr::mutate(Annotated_Sequence=toupper(Annotated_Sequence))
  }else if(any(stringr::str_detect(names(df2),"Annotated"))){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"=names(df2)[stringr::str_detect(names(df2),"Annotated")])
}
  if (any(names(df2)=="Master.Protein.Accessions")){
    df2<-df2 %>% dplyr::rename("Accession"="Master.Protein.Accessions")
  }else if(any(stringr::str_detect(names(df2),"Master"))){
    df2<-df2 %>% dplyr::rename("Accession"=names(df2)[stringr::str_detect(names(df2),"Master")])
}else if(any(stringr::str_detect(names(df2),"Accession"))){
    df2<-df2 %>% dplyr::rename("Accession"=names(df2)[stringr::str_detect(names(df2),"Accession")])
}
  if(any(names(df2)=="Average.Reporter.S.N")){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"="Average.Reporter.S.N")
  }else if(any(stringr::str_detect(names(df2),"Reporter"))){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"=names(df2)[stringr::str_detect(names(df2),"Reporter")])
}
  if(any(names(df2)=="Isolation.Interference.")){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"="Isolation.Interference.")
  }else if(any(stringr::str_detect(names(df2),"Isolation"))){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"=names(df2)[stringr::str_detect(names(df2),"Isolation")])
}
  if(any(names(df2)=="Ion.Inject.Time.ms.")){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"="Ion.Inject.Time.ms.")
  }else if(any(stringr::str_detect(names(df2),"Inject"))){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"=names(df2)[stringr::str_detect(names(df2),"Inject")])
}
  if(any(stringr::str_detect(names(df2),".PSMs"))){
    df2<-df2 %>% dplyr::rename("Num_PSMs"=".PSMs")
  }else if (any(stringr::str_detect(names(df2),"PSMs"))){
    df2<-df2 %>% dplyr::rename("Num_PSMs"=names(df2)[stringr::str_detect(names(df2),"PSM")])
    }
  if(any(stringr::str_detect(colnames(df2),"Spectrum.File"))){#if this is a PSM file
    if(any(stringr::str_detect(df2$Spectrum.File,"NOcarrier"))){#if this experiment has Carrier FAIMS and PhiSDM involved, shorten the file names
      df2<-df2 %>% dplyr::mutate(sample_name=paste0(ifelse(stringr::str_detect(Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                                                    ifelse(stringr::str_detect(Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                                                    ifelse(stringr::str_detect(Spectrum.File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(Spectrum.File,"S_Phi")==TRUE,"S",NA))),
                                 dataset=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
                                 CC=ifelse(stringr::str_detect(Spectrum.File,solvent),0,1))
      if(any(stringr::str_detect(df2$File.ID,"."))){
        df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                   sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
      }else{#if theres no fractionation
        df2<-df2 %>% dplyr::mutate(sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
      }
      
    }else{#if this experiment does not have Carrier FAIMS and PhiSDM involved shorten he file names
      df2<-df2 %>% dplyr::mutate(
        dataset=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
        sample_name = ifelse(stringr::str_detect(Spectrum.File,solvent),solvent,stringr::str_extract(stringr::str_to_lower(Spectrum.File),"[[:lower:]]+_[[:digit:]]+_")),
        CC=ifelse(stringr::str_detect(Spectrum.File,solvent),0,1))
      if(any(stringr::str_detect(df2$File.ID,"."))){
        df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                   sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
      }else{
        df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                   sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
      }
      
    }
    
  }
  return(df2)
}


read_PD<-function(x,solvent){#rename data for peptide groups

  df2<-x %>% as.data.frame(.) %>% 
    dplyr::select(names(x)[stringr::str_detect(names(x),'Master')],
                  tidyselect::starts_with('File'),
                  tidyselect::starts_with('Abundance.')|tidyselect::starts_with('1')|tidyselect::starts_with('Abundances.'),
                  tidyselect::starts_with('Annotated'),
                  tidyselect::contains('Isolation'),
                  tidyselect::contains('Ion'),
                  tidyselect::starts_with('Charge'),
                  tidyselect::contains('PEP'),
                  tidyselect::starts_with('Percolator'),
                  tidyselect::starts_with('.Missed'),
                  tidyselect::contains('Modifications'),
                  tidyselect::contains('Cleavages'),
                  tidyselect::starts_with('XCorr'),
                  tidyselect::contains('Delta'),
                  tidyselect::contains('File'),
                  tidyselect::contains('S/N'),
                  tidyselect::starts_with('Average'),
                  tidyselect::contains('Spectrum'),
                  tidyselect::contains('.'),
                  tidyselect::ends_with('PSMs'),
                  -tidyselect::contains('Grouped'),
                  -tidyselect::starts_with('Found')) 
  if(any(names(df2)=="Percolator.PEP.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"=names(df2)[stringr::str_detect(names(df2),"PEP")])
  }else if(any(stringr::str_detect(names(df2),"PEP"))){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"=names(df2)[stringr::str_detect(names(df2),"PEP")])
  }
  if(any(names(df2)=="Charge.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Charge"=names(df2)[stringr::str_detect(names(df2),"Charge")])
  }else if (any(stringr::str_detect(names(df2),"Charge"))){
     df2<-df2 %>% dplyr::rename("Charge"=names(df2)[stringr::str_detect(names(df2),"Charge")])
    }
  
  if(any(names(df2)=="DeltaM.ppm.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.by.Search.Engine.Sequest.HT")
    
  }else if(any(stringr::str_detect(names(df2),"DeltaM"))){
    df2<-df2 %>% dplyr::rename("DeltaM"=names(df2)[stringr::str_detect(names(df2),"DeltaM")])
  }
  if(any(names(df2)=="XCorr.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("XCorr"="XCorr.by.Search.Engine.Sequest.HT")
  }else if(any(stringr::str_detect(names(df2),"XCorr"))){
    df2<-df2 %>% dplyr::rename("XCorr"=names(df2)[stringr::str_detect(names(df2),"XCorr")])
}
  if(any(names(df2)==".Missed.Cleavages")){
    df2<-df2 %>% dplyr::rename("MissedCleavages"=".Missed.Cleavages")
  }else if(any(stringr::str_detect(names(df2),"Missed"))){
    df2<-df2 %>% dplyr::rename("XCorr"=names(df2)[stringr::str_detect(names(df2),"XCorr")])
}
  if(any(names(df2)=="Annotated.Sequence")){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"="Annotated.Sequence") %>% 
      dplyr::mutate(Annotated_Sequence=toupper(Annotated_Sequence))
  }else if(any(stringr::str_detect(names(df2),"Annotated"))){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"=names(df2)[stringr::str_detect(names(df2),"Annotated")])
}
  if (any(names(df2)=="Master.Protein.Accessions")){
    df2<-df2 %>% dplyr::rename("Accession"="Master.Protein.Accessions")
  }else if(any(stringr::str_detect(names(df2),"Master"))){
    df2<-df2 %>% dplyr::rename("Accession"=names(df2)[stringr::str_detect(names(df2),"Master")])
}else if(any(stringr::str_detect(names(df2),"Accession"))){
    df2<-df2 %>% dplyr::rename("Accession"=names(df2)[stringr::str_detect(names(df2),"Accession")])
}
  if(any(names(df2)=="Average.Reporter.S.N")){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"="Average.Reporter.S.N")
  }else if(any(stringr::str_detect(names(df2),"Reporter"))){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"=names(df2)[stringr::str_detect(names(df2),"Reporter")])
}
  if(any(names(df2)=="Isolation.Interference.")){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"="Isolation.Interference.")
  }else if(any(stringr::str_detect(names(df2),"Isolation"))){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"=names(df2)[stringr::str_detect(names(df2),"Isolation")])
}
  if(any(names(df2)=="Ion.Inject.Time.ms.")){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"="Ion.Inject.Time.ms.")
  }else if(any(stringr::str_detect(names(df2),"Inject"))){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"=names(df2)[stringr::str_detect(names(df2),"Inject")])
}
  if(any(stringr::str_detect(names(df2),".PSMs"))){
    df2<-df2 %>% dplyr::rename("Num_PSMs"=".PSMs")
  }else if (any(stringr::str_detect(names(df2),"PSMs"))){
    df2<-df2 %>% dplyr::rename("Num_PSMs"=names(df2)[stringr::str_detect(names(df2),"PSM")])
    }
  df2<-df2 %>% 
    tidylog::pivot_longer(cols=colnames(df2)[stringr::str_detect(colnames(df2),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                          names_to = "id",
                          values_to ="value") %>% 
    dplyr::mutate(dataset=ifelse(stringr::str_detect(.$id,solvent),"vehicle","treated"),
                  CC=ifelse(stringr::str_detect(.$id,solvent),0,1),
                  sample_id = stringr::str_extract(.$id,"F[[:digit:]]+"),
                  temp_ref = unlist(stringr::str_extract(.$id,"[:digit:][:digit:][:digit:][N|C]|126|131")),
                  value = as.numeric(value))
  
  return(df2)
}

