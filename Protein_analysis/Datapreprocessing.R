#Renaming the protein, peptide, psms


read_xl<-function(x,solvent,CFS) {
  df2<-readxl::read_xlsx(x,trim_ws=TRUE,.name_repair = function(x) gsub("[[:punct:][:blank:]]+",".",x))
  df2<-df2 %>% dplyr::select(-tidyr::starts_with("Found"))
  
  if(any(names(df2)=="Percolator.PEP.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"="Percolator.PEP.by.Search.Engine.Sequest.HT")
  }else if(any(names(df2)=="Percolator.PEP")){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"="Percolator.PEP")
  }
  if(any(names(df2)=="Charge.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Charge"="Charge.by.Search.Engine.Sequest.HT")
    
  }
  
  if(any(names(df2)=="DeltaM.ppm.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.by.Search.Engine.Sequest.HT")
    
  }else if(any(names(df2)=="DeltaM.ppm.")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.")
  }
  if(any(names(df2)=="XCorr.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("XCorr"="XCorr.by.Search.Engine.Sequest.HT")
  }
  if(any(names(df2)==".Missed.Cleavages")){
    df2<-df2 %>% dplyr::rename("MissedCleavages"=".Missed.Cleavages")
  }
  if(any(names(df2)=="Annotated.Sequence")){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"="Annotated.Sequence") %>% 
      dplyr::mutate(Annotated_Sequence=toupper(Annotated_Sequence))
  }
  if (any(names(df2)=="Master.Protein.Accessions")){
    df2<-df2 %>% dplyr::rename("Accession"="Master.Protein.Accessions")
  }
  if(any(names(df2)=="Average.Reporter.S.N")){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"="Average.Reporter.S.N")
  }
  if(any(names(df2)=="Isolation.Interference.")){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"="Isolation.Interference.")
  }
  if(any(names(df2)=="Ion.Inject.Time.ms.")){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"="Ion.Inject.Time.ms.")
  }
  if(any(stringr::str_detect(names(df2),".PSMs"))){
    df2<-df2 %>% dplyr::rename("Num_PSMs"=".PSMs")
  }
  if(any(stringr::str_detect(colnames(df2),"Spectrum.File"))){
    if(any(stringr::str_detect(df2$Spectrum.File,"NOcarrier"))|isTRUE(CFS)){
      df2<-df2 %>% dplyr::mutate(sample_name=paste0(ifelse(stringr::str_detect(Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                                                    ifelse(stringr::str_detect(Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                                                    ifelse(stringr::str_detect(Spectrum.File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(Spectrum.File,"S_Phi")==TRUE,"S",NA))),
                                 treatment=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
                                 CC=ifelse(stringr::str_detect(Spectrum.File,solvent),0,1))
      if(any(stringr::str_detect(df2$File.ID,"."))){
        df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                   sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
      }else{
        df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                   sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"),
                                   sample_name = ifelse(stringr::str_detect(Spectrum.File,solvent),solvent,stringr::str_extract(stringr::str_to_lower(Spectrum.File),"[[:lower:]]+_[[:digit:]]+_")))
      }
      
    }else{
      df2<-df2 %>% dplyr::mutate(
        treatment=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
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

read_PD<-function(x,solvent){#rename data
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
    df2<-df2 %>% dplyr::rename("Percolator_PEP"="Percolator.PEP.by.Search.Engine.Sequest.HT")
  }else if(any(names(df2)=="Percolator.PEP")){
    df2<-df2 %>% dplyr::rename("Percolator_PEP"="Percolator.PEP")
  }
  if(any(names(df2)=="Charge.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("Charge"="Charge.by.Search.Engine.Sequest.HT")
    
  }
  
  if(any(names(df2)=="DeltaM.ppm.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.by.Search.Engine.Sequest.HT")
    
  }else if(any(names(df2)=="DeltaM.ppm.")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.")
  }
  if(any(names(df2)=="XCorr.by.Search.Engine.Sequest.HT")){
    df2<-df2 %>% dplyr::rename("XCorr"="XCorr.by.Search.Engine.Sequest.HT")
  }
  if(any(names(df2)==".Missed.Cleavages")){
    df2<-df2 %>% dplyr::rename("MissedCleavages"=".Missed.Cleavages")
  }
  if(any(names(df2)=="Annotated.Sequence")){
    df2<-df2 %>% dplyr::rename("Annotated_Sequence"="Annotated.Sequence")%>% 
      dplyr::mutate(Annotated_Sequence=toupper(Annotated_Sequence))
  }
  if (any(names(df2)=="Master.Protein.Accessions")){
    df2<-df2 %>% dplyr::rename("Accession"="Master.Protein.Accessions")
  }
  if(any(names(df2)=="Average.Reporter.S.N")){
    df2<-df2 %>% dplyr::rename("Average_Reporter_S/N"="Average.Reporter.S.N")
  }
  if(any(names(df2)=="Isolation.Interference.")){
    df2<-df2 %>% dplyr::rename("Isolation_Interference_[%]"="Isolation.Interference.")
  }
  if(any(names(df2)=="Ion.Inject.Time.ms.")){
    df2<-df2 %>% dplyr::rename("Ion_Inject_Time_[ms]"="Ion.Inject.Time.ms.")
  }
  if(any(names(df2)=="DeltaM.ppm.by.Search.Engine.Sequest")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.by.Search.Engine.Sequest.HT")
    
  }else if(any(names(df2)=="DeltaM.ppm.")){
    df2<-df2 %>% dplyr::rename("DeltaM"="DeltaM.ppm.")
  }
  df2<-df2 %>% 
    tidylog::pivot_longer(cols=colnames(df2)[stringr::str_detect(colnames(df2),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                          names_to = "id",
                          values_to ="value") %>% 
    dplyr::mutate(treatment=ifelse(stringr::str_detect(.$id,solvent),"vehicle","treated"),
                  CC=ifelse(stringr::str_detect(.$id,solvent),0,1),
                  sample_id = stringr::str_extract(.$id,"F[[:digit:]]+"),
                  temp_ref = unlist(stringr::str_extract(.$id,"[:digit:][:digit:][:digit:][N|C]|126|131")),
                  value = as.numeric(value))
  
  return(df2)
}
read_cetsa <- function(protein_path,peptide_path,Prot_Pattern,Peptide=FALSE,Frac=TRUE,CFS=TRUE,solvent="DMSO",CARRIER=TRUE,rank=TRUE,sub=NA,temperatures=temps,baseline="min",NORM="QUANTILE"){
  file.list<-protein_path
  i=1
  peptide_path<-as.character(peptide_path)
  protein_path<-as.character(protein_path)
  
  find<-c('[:digit:][:digit:][:digit:][N|C]|126|131')
  df<-list()
  df1<-list()
  df3<-data.frame()
  df2<-list()
  df.raw<-list(data.frame())
  
  setwd(protein_path)
  f<-list.files(protein_path,pattern="Proteins.xlsx")
  setwd(peptide_path)
  h<-list.files(peptide_path,pattern="PSMs.xlsx")
  
  
  #read_PSMs and proteins
  if (.Platform$OS.type=="windows"){
    Proteins<-parallel::mclapply(f,read_xl,solvent=solvent)
    PSMs<-parallel::mclapply(h,read_xl,solvent=solvent,CFS=CFS)
  }else{
    PSMs<-parallel::mclapply(h,read_xl,solvent=solvent,CFS=CFS,mc.cores = availableCores())
  }
  
  if(any(stringr::str_detect(Peptide,c("PG","PSMs")))){
    
    if(!is.na(sub)){#if there's only a subset of PSMs
      n<-as.numeric(sub)
      #subset a set number of PSMs
    
    }
    
    if (Peptide=="PSMs"){ #if this is a PSM file
      #split PSM subset into lists by way of sample_id name
      PSMs<-PSMs[1:n]
      PSMs<-dplyr::bind_rows(PSMs) %>% dplyr::group_split(sample_name)
      
      #only select abundance columns
      mat_norm<-which(stringr::str_detect(names(PSMs[[1]]),"Abundance"))
      if(isTRUE(CARRIER)){
        mat_norm<-which(stringr::str_detect(names(PSMs[[1]]),"Abundance")&!stringr::str_detect(names(PSMs[[1]]),"Abundance.131C"))
      }
      #calculate scaling factor
      df2<-furrr::future_map(PSMs,function(x)
      { 
        if(any(names(x)=="Fraction")|any(stringr::str_detect(x$File.ID[1],"."))){
          x<-dplyr::bind_rows(x) %>%
            dplyr::mutate(Fraction=stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."),
                          File.ID =stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
          
          #}
          x<-x %>%
            dplyr::ungroup() %>%
            dplyr::group_by(Accession,Annotated_Sequence,treatment,Fraction,sample_id) %>% 
            dplyr::group_split()
        }else{
          x<-x %>%
            dplyr::ungroup() %>%
            dplyr::group_by(Accession,Annotated_Sequence,treatment,sample_id)%>% 
            dplyr::group_split()
        }
        x<-lapply(x,function(y)tryCatch(
          {
            if(baseline=="min"){
              
              x<-y%>%
                dplyr::mutate(scal_fac = list(as.vector(colMeans(y[,mat_norm],na.rm=TRUE)/colMeans(y[,mat_norm[1]],na.rm=TRUE))))
            }else{
              
              x<-y %>%
                dplyr::mutate(scal_fac = list(as.vector(colMeans(y[,mat_norm],na.rm=TRUE)/colMeans(y[,mat_norm[length(mat_norm)]],na.rm=TRUE))))
            }
            
            return(x)
          },
          error=function(cond){
            return(NA)
          },
          finally={
          }
        )
        )
        x<-dplyr::bind_rows(x)
        return(x)
      }
      )
      #if any data has the bioreplicate number, truncate the sample_id name contents
      if(any(stringr::str_detect(unique(dplyr::bind_rows(df2)$sample_name),"_[[:digit:]]+"))){
        df2<-dplyr::bind_rows(df2) %>%
          dplyr::mutate(
            sample_name=stringr::str_remove(.$sample_name,"_[[:digit:]]+_"))
      }
      df2<-dplyr::bind_rows(df2) %>% 
        dplyr::group_by(treatment,sample_name) %>% 
        dplyr::group_split()
      
      df2<-purrr::map(df2,function(x) choose_PSM(x,
                                                 Frac=Frac,
                                                 NORM=NORM,
                                                 CARRIER=CARRIER,
                                                 sub=sub,
                                                 baseline=baseline,
                                                 Prot_Pattern=Prot_Pattern,
                                                 Peptide=TRUE
      )
      )
      
      df2<-df2 %>% purrr::keep(function(x) any(class(x)=="data.frame"))
      #make sure the data is ready to be processed by sample_id name
      if(any(names(df2)=="sample_name")){
        df2<-dplyr::bind_rows(df2) %>% 
          dplyr::group_by(sample_name) %>%
          dplyr::group_split(.)
      }
      if(any(stringr::str_detect(names(df2[[1]]),"Accession"))){
        df2<-dplyr::bind_rows(df2) %>% dplyr::rename("uniqueID"="Accession")
      }
      if(isTRUE(Frac)){
        PSMs<-dplyr::bind_rows(df2) %>% dplyr::select(uniqueID,Spectrum.File,sample_name,sample_id,Fraction)
      }else{
        PSMs<-dplyr::bind_rows(df2) %>% dplyr::select(uniqueID,Spectrum.File,sample_name,sample_id)
      }
      
      if(any(stringr::str_detect(names(df2),"Protein.Accessions"))&!any(stringr::str_detect(names(df2),"uniqueID"))){
        df2<-df2 %>% dplyr::rename("uniqueID"="Protein.Accessions")
      }
      #df2 is PSM data
      df2<-df2 %>%
        dplyr::group_by(sample_name) %>%
        dplyr::group_split()
      
      
      df2<-furrr::future_map(df2,function(x) mutate_missing(x))
      df2<-dplyr::bind_rows(df2)
      df2<-dplyr::bind_rows(df2) %>% 
        tidylog::pivot_longer(cols=colnames(df2)[stringr::str_detect(colnames(df2),"Abundances")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(treatment= ifelse(stringr::str_detect(.$Spectrum.File,solvent),"vehicle","treated"),
                      CC= ifelse(stringr::str_detect(.$Spectrum.File,solvent),0,1),
                      sample_id = stringr::str_extract(.$id,"[[:upper:]]+[[:digit:]]+"),
                      temp_ref = stringr::str_extract(.$id,"[:digit:][:digit:][:digit:][N|C]|126|131"),
                      value = as.numeric(value),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))) %>% 
        dplyr::group_by(sample_name) %>% 
        dplyr::group_split()
      if(!any(stringr::str_detect(names(dplyr::bind_rows(df2)),"temperature"))){
        df2<-df2 %>% dplyr::right_join(temperatures)
      }
      return(df2)
    }else if(Peptide=="PG"){
      if(length(list.files(peptide_path,pattern="PeptideGroups.xlsx"))>0){
        f<-list.files(peptide_path,pattern="PeptideGroups.xlsx")
      }
      #first get Peptide group file read
      if (.Platform$OS.type=="windows"){
        df.raw<-parallel::mclapply(f,read_xl,solvent=solvent,CFS=CFS)
      }else{
        df.raw<-parallel::mclapply(f,read_xl, solvent=solvent,CFS=CFS,mc.cores = availableCores())
      }
      #get row number for peptide groups in case of batch files present
      df.raw1<-purrr::map2(df.raw,seq(df.raw),function(x,y) x %>% dplyr::mutate(n=y))
      #names<-purrr::map2(df.raw,f,function(x,y) ifelse(length(names(x))<45,warning(paste0("Please check the columns on file names",y)),print("All files have all necessary columns")))
      df.raw<-dplyr::bind_rows(df.raw1)
      
      
      
      PG<-read_PD(df.raw,solvent=solvent)#read in peptide groups
      
      PSMs<- furrr::future_map(PSMs,function(x) x %>%
                                 dplyr::select(Accession,File.ID,Spectrum.File) %>%
                                 distinct())
      #no need to look at fractions since peptide groups are already summarized
      PSMs<-dplyr::bind_rows(PSMs) %>%
        dplyr::mutate(File.ID = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
      
      PG<-dplyr::bind_rows(PG) %>%
        dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"))
      PG<-data.table::data.table(dplyr::bind_rows(PG))
      PSMs<-data.table::data.table(dplyr::bind_rows(PSMs))
      
      name<-dplyr::intersect(names(PSMs),names(PG))#dfP is PSMs df3 is proteins
      
      #Join Peptide Group and PSM file
      PG<-PG %>%
        distinct(.) %>%
        dplyr::right_join(PSMs,by=name)
      PG<-PG %>% dplyr::filter(!is.na(sample_id))
      PG <- PG %>%
        dplyr::mutate(sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+")))
      
      
      df3<-dplyr::bind_rows(PG) %>%
        as.data.frame()
      if(isTRUE(CFS)){
        df3$sample_name<-paste0(ifelse(stringr::str_detect(df3$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(df3$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                                ifelse(stringr::str_detect(df3$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(df3$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                                ifelse(stringr::str_detect(df3$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(df3$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
      }
      
      
      df3<-dplyr::bind_rows(df3) %>%
        dplyr::group_split(sample_name)#peptide groups
      
      return(df3)
      
    }
  }else{ #if this is a protein file
    
    PSMs<- furrr::future_map(PSMs,function(x) x %>%
                               dplyr::select(Accession,File.ID,Spectrum.File) %>%
                               dplyr::mutate(File.ID=stringr::str_extract(File.ID,"[:upper:][[:digit:]]+")) %>% 
                               distinct())
    Proteins<-dplyr::bind_rows(Proteins)
    if(isTRUE(Frac)){
      Proteins <- Proteins %>% 
        tidylog::pivot_longer(cols=colnames(Proteins)[stringr::str_detect(colnames(Proteins),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"),
                      Fraction=stringr::str_remove(id,"[:upper:][[:digit:]]+."),
                      temp_ref = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|126|131'),
                      value = as.numeric(value),
                      Fraction = stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."),
                      sample_id = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
      Proteins<-Proteins %>% dplyr::right_join(temperatures, by ="temp_ref")
      
    }else{
      Proteins <- Proteins %>% 
        tidylog::pivot_longer(cols=colnames(Proteins)[stringr::str_detect(colnames(Proteins),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(File.ID = stringr::str_extract(id,"[:upper:][[:digit:]]+"),
                      temp_ref = stringr::str_extract(id,'[:digit:][:digit:][:digit:][N|C]|126|131'),
                      value = as.numeric(value),
                      Fraction = stringr::str_remove(File.ID,"[:upper:][[:digit:]]+."),
                      sample_id = stringr::str_extract(File.ID,"[:upper:][[:digit:]]+"))
      Proteins<-Proteins %>% dplyr::right_join(temperatures, by ="temp_ref")
    }
    Proteins<-data.table::data.table(dplyr::bind_rows(Proteins))
    PSMs<-data.table::data.table(dplyr::bind_rows(PSMs))
    
    name<-dplyr::intersect(names(Proteins),names(PSMs))#dfP is PSMs df3 is proteins
    #Join protein and PSM file
    Proteins<-Proteins  %>%
      dplyr::right_join(PSMs,by=name) %>%  #Add spectrum_File values
      distinct(.)
    Proteins<-Proteins %>% tibble::as_tibble()
    #Select protein columns
    if(any(names(Proteins)=="Biological.Process")&any(names(Proteins[[1]]=="MW.kDa."))){
      Proteins<-dplyr::bind_rows(Proteins) %>% 
        dplyr::rename("Cell_Component"="Cellular.Component","MW_kDa"="MW.kDa.","Bio_Process"="Biological.Process","Coverage"="Coverage.","Molecular_Function"="Molecular.Function")
    }else if(any(names(Proteins)=="MW.kDa.")){
      Proteins<-dplyr::bind_rows(Proteins) %>% 
        dplyr::rename("MW_kDa"="MW.kDa.","Coverage"="Coverage.")
    }else if (any(stringr::str_detect(names(Proteins),"Coverage."))){
      Proteins<-dplyr::bind_rows(Proteins) %>% 
        dplyr::rename("Coverage"="Coverage.")
    }
    #Pivot longer
    if(isTRUE(Frac)){
      Proteins <- Proteins %>% 
        dplyr::mutate(CC = ifelse(stringr::str_detect(Spectrum.File,solvent),0,1),
                      treatment = ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+")))
      Proteins$treatment<-as.factor(Proteins$treatment)
    }else{#if this isnt fractionated
      Proteins <- Proteins %>%  
        dplyr::mutate(CC = ifelse(stringr::str_detect(Spectrum.File,solvent),0,1),
                      treatment = ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+")))
      
      Proteins$treatment<-as.factor(Proteins$treatment)
    }
    
    
    Proteins<-dplyr::bind_rows(Proteins) %>% dplyr::group_split(Accession,sample_id,sample_name,treatment) 
    
    if (.Platform$OS.type=="windows"){
      Proteins<-parallel::mclapply(Proteins,missing_label)
    }else{
      Proteins<-parallel::mclapply(Proteins,missing_label,mc.cores=availableCores())
    }
    
    if(isTRUE(rank)){
      if (.Platform$OS.type=="windows"){
        df_raw_D_R<-parallel::mclapply(Proteins,rank_label,baseline=baseline,Frac=Frac)
      }else{
        df_raw_D_R<-parallel::mclapply(Proteins,rank_label,baseline=baseline,Frac=Frac,mc.cores=availableCores())
      }
      if(any(names(Proteins[[1]]$Spectrum.File)=="Fraction")){
        Proteins<-dplyr::bind_rows(Proteins)
        Proteins$Fraction<-stringr::str_extract(Proteins$Spectrum.File,"Fraction_[[:digit:]]+")
        Proteins$Fraction<-stringr::str_remove(Proteins$Fraction,"Fraction_")
        
        df_raw_D_R<-dplyr::bind_rows(df_raw_D_R)
        df_raw_D_R$Fraction<-stringr::str_extract(df_raw_D_R$Spectrum.File,"Fraction_[[:digit:]]+")
        df_raw_D_R$Fraction<-stringr::str_remove(df_raw_D_R$Fraction,"Fraction_")
      }
      #ID will contain accessions with high-intensity values at reference channel
      ID<-dplyr::intersect(unique(Proteins$Accession),unique(df_raw_D_R$Accession)) %>%
        na.omit(.)
      
      Proteins<-dplyr::bind_rows(Proteins)%>% 
        dplyr::filter(Accession %in% ID) %>% distinct(.)
      df_raw_D_R<-dplyr::bind_rows(df_raw_D_R)%>% 
        dplyr::filter(Accession %in% ID) %>% distinct(.)
      
      name<-dplyr::intersect(names(df_raw_D_R),names(Proteins))
      #df_raw_D_R is ranked protein data frame
      df_raw_D_R<-data.table::data.table(df_raw_D_R)
      Proteins<-data.table::data.table(Proteins)
      #set key to join ranked protein data frame to protein file
      data.table::setkeyv(df_raw_D_R,cols=name)
      data.table::setkeyv(Proteins,cols=name)
      df2<-data.table::merge.data.table(Proteins,df_raw_D_R,by=name)
      df2<-df2 %>% as.data.frame(.) %>% distinct(.)
      
    }else{
      if(any(names(Proteins)=="Fraction")){
        Proteins<-dplyr::bind_rows(Proteins)
        Proteins$Fraction<-stringr::str_extract(Proteins$Spectrum.File,"Fraction_[[:digit:]]+")
        Proteins$Fraction<-stringr::str_remove(Proteins$Fraction,"Fraction_")
        
        df_raw_D_R<-dplyr::bind_rows(df_raw_D_R)
        df_raw_D_R$Fraction<-stringr::str_extract(df_raw_D_R$Spectrum.File,"Fraction_[[:digit:]]+")
        df_raw_D_R$Fraction<-stringr::str_remove(df_raw_D_R$Fraction,"Fraction_")
      }
      df2<-dplyr::bind_rows(Proteins)
      if(isTRUE(CFS)){
        df2$sample_name<-paste0(ifelse(stringr::str_detect(df2$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(df2$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                                ifelse(stringr::str_detect(df2$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(df2$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                                ifelse(stringr::str_detect(df2$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(df2$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
      }
      
    }
  }
  df2<-dplyr::bind_rows(df2) %>% distinct(.)
  return(df2)
}

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
      dplyr::group_by(Accession,Annotated_Sequence,sample_name,treatment,Fraction) %>% 
      dplyr::group_split()
  }else{
    df2<-dplyr::bind_rows(df2) %>%
      dplyr::filter(Accession %in%united1$Accession & Annotated_Sequence %in% united1$Annotated_Sequence) %>% 
      dplyr::group_by(Accession,Annotated_Sequence,sample_name,treatment) %>% 
      dplyr::group_split()
  }
  mat_norm<-which(stringr::str_detect(names(df2[[1]]),"Abundance"))
  
  #transform abundance to log2 abundance
  df2_log<-furrr::future_map(df2,function(x) log2(x[,mat_norm]))
  #apply selected normalization method
  if(is.na(NORM)){
    check<-(df2_log)
  }
  else if(NORM=="QUANTILE"){
    print("norm",Norm)
    check<-rowmeanquantile(df2_log)
  }
  else if(NORM=="EQ_Median"){
    check<-furrr::future_map(tryCatch(DEqMS::equalMedianNormalization(as.matrix(x)),
                                      error = function(e) {print(e)
                                        return(NA)},
                                      finally = {print("Equal median normalization finished")})
    )
  }
  #convert log2check are now abundance values 
  df2<-purrr::map2(df2,check,tryCatch({function(x,y)
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
  
  #group data by sample_id name
  df2<-dplyr::bind_rows(df2) %>%
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
    },mc.cores = future::availableCores())
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
    dplyr::group_by(Accession,Annotated_Sequence,sample_name,sample_id,treatment) %>%
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
    dplyr::group_by(Accession,sample_name,sample_id,treatment)
  #data is grouped by  protein, bioreplicate and condition to aggregate PSMs to proteins
  united<-df2 %>% 
    dplyr::select(-id) %>% 
    dplyr::filter(!is.na(temp_ref)) %>%
    dplyr::distinct(.) %>% 
    dplyr::group_by(Accession,Annotated_Sequence,treatment,sample_name,sample_id) %>% 
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
    dplyr::group_by(Accession,sample_name,treatment) %>%
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
  #Join all the PSM data together with aggregated abundances and split by method or sample_id name
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
  #   dplyr::group_split(sample_name)
  # ggplot(df2[[1]],mapping=aes(x=temp_ref,y=value))+geom_boxplot()+ylim(0,5)
  # #
  
  if(any(stringr::str_detect(names(df2[[1]]),"Accession"))){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("uniqueID"="Accession"))
  }
  return(df2)
}
clean_cetsa <- function(df, temperatures = NULL,samples = NA,Peptide=FALSE,solvent,CFS=TRUE,CARRIER=TRUE,baseline="min"){
  df<-df %>% as.data.frame()
  if(isTRUE(CARRIER)){#if the carrier channel is present Default 131C
    df<-df %>% dplyr::filter(!temp_ref=="131C")
    temperatures<-temperatures %>% dplyr::filter(!temp_ref=="131C")
  }
  if(!any(names(df)=="sample_id")&any(names(df)=="sample_id")){
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
      hi<-df %>% dplyr::select(sample_name) %>% dplyr::mutate(sample_id=stringr::str_extract(sample_name,"_S[[:digit:]]+_"))
      hi<-hi %>% dplyr::mutate(sample_id=stringr::str_remove_all(sample_id,"_")) %>% distinct(.)
      df<-df %>% dplyr::right_join(hi,by="sample_name")
      temperatures<-temperatures %>% dplyr::select(-temp_ref) %>% 
        dplyr::mutate(treatment=ifelse(stringr::str_detect(.$sample_id,"no"),"vehicle","treated")) %>% dplyr::select(-sample_id)
      temperatures<-temperatures %>% dplyr::mutate(sample_id=.$sample_name,
                                                   replicate=as.character(.$Replicate)) %>% dplyr::select(-sample_name,-Replicate)
    }else if(stringr::str_detect(unique(df$sample_name),temperatures$sample_name[1])){
      hi<-df %>%
        dplyr::select(sample_name) %>%
        dplyr::mutate(sample_id=stringr::str_extract(sample_name,"_S[[:digit:]]+_"))
      hi<-hi %>% dplyr::mutate(sample_id=stringr::str_remove_all(sample_id,"_")) %>% distinct(.)
      
      df<-df %>% dplyr::right_join(hi,by="sample_id")
      temperatures<-temperatures %>% dplyr::select(-temp_ref) %>% 
        dplyr::mutate(treatment=ifelse(stringr::str_detect(.$sample_id,"no"),"vehicle","treated")) %>% dplyr::select(-sample_id)
      temperatures<-temperatures %>% dplyr::mutate(sample_id=.$sample_name,
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
    if(!any(names(df)=="sample_name")&!any(names(df)=="treatment")){
      df$sample_name<-df$Spectrum.File
      df_<-df
      df_<-df_ %>% dplyr::mutate(treatment=ifelse(stringr::str_detect(df_$sample_name,solvent),"vehicle","treated"))
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
      dplyr::group_by(Accession,Annotated_Sequence,treatment,sample_name,Fraction) %>%
      dplyr::group_split()
  }else{
    df <- df %>%
      dplyr::filter(!is.na(.$Accession),
                    !is.na(.$temperature),
                    !is.na(.$value)) %>%
      dplyr::group_by(Accession,Annotated_Sequence,treatment,sample_name) %>%
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
    dplyr::select(Accession,treatment,rank,sample_name) %>%
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
mutate_missing<-function(x){
  if(any(names(x)=="replicate")){
    x<-x %>%  dplyr::group_by(uniqueID,Annotated_Sequence,treatment,sample_id,replicate) %>%
      dplyr::mutate(missing_pct=ifelse(all(is.na(.$value)),100,(100*(sum(as.numeric(is.na(.$value))))/nrow(.))))%>% ungroup(.)
    
  }else{
    x<-x %>%  dplyr::group_by(uniqueID,Annotated_Sequence,treatment,sample_id,Fraction) %>%
      dplyr::mutate(missing_pct=ifelse(all(is.na(.$value)),100,(100*(sum(as.numeric(is.na(.$value))))/nrow(.))))%>% ungroup(.)
  }
}

missing_label<-function(x) {
  x<-x %>% dplyr::group_by(Accession,sample_id,sample_name,treatment) %>%
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




