set.seed(123)

theme_set(theme_bw())

#' Read in data
#'
#' Read in Excel file and apply minimal pre-processing
#'
#'
#' @param  f.  Path of Excel file output from Proteome Discoverer
#' @return a dataframe containing extracted information
#'
#' @importFrom readxl read_excel
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom stringr stringr::str_extract
#' @export
#' 
read_cetsa <- function(protein_path,peptide_path,Prot_Pattern,Peptide=FALSE,Frac=TRUE,CFS=TRUE,solvent="DMSO",CARRIER=TRUE,rank=TRUE,sub=NA,temperatures=temps,baseline="min",NORM="QUANTILE"){
  #print("printing solvent",solvent)
  file.list<-protein_path
  i=1
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
    PSMs<-parallel::mclapply(h,read_xl,solvent=solvent)
  }else{
    PSMs<-parallel::mclapply(h,read_xl,solvent=solvent,mc.cores = availableCores())
  }
  
  
  if(any(stringr::str_detect(Peptide,c("PG","PSMs")))){
    
    if(!is.na(sub)){#if there's only a subset of PSMs
      n<-as.numeric(sub)
      #subset a set number of PSMs
      PSMs<-PSMs[1:n]
    }
    
    if (Peptide=="PSMs"){ #if this is a PSM file
      #split PSM subset into lists by way of sample name
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
            dplyr::group_by(Accession,Annotated_Sequence,dataset,Fraction,sample_id) %>% 
            dplyr::group_split()
        }else{
          x<-x %>%
            dplyr::ungroup() %>%
            dplyr::group_by(Accession,Annotated_Sequence,dataset,sample_id)%>% 
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
            message("Scaling factor calculated per protein,annotated sequence and dataset")
          }
        )
        )
        x<-dplyr::bind_rows(x)
        return(x)
      }
      )
      #if any data has the bioreplicate number, truncate the sample name contents
      if(any(stringr::str_detect(unique(dplyr::bind_rows(df2)$sample_name),"_[[:digit:]]+"))){
        df2<-dplyr::bind_rows(df2) %>%
          dplyr::mutate(
            sample_name=stringr::str_remove(.$sample_name,"_[[:digit:]]+_"))
      }
      df2<-dplyr::bind_rows(df2) %>% 
        dplyr::group_by(dataset,sample_name) %>% 
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
      #make sure the data is ready to be processed by sample name
      df2<-dplyr::bind_rows(df2) %>% 
        dplyr::group_by(sample_name) %>%
        dplyr::group_split(.)
      
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
        dplyr::mutate(dataset= ifelse(stringr::str_detect(.$Spectrum.File,"DMSO"),"vehicle","treated"),
                      CC= ifelse(stringr::str_detect(.$Spectrum.File,"DMSO"),0,1),
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
        df.raw<-parallel::mclapply(f,read_xl,solvent=solvent)
      }else{
        df.raw<-parallel::mclapply(f,read_xl, solvent=solvent,mc.cores = availableCores())
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
        dplyr::mutate(CC = ifelse(stringr::str_detect(Spectrum.File,"DMSO"),0,1),
                      dataset = ifelse(stringr::str_detect(Spectrum.File,"DMSO"),"vehicle","treated"),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+")))
      Proteins$treatment<-as.factor(Proteins$dataset)
    }else{#if this isnt fractionated
      Proteins <- Proteins %>%  
        dplyr::mutate(CC = ifelse(stringr::str_detect(Spectrum.File,"DMSO"),0,1),
                      dataset = ifelse(stringr::str_detect(Spectrum.File,"DMSO"),"vehicle","treated"),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+")))
      
      Proteins$treatment<-as.factor(Proteins$dataset)
    }
    
    
    Proteins<-dplyr::bind_rows(Proteins) %>% dplyr::group_split(Accession,sample_id,sample_name,dataset) 
  
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