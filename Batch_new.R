library(minpack.lm)
library(rlist)
library(data.table)
library(knitr)
library(ggthemes)
library(gridExtra)
library(grid)
library(readxl)
library(nls2)
library(stats)
library(ggplot2)
library(pkgcond)
library(rlist)
library(pracma)
library(fs)
library(tidyverse)
library(splines)
library(mgcv)
library(purrr)
library(nlstools)
library(stringr)
library(stringi)
library(mice)
library(DBI)
library(furrr)
library(tibble)
library(ComplexUpset)
library(patchwork)
library(caret)
library(ggpubr)
library(furrr)
library(parallel)
library(janitor)
library(ggrepel)
library(ggpubr)
requireNamespace("plyr")


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
read_cetsa <- function(protein_path,peptide_path,Prot_Pattern,Peptide=FALSE,Batch=TRUE,CFS=TRUE,solvent="DMSO"){
  
  file.list<-protein_path
  i=1
  peptide_path<-as.character(peptide_path)
  protein_path<-as.character(protein_path)
  find<-c('[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]')
  df<-list()
  df1<-list()
  df3<-data.frame()
  df2<-list()
  df.raw<-list(data.frame())
  setwd(protein_path)
  f<-list.files(protein_path,pattern="Proteins.xlsx")
  setwd(peptide_path)
  h<-list.files(peptide_path,pattern="PSMs.xlsx")
  read_xl<-function(x) {
    readxl::read_xlsx(x,trim_ws=TRUE,.name_repair = function(x) gsub("[[:punct:][:blank:]]+",".",x))
  }
  #read_PSMS and proteins
  Proteins<-mclapply(f,read_xl,mc.cores = availableCores())
  PSMs<-mclapply(h,read_xl,mc.cores = availableCores())
  
  if (isTRUE(Peptide)){ #if this is a PSM of peptide group file
    if(length(list.files(peptide_path,pattern="PeptideGroups.xlsx"))>0){
      f<-list.files(peptide_path,pattern="PeptideGroups.xlsx")
    }
    #first get Peptide file read
    df.raw<-mclapply(f, read_xl, mc.cores = availableCores())
    #get row number for peptide groups in case of batch files present
    df.raw<-purrr::map2(df.raw,seq(df.raw),function(x,y) x %>% dplyr::mutate(n=y))
    PSMs<-purrr::map2(PSMs,seq(PSMs),function(x,y) x %>% dplyr::mutate(n=y))
    #names<-purrr::map2(df.raw,f,function(x,y) ifelse(length(names(x))<45,warning(paste0("Please check the columns on file names",y)),print("All files have all necessary columns")))
    
    df.raw<-dplyr::bind_rows(df.raw)
    PSMs<-dplyr::bind_rows(PSMs)
    
    read_PD<-function(x) {#rename data
      df2<-x %>% as.data.frame(.) %>% 
        dplyr::select(names(x)[stringr::str_detect(names(x),'Master')],
                      tidyselect::starts_with('File'),
                      tidyselect::starts_with('Abundance')|tidyselect::starts_with('1'),
                      tidyselect::starts_with('Annotated'),
                      tidyselect::contains('Isolation'),
                      tidyselect::contains('Ion'),
                      tidyselect::starts_with('Charge'),
                      tidyselect::contains('PEP'),
                      tidyselect::contains('Modifications'),
                      tidyselect::contains('Cleavages'),
                      tidyselect::starts_with('XCorr'),
                      tidyselect::contains('Delta'),
                      tidyselect::contains('File'),
                      tidyselect::contains('S/N'),
                      tidyselect::contains('Spectrum'),
                      -tidyselect::contains('Grouped'),
                      -tidyselect::starts_with('Found'),n) 
      if(!(any(names(df2)=="Percolator.PEP.by.Search.Engine.Sequest.HT"))){
        df2$Percolator.PEP.by.Search.Engine.Sequest.HT<-0
      }
      df2$Charge<-df2$Charge.by.Search.Engine.Sequest.HT
      df2$Percolator_PEP<-df2$Percolator.PEP.by.Search.Engine.Sequest.HT
      df2$DeltaM<-df2$DeltaM.ppm.by.Search.Engine.Sequest.HT
      df2$XCorr<-df2$XCorr.by.Search.Engine.Sequest.HT
      df2$MissedCleavages<-df2$.Missed.Cleavages
      df2$Annotated_Sequence<-df2$Annotated.Sequence
      df2<-df2 %>% dplyr::rename("Accession"="Master.Protein.Accessions")
      
    }
    df.raw<-read_PD(df.raw)
    PSMs<-read_PD(PSMs)
    
    
    
    if(any(names(df.raw)=="Spectrum.File")){
      df2<-dplyr::bind_rows(df.raw)
      df2<-df2 %>% dplyr::mutate(Spectrum_File=str_replace(df2$Spectrum.File,"[:punct:][[:digit:]]+.raw",paste0(Prot_Pattern,".xlsx")))
      df2<-df2 %>% dplyr::group_split(Spectrum_File)
      
      if(any(names(df.raw)=="File.ID")){#if this is fractionated
        df2<-furrr::future_map(df2,.f=function(x) x %>%  
                                 dplyr::select(Accession,
                                               tidyselect::starts_with('Abundance'),
                                               tidyselect::starts_with('File'),
                                               Spectrum_File,Annotated_Sequence,
                                               tidyselect::contains('Isolation'),
                                               tidyselect::contains('Ion'),
                                               tidyselect::contains('Charge'),
                                               tidyselect::contains('PEP'),
                                               tidyselect::contains('Modifications'),
                                               tidyselect::contains('Cleavages'),
                                               tidyselect::starts_with('XCorr'),
                                               tidyselect::contains('Delta'),
                                               tidyselect::contains('File'),
                                               tidyselect::contains('S/N')) %>% 
                                 tidyr::gather('id', 'value', -Accession,-File.ID,-Spectrum_File,-Annotated_Sequence,-Modifications,-Charge,-"Average_Reporter_S/N",-XCorr,-"Isolation_Interference_[%]",-"MissedCleavages",-"Percolator_PEP",-"Ion_Inject_Time_[ms]",-"DeltaM") %>% 
                                 dplyr::mutate(temp_ref = stringr::str_extract_all(id,find),
                                               sample_id= stringr::str_extract_all(id,"F[[:digit:]]+")) %>% 
                                 tidyr::unnest(cols=c("temp_ref","sample_id")) %>% 
                                 dplyr::select(-id,Accession,File.ID,Spectrum_File,tidyr::contains("Sequence"),temp_ref,value,Modifications,Charge,"Average_Reporter_S/N",XCorr,"Isolation_Interference_[%]","MissedCleavages","Percolator_PEP","Ion_Inject_Time_[ms]","DeltaM"))
      }else{#if this is unfractionated
        df2<-furrr::future_map(df2,.f=function(x) x %>%  
                                 dplyr::select(Accession,
                                               tidyselect::starts_with('Abundance'),
                                               Spectrum_File,Annotated_Sequence,
                                               tidyselect::contains('Isolation'),
                                               tidyselect::contains('Ion'),
                                               tidyselect::contains('Charge'),
                                               tidyselect::contains('PEP'),
                                               tidyselect::contains('Modifications'),
                                               tidyselect::contains('Cleavages'),
                                               tidyselect::starts_with('XCorr'),
                                               tidyselect::contains('Delta'),
                                               tidyselect::contains('File'),
                                               tidyselect::contains('S/N')) %>% 
                                 tidyr::gather('id', 'value', -Accession,-Spectrum_File,-Annotated_Sequence,-Modifications,-Charge,-"Average_Reporter_S/N",-XCorr,-"Isolation_Interference_[%]",-"MissedCleavages",-"Percolator_PEP",-"Ion_Inject_Time_[ms]",-"DeltaM") %>% 
                                 dplyr::mutate(temp_ref = stringr::str_extract_all(id,find),
                                               sample_id= stringr::str_extract_all(id,"F[[:digit:]]+")) %>% 
                                 tidyr::unnest(cols=c("temp_ref","sample_id")) %>% 
                                 dplyr::select(-id,Accession,Spectrum_File,tidyr::contains("Sequence"),temp_ref,value,Modifications,Charge,"Average_Reporter_S/N",XCorr,"Isolation_Interference_[%]","MissedCleavages","Percolator_PEP","Ion_Inject_Time_[ms]","DeltaM"))
        
        df2<-dplyr::bind_rows(df2)
      }
    }else{
      
      if(all(c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n")%in% names(df.raw) )==FALSE){
        check<-names(df.raw)[which(names(df.raw) %in% c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n"))]
        warning(paste0("Please check columns from Protein Group PD output:",setdiff(c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n"),check)))
      }
      df2<-data.table(df.raw) 
      df3<-melt(data = df2, 
                id.vars = c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n"),
                measure.vars = names(df2)[grepl( "Abundance" , names( df2 ) )],
                variable.name = "id",
                value.name = "value")
      df3<-df3[, id:=as.character(id)]
      df3<-df3 %>% dplyr::group_split(n)
      #extract temperature channel and sample_id from abundance columns
      cols_PD<-function(x){x %>% dplyr::mutate(temp_ref = unlist(str_extract(x$id,find)),
                                               sample_id= unlist(str_extract(x$id,"F[[:digit:]]+"))) %>% 
          dplyr::select(-id)
      }
      df3<-mclapply(df3,cols_PD,mc.cores=availableCores())
      
      if(isTRUE(Batch)){
        df2<-purrr::map2(df3,f,function(x,y)x %>% dplyr::mutate(sample_name=str_remove(y,"[[:digit:]]+[:punct:]PeptideGroups.xlsx")))
      }else{
        
        PSMs<-PSMs %>% dplyr::rename("sample_name"="Spectrum.File") 
        if(any(names(PSMs)=="File.ID")){
          PSMs<-PSMs %>% dplyr::mutate(sample_id=str_extract(File.ID,"F[[:digit:]]+"))
          PSMs<-PSMs%>% dplyr::select(sample_id,sample_name) %>% unique(.)
        }else{
          PSMs<-PSMs%>% dplyr::select(sample_id,sample_name) %>% unique(.)
        }
        df3<-dplyr::bind_rows(df3)
        df2<-df3 %>% dplyr::right_join(PSMs,by="sample_id")
        df2<-df2 %>% dplyr::group_by(Accession,sample_id) %>%  dplyr::mutate(missing_pct=(100*sum(as.numeric(is.na(value)))/length(value)))%>% ungroup(.)
        
        
      }
      
      
      check_PD<-function(x) {
        if(any("Grouped" %in% names(x))==TRUE){
          warning(paste0("please check sample file ",x$sample_name," for Abundance columns")) 
        }
      }
      ch<-mclapply(df2,check_PD,mc.cores=availableCores())
    }
    
    #return(df2)
  }else{ #if this is a protein file
    i<-1
    #Select protein columns
    if(any(names(Proteins)=="File.ID")){
      Proteins<-dplyr::bind_rows(Proteins) %>% dplyr::select("Accession","File.ID","MW.kDa.","Biological.Process","Molecular.Function","Cellular.Component","Coverage.",tidyselect::starts_with('Abundance')|tidyselect::starts_with('1')) %>% 
        dplyr::rename("Cell_Component"="Cellular.Component","MW_kDa"="MW.kDa.","Bio_Process"="Biological.Process","Coverage"="Coverage.")
    }else{
      Proteins<-dplyr::bind_rows(Proteins) %>% dplyr::select("Accession","MW.kDa.","Biological.Process","Molecular.Function","Cellular.Component","Coverage.",tidyselect::starts_with('Abundance')|tidyselect::starts_with('1')) %>% 
        dplyr::rename("Cell_Component"="Cellular.Component","MW_kDa"="MW.kDa.","Bio_Process"="Biological.Process","Coverage"="Coverage.")
      
    }
    #Wide to long for proteins
    df2<-data.table(dplyr::bind_rows(Proteins))
    df3<-melt(data = df2, 
              id.vars = c("Accession","MW_kDa","Cell_Component","Bio_Process","Coverage"),
              measure.vars = names(df2)[grepl( "Abundance" , names( df2 ) )],
              variable.name = "id",
              value.name = "value")
    df3<-df3[, id:=as.character(id)]
    df3<-data.frame(df3) %>% dplyr::mutate(sample_id=stringr::str_extract(id,"F[[:digit:]]+"),
                                           temp_ref = stringr::str_extract(id,find))
    
    #Wide to long for peptides
    if(any(names(dplyr::bind_rows(PSMs))=="File.ID")){
      dfP<-data.table(dplyr::bind_rows(PSMs) %>% dplyr::rename("Accession"="Protein.Accessions","sample_name"="Spectrum.File","sample_id"="File.ID"))
    }else{
      dfP<-data.table(dplyr::bind_rows(PSMs) %>% dplyr::rename("Accession"="Protein.Accessions","sample_name"="Spectrum.File","sample_id"="sample"))
    }
    dfP<-melt(data = dfP, 
              id.vars = c("Accession","sample_name","sample_id"),
              measure.vars = names(dfP)[grepl( "Abundance" , names( dfP ) )],
              variable.name = "id",
              value.name = "value")
    dfP<-dfP[, id:=as.character(id)]
    dfP<-data.frame(dfP) %>% dplyr::mutate(temp_ref=stringr::str_extract(id,find),
                                           sample_id=stringr::str_extract(sample_id,"F[[:digit:]]+")) %>%
      dplyr::select(-id,-value)
    #Join protein and PSM file
    Protein<-df3 %>% dplyr::right_join(dfP,by=c("Accession","sample_id","temp_ref"))
    Protein<-Protein %>% dplyr::mutate(sample_name=str_remove(sample_name,"[[:digit]]+[:punct:][[:digit]]+.raw"))
    Protein<-Protein %>% dplyr::group_split(Accession,sample_id,sample_name)
    Protein<-purrr::map(Protein,function(x){ x %>%  dplyr::mutate(missing=is.na(value),missing_pct=100*sum(is.na(value),na.rm=TRUE)/length(value))})
    
    df_raw_D_R<-dplyr::bind_rows(Protein) %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(value,3))%>% dplyr::select(sample_id,Accession,rank)
    df2<-dplyr::bind_rows(Protein) %>% dplyr::right_join(df_raw_D_R,by=c('sample_id','Accession'))
    df2<-df2 %>% dplyr::group_by(Accession,sample_id)
    df2<-df2 %>%  dplyr::mutate(missing_pct=100*sum(as.numeric(is.na(value)))/length(value)) %>% ungroup(.)
    
  }
  
  
  if(isTRUE(Batch)){
    
    i=1
    file.list<-list.files(protein_path,pattern=Prot_Pattern)
    df<-vector("list",length(file.list))
    df[[i]]<-data.frame()
    for( i in seq(file.list)){
      df.raw <- readxl::read_excel(file.list[i])
      
      check<-any(stringr::str_detect(names(df.raw),'Grouped'))
      
      if (isTRUE(check)){
        df[[i]] <- df.raw %>%
          dplyr::select(Accession,tidyselect::starts_with('Abundance'),-tidyselect::contains('Grouped')) %>%
          tidyr::gather('id', 'value', -Accession) %>%
          dplyr::mutate(sample_id = as.factor(paste("F",as.character(i),sep="")),
                        temp_ref = stringr::str_extract_all(id,find),
                        missing=ifelse(is.na(value),1,0))
      }
    }
    df<-lapply(df,function(x) x %>% tidyr::unnest(cols="sample_id"))
    df<-dplyr::bind_rows(df)
    
    df_raw_D_R<-df %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(value,3))%>% dplyr::select(sample_id,Accession,rank)
    df<-df %>% dplyr::left_join(df_raw_D_R,by=c('sample_id','Accession'))
    df2<-df %>% dplyr::group_by(Accession,sample_id)
    df2<-df2 %>%  dplyr::mutate(missing_pct=100*sum(as.numeric(is.na(value)))/length(value)) %>% ungroup(.)
    
    return(df2)
    
  }else{
    if(any(stringr::str_detect(df2$sample_name,solvent)==TRUE) & isTRUE(CFS)){#if the preset solvent is present under sample name
      df2$CC<-ifelse(stringr::str_detect(df2$sample_name,solvent)==TRUE,0,1)
      df2$dataset<-ifelse(df2$CC==0,"vehicle","treated")
    }else{#if doing a time_point experiment
      df2$CC<-as.numeric(df2$sample_name)
      
    }
    if(isTRUE(CFS)){
      df2$sample_name<-paste0(ifelse(str_detect(df2$sample_name,"NOcarrier")==TRUE,"nC",ifelse(str_detect(df2$sample_name,"carrier")==TRUE,"C",NA)),'_',
                              ifelse(str_detect(df2$sample_name,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(df2$sample_name,"r_FAIMS")==TRUE,"F",NA)),'_',
                              ifelse(str_detect(df2$sample_name,"S_eFT")==TRUE,"E",ifelse(str_detect(df2$sample_name,"S_Phi")==TRUE,"S",NA)))#oncentration values are defined in uM
    }
    return(df2)
    
  }
  
  
}
#' Clean data
#'
#' Clean CETSA data
#'
#' @param df.  Data frame returned by read_cetsa
#' @param temperatures.  Data frame of temperatures related to TMT tags - columns =
#'     temp_ref, temperature
#' @param samples.  Optional data frame of sample names - columns = sample, sample_name
#' @param separator.  Character used to separate parts of sample name.  Sample name is constructed as
#'     sampleroot_A_B where sampleroot is the name of the sample, A is the number of the biological
#'     replicate and B is the number of the technical replicate.
#'
#' @return a data frame of clean data
#'
#' @import dplyr
#'
#' @export
clean_cetsa <- function(df, temperatures = NULL, samples = NULL,Peptide=FALSE,solvent,CFS=TRUE,CARRIER=TRUE){
  if(isTRUE(CARRIER)){
    df<-df %>% dplyr::filter(!temp_ref=="131C")
  }
  if (is.null(temperatures)) {return(warning('No temperature data'))}
  if (is.null(samples)) {return(warning('No sample data'))}
  if(any(names(df)=="I")){
    df<-df %>% dplyr::rename("value"="I")
  }
  if(any(names(df)=="uniqueID")){
    df<-df %>% dplyr::rename("Accession"="uniqueID")
  }
  if(any(names(df)=="sample")){
    df<-df %>% dplyr::rename("sample_id"="sample")
  }
  df <- df %>% dplyr::ungroup(.)
  if(any(ncol(df.temps)>2)){
    if(any(names(df)=="dataset")){
      df <- df %>% dplyr::select(-dataset,-sample_name)
    }else{
      df <- df %>% dplyr::select(-sample_name)
    }
    col_n<-dplyr::intersect(names(df),names(temperatures))
    df<-data.table(df)
    temperatures<-data.table(temperatures)
    setkeyv(df,cols=col_n)
    setkeyv(temperatures,cols=col_n)
    #right_join
    df_ <- data.frame(merge(df,temperatures, all.y=TRUE))
    
    df$sample_name<-df$dataset
    df<-df %>% dplyr::mutate(dataset=ifelse(df$sample_name==0,"vehicle","treated"))
  }else if(!is.na(samples)){
    df <- df %>%
      dplyr::left_join(samples, by = intersect(names(df),names(samples)))
  }else{#if samples data is missing
    df <- df %>%
      dplyr::right_join(temperatures, by = intersect(names(df),names(temperatures)))
  }
  
  
  if(isTRUE(Peptide) & any(str_detect(names(df),"Fraction")=="TRUE")){#if this is a peptide file and is fractionated
    df <- df %>%
      dplyr::group_split(Accession, sample_id, temperature)
    df<-purrr::map(df,function(x) x %>% dplyr::mutate(value=sum(value,na.rm=TRUE)) %>% ungroup(.))
    df<-dplyr::bind_rows(df)
    
    df <- df %>%
      dplyr::select(Accession, sample_id, temperature,value,rank,Fraction,missing,missing_pct,sample_name) %>% 
      dplyr::filter(!is.na(temperature),!is.na(value)) %>%
      dplyr::group_by(Accession,sample_id) %>%
      dplyr::mutate(rank=rank,missing=missing,value = value / value[temperature == min(temperature,na.rm=TRUE)]) %>% distinct(.) %>% ungroup(.)
    
    if(any(names(df)=="sample_id" & !any(names(df)=="sample"))){
      df<-df %>% dplyr::rename("sample"="sample_id")
    }
    df<-df[!is.na(df$Accession),]
    
  }else if(isTRUE(Peptide)){#if this is an unfractionated peptide file
    if(any(names(df)=="sample" & !any(names(df)=="sample_id"))){
      df<-df %>% dplyr::rename("sample_id"="sample")
    }
    df <- df %>%
      #dplyr::select(Accession, sample_id, temperature,value,rank,missing,missing_pct,sample_name) %>% 
      dplyr::filter(!is.na(temperature),!is.na(value)) %>%
      dplyr::group_split(Accession,sample_id)
    df<-purrr::map(df,function(x) x%>%
                     dplyr::mutate(value = value / value[temperature == min(temperature,na.rm=TRUE)]) %>% distinct(.) %>% ungroup(.))
    df<-dplyr::bind_rows(df)
    if(any(names(df)=="sample_id" & !any(names(df)=="sample"))){
      df<-df %>% dplyr::rename("sample"="sample_id")
    }
    df<-df[!is.na(df$Accession),]
    
  }else{#if this is a protein file
    
    if(any(names(df)=="sample") & !isTRUE(any(names(df)=="sample_id"))){
      df<-df %>% dplyr::rename("sample_id"="sample")
    }
    df <- df %>%
      #dplyr::select(Accession, sample_id, temperature,value,rank,missing,missing_pct,sample_name) %>% 
      dplyr::filter(!is.na(temperature),!is.na(value)) %>%
      dplyr::group_by(Accession,sample_id) %>%
      dplyr::mutate(rank=rank,missing=missing,value = value / value[temperature == min(temperature,na.rm=TRUE)]) %>% distinct(.) %>%  ungroup(.)
    
    if(any(names(df)=="sample_id") & !isTRUE(any(names(df)=="sample"))){
      df<-df %>% dplyr::rename("sample"="sample_id")
    }
    df<-df[!is.na(df$Accession),]
    
  }
  
  if(any(names(df)=="uniqueID")){
    df<-df %>% dplyr::rename("Accession"="uniqueID")
  }
  df<-df %>% dplyr::rename("I"="value")
  return(df)
  
}
# if(any(names(df)=="sample")){
#   df<-df %>% dplyr::rename("uniqueID"="Accession","sample_id"="sample")
# }
# df$I<-as.numeric(df$I)
# df <- df %>%
#   #dplyr::select(Accession, sample_id, temperature,value,rank,missing,missing_pct,sample_name) %>% 
#   dplyr::filter(!is.na(temperature),!is.na(I)) %>%
#   dplyr::group_by(uniqueID,sample_id) %>%
#   dplyr::mutate(I = I / I[temperature == min(temperature,na.rm=TRUE)]) %>% unique(.) %>% 
#   dplyr::rename("sample"="sample_id")
# df<-df[!is.na(df$uniqueID),]
#only for proteins
#   if(isTRUE(CFS)&any(str_detect(df$sample_name,solvent)==TRUE)){
#     df$CC<-ifelse(stringr::str_detect(df$sample_name,solvent)==TRUE,0,1)
#     df$sample_name<-paste0(ifelse(str_detect(df$sample_name,"NOcarrier")==TRUE,"nC",ifelse(str_detect(df$sample_name,"carrier")==TRUE,"C",NA)),'_',
#                            ifelse(str_detect(df$sample_name,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(df$sample_name,"r_FAIMS")==TRUE,"F",NA)),'_',
#                            ifelse(str_detect(df$sample_name,"S_eFT")==TRUE,"E",ifelse(str_detect(df$sample_name,"S_Phi")==TRUE,"S",NA)))#oncentration values are defined in uM
#     
#     
#     df$dataset<-ifelse(df$CC==0,"vehicle","treated")
#   }
#   return(df)
#   
# }


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
normalize_cetsa <- function(df, temperatures,Peptide=FALSE,filters=FALSE,CARRIER=TRUE) {
  if(isTRUE(CARRIER)){
    df<-df %>% dplyr::filter(!temp_ref=="131C")
  }
  if(isTRUE(Peptide)){
    temperatures <- sort(temperatures)
    if(any(names(df)=="uniqueID")){
      df<-df %>% dplyr::rename("Accession"="uniqueID")
    }
    if(any(names(df)=="I")){
      df<-df %>% dplyr::rename("value"="I")
    }
    df$Accession<-as.factor(df$Accession)
    df$dataset<-as.factor(df$dataset)
    df$temperature<-as.numeric(as.character(df$temperature))
    #select top3 top5 and top10 peptides according to rank
    df_filt3 <- df %>%dplyr::filter(!is.na(value),temperature==37) %>%  
      dplyr::group_by(rank) %>% slice(1:(nrow(df)/4)) %>% dplyr::select(Annotated_Sequence,Accession,Modifications,sample,sample_name,dataset,rank) %>% ungroup(.)
    if(any(names(df_filt3)=="id")){
      df_filt3<-df_filt3 %>% dplyr::select(-id)
    }
    df_filt5 <- df %>%dplyr::filter(!is.na(value),temperature==37) %>%  
      arrange(desc(value)) %>% 
      group_by(rank) %>% slice(1:(nrow(df)/3)) %>%
      dplyr::select(Annotated_Sequence,Accession,Modifications,sample,sample_name,dataset,rank) %>% ungroup(.)
    if(any(names(df_filt5)=="id")){
      df_filt5<-df_filt5 %>% dplyr::select(-id)
    }
    
    df_filt10 <- df %>% dplyr::filter(!is.na(value),temperature==37) %>%  
      arrange(desc(value)) %>% 
      group_by(rank) %>% slice(1:(nrow(df)/2)) %>% dplyr::select(Annotated_Sequence,Accession,Modifications,sample,sample_name,dataset,rank) %>% ungroup(.)
    if(any(names(df_filt10)=="id")){
      df_filt10<-df_filt10 %>% dplyr::select(-id)
    }
    #preserve original dataset
    df1<-df %>%dplyr::filter(temperature<68)
    df<-df %>% group_by(sample_name)
    name<-intersect(names(df1),names(df_filt3))
    #Only keep curves with the topN values
    df3<-df1 %>% dplyr::right_join(df_filt3,by=name)
    df5<-df1 %>% dplyr::right_join(df_filt5,by=name)
    df10<-df1 %>% dplyr::right_join(df_filt10,by=name)
    #remove missing values 
    df3<-df3[!is.na(df3$value),]
    df5<-df5[!is.na(df5$value),]
    df10<-df10[!is.na(df10$value),]
    
    #Calculate fold changes acros temperatures 7, 9 and 10
    df.jointP3 <- suppressWarnings(df3 %>%
                                     dplyr::group_split(Accession,sample) %>% 
                                     purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                                  dplyr::mutate(.,T7 = try(mean(value[temperature == temperatures[7]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                T9 = try(mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                T10 = try(mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]],na.rm=TRUE)))))
    
    df.jointP5 <- suppressWarnings(df5 %>%
                                     dplyr::group_split(Accession,sample) %>% 
                                     purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                                  dplyr::mutate(T7 = try(mean(x$value[temperature == temperatures[7]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                T9 = try(mean(x$value[temperature == temperatures[9]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                T10 = try(mean(x$value[temperature == temperatures[10]]/value[temperature == temperatures[1]],na.rm=TRUE)))))
    
    df.jointP10 <- suppressWarnings(df10 %>%
                                      dplyr::group_split(Accession,sample) %>% 
                                      purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                                   dplyr::mutate(.,T7 = try(mean(value[temperature == temperatures[7]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                 T9 = try(mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                                 T10 = try(mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]],na.rm=TRUE)))))
    
    df.jointP3<- dplyr::bind_rows(df.jointP3)
    df.jointP5<- dplyr::bind_rows(df.jointP5)
    df.jointP10<- dplyr::bind_rows(df.jointP10)
    
    if(isTRUE(filters)){
      #top3
      df.jointP3<-df.jointP3 %>% dplyr::filter(T7 >= 0.4, T7 <= 0.6)
      df.jointP3<-df.jointP3 %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
      if(any(names(df.jointP3)=="T10" & all(!is.na(df.jointP3$T10)))){
        df.jointP3<- df.jointP3 %>% subset(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
      }
      #top 5
      df.jointP5<-df.jointP5 %>% dplyr::filter(T7 >= 0.4 & T7 <= 0.6)
      df.jointP5<-df.jointP5 %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
      if(any(names(df.jointP5)=="T10"& all(!is.na(df.jointP5$T10)))){
        df.jointP5<- df.jointP5 %>% dplyr::filter(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
      }
      #top 10
      df.jointP10<-df.jointP10 %>% dplyr::filter(T7 >= 0.4 & T7 <= 0.6)
      df.jointP10<-df.jointP10 %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
      if(any(names(df.jointP10)=="T10")& all(!is.na(df.jointP3$T10))){
        df.jointP10<- df.jointP10 %>% dplyr::filter(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
      }
      #also filter original data
      df1<-dplyr::bind_rows(df1)
      df1 <- suppressWarnings(df1 %>% dplyr::group_split(Accession,sample) %>% 
                                purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                             dplyr::mutate(.,T7 = try(mean(value[temperature == temperatures[7]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                           T9 = try(mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                           T10 = try(mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]],na.rm=TRUE)))))
      df1<-dplyr::bind_rows(df1)
      df1<-df1 %>% dplyr::filter(T7 >= 0.4 & T7 <= 0.6)
      df1<-df1 %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
      if(any(names(df1)=="T10") & all(!is.na(df1$T10))){
        df1<- df1 %>% dplyr::filter(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
      }
    }
    if(nrow(df.jointP3)==0){
      return(warning("Please disable filters, all data was filtered out for top3."))
    }
    if(nrow(df.jointP5)==0){
      return(warning("Please disable filters, all data was filtered out for top5."))
    }
    if(nrow(df.jointP10)==0){
      return(warning("Please disable filters, all data was filtered out for top10."))
    }
    
    ## split[[i]] by sample group and filter
    l.bytype3 <- split.data.frame(df.jointP3, df.jointP3$sample)
    l.bytype5 <- split.data.frame(df.jointP5, df.jointP5$sample)
    l.bytype10 <- split.data.frame(df.jointP10, df.jointP10$sample)
    
    ## determine which sample (F1 through FN) contains the greatest number of PSM curves and use this for normalization
    n.filter3 <- lapply(l.bytype3, nrow)
    n.filter5 <- lapply(l.bytype5, nrow)
    n.filter10 <- lapply(l.bytype10, nrow)
    #which sample has the greatest # of PSMs
    df.normP3 <- l.bytype3[[which.max(n.filter3)]]
    df.normP5 <- l.bytype5[[which.max(n.filter5)]]
    df.normP10<- l.bytype10[[which.max(n.filter10)]]
    #choose the accessions
    norm.accessions3 <- df.normP3$Accession
    norm.accessions5 <- df.normP5$Accession
    norm.accessions10 <- df.normP10$Accession
    ## calculate median for each sample group
    
    df.jointP3<-dplyr::bind_rows(df.jointP3)
    df.jointP5<-dplyr::bind_rows(df.jointP5)
    df.jointP10<-dplyr::bind_rows(df.jointP10)
    #calculate median values for TopN
    df.median3 <- df.jointP3 %>%
      dplyr::group_by(sample,temperature,dataset) %>%
      dplyr::mutate(value = median(value,na.rm=TRUE)) %>% ungroup(.)
    
    df.median5 <- df.jointP5 %>%
      dplyr::group_by(sample,temperature,dataset) %>%
      dplyr::mutate(value = median(value,na.rm=TRUE)) %>% ungroup(.)
    
    df.median10 <- df.jointP10 %>%
      dplyr::group_by(sample,temperature,dataset) %>%
      dplyr::mutate(value = median(value,na.rm=TRUE)) %>% ungroup(.)
    
    df.median3$temperature<-as.numeric(df.median3$temperature)
    df.median5$temperature<-as.numeric(df.median5$temperature)
    df.median10$temperature<-as.numeric(df.median10$temperature)
    
    df.median3<-df.median3[!is.na(df.median3$value),]
    df.median5<-df.median5[!is.na(df.median5$value),]
    df.median10<-df.median10[!is.na(df.median10$value),]
    ## fit curves to the median data for each sample (F1 through FN)
    
    df.median3<-df.median3 %>% dplyr::filter(temperature<68)
    df.median5<-df.median5 %>% dplyr::filter(temperature<68)
    df.median10<-df.median10 %>% dplyr::filter(temperature<68)
    #
    df.fit3 <- df.median3 %>%
      dplyr::group_by(sample,dataset) %>% 
      dplyr::mutate(fit = try(list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                           start = c(Pl=0, a = 550, b = 10),
                                           data = list(x=temperature,y=value),
                                           na.action = na.exclude,
                                           algorithm = "port",
                                           lower = c(0.0,1e-5,1e-5),
                                           upper = c(1.5,15000,300),
                                           control = nls.control(maxiter = 50)),
                                       silent = TRUE)))) 
    
    df.fit3<- df.fit3%>% 
      dplyr::mutate(fitted_values3 = ifelse(!is.logical(fit[[1]]),list(data.frame(fitted_values=predict(fit[[1]]))),NA)) %>% 
      dplyr::select(sample,fitted_values3,temperature,dataset) %>% ungroup(.)
    
    df.fit5 <- df.median5 %>%
      dplyr::group_by(sample,dataset) %>% 
      dplyr::mutate(fit = try(list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                           start = c(Pl=0, a = 550, b = 10),
                                           data = list(x=temperature,y=value),
                                           na.action = na.exclude,
                                           algorithm = "port",
                                           lower = c(0.0,1e-5,1e-5),
                                           upper = c(1.5,15000,300),
                                           control = nls.control(maxiter = 50)),
                                       silent = TRUE)))) 
    
    df.fit5<- df.fit5%>% 
      dplyr::mutate(fitted_values5 = ifelse(!is.logical(fit[[1]]),list(data.frame(fitted_values=predict(fit[[1]]))),NA)) %>% 
      dplyr::select(sample,fitted_values5,temperature,dataset) %>% ungroup(.)
    
    df.fit10 <- df.median10 %>%
      dplyr::group_by(sample,dataset) %>% 
      dplyr::mutate(fit = try(list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                           start = c(Pl=0, a = 550, b = 10),
                                           data = list(x=temperature,y=value),
                                           na.action = na.exclude,
                                           algorithm = "port",
                                           lower = c(0.0,1e-5,1e-5),
                                           upper = c(1.5,15000,300),
                                           control = nls.control(maxiter = 50)),
                                       silent = TRUE))))  
    
    df.fit10<- df.fit10%>% 
      dplyr::mutate(fitted_values10 = ifelse(!is.logical(fit[[1]]),list(data.frame(fitted_values=predict(fit[[1]]))),NA)) %>% 
      dplyr::select(sample,fitted_values10,temperature,dataset) %>% ungroup(.)
    
    ## calculate the fitted values
    # d3<-length(dplyr::bind_rows(df.fit3$fitted_values3))
    # d5<-length(dplyr::bind_rows(df.fit5$fitted_values5))
    # d10<-length(dplyr::bind_rows(df.fit10$fitted_values10))
    #unnest fitted values from list and name value column and keep fitted values and temps
    
    #unnest fitted values from list and name value column
    check3<-data.frame(fitted_values3=unique(df.fit3$fitted_values3[df.fit3$dataset=="treated"][[1]]),dataset="treated",temperature=unique(df.fit3$temperature))
    check31<-data.frame(fitted_values3=unique(df.fit3$fitted_values3[df.fit3$dataset=="vehicle"][[1]]),dataset="vehicle",temperature=unique(df.fit3$temperature))
    check3<-rbind(check3,check31)
    check_3<-df.fit3 %>% dplyr::select(-sample,-fitted_values3) %>% unique
    check3<-check_3 %>% dplyr::right_join(check3,by=c("temperature","dataset"))
    
    check5<-data.frame(fitted_values5=unique(df.fit5$fitted_values5[df.fit5$dataset=="treated"][[1]]),dataset="treated",temperature=unique(df.fit5$temperature))
    check51<-data.frame(fitted_values5=unique(df.fit5$fitted_values5[df.fit5$dataset=="vehicle"][[1]]),dataset="vehicle",temperature=unique(df.fit5$temperature))
    check5<-rbind(check5,check51)
    check_5<-df.fit5 %>% dplyr::select(-sample,-fitted_values5) %>% unique
    check5<-check_5 %>% dplyr::right_join(check5,by=c("temperature","dataset"))
    
    check10<-data.frame(fitted_values10=unique(df.fit10$fitted_values10[df.fit10$dataset=="treated"][[1]]),dataset="treated",temperature=unique(df.fit10$temperature))
    check101<-data.frame(fitted_values10=unique(df.fit10$fitted_values10[df.fit10$dataset=="vehicle"][[1]]),dataset="vehicle",temperature=unique(df.fit10$temperature))
    check10<-rbind(check10,check101)
    check_10<-df.fit10 %>% dplyr::select(-sample,-fitted_values10) %>% unique
    check10<-check_10 %>% dplyr::right_join(check10,by=c("temperature","dataset"))
    
    df1<-df1 %>% dplyr::filter(temperature<68)
    col_n<-dplyr::intersect(names(df1),names(check3))
    check3<-data.table(check3)
    check5<-data.table(check5)
    check10<-data.table(check10)
    df1<-data.table(df1)
    
    setkeyv(df1,cols=col_n)
    setkeyv(check3,cols=col_n)
    setkeyv(check5,cols=col_n)
    setkeyv(check10,cols=col_n)
    #right_join
    test3 <- data.frame(merge(df1,check3, all.y=TRUE))
    test5 <- data.frame(merge(df1,check5, all.y=TRUE))
    test10 <- data.frame(merge(df1,check10, all.y=TRUE))
    
    ## calculate ratios between the fitted curves and the median values
    df.out3 <- test3 %>% as.data.frame() %>% 
      dplyr::mutate(correction3 = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
      dplyr::select('sample','temperature','correction3','dataset')
    
    
    df.out5 <- test5 %>%as.data.frame() %>% 
      dplyr::mutate(correction5 = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
      dplyr::select('sample','temperature','correction5','dataset')
    
    
    df.out10 <- test10 %>%as.data.frame() %>% 
      dplyr::mutate(correction10 = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
      dplyr::select('sample','temperature','correction10','dataset')
    
    ## join correction factor to data
    df1$temperature<-as.factor(df1$temperature)
    df1<-df1 %>% dplyr::group_split(temperature)
    #group split data by temperature
    df.out3$temperature<-as.factor(df.out3$temperature)
    df.out3<-df.out3 %>% dplyr::group_by(temperature,dataset) %>% dplyr::group_split(.)
    
    df.out5$temperature<-as.factor(df.out5$temperature)
    df.out5<-df.out5 %>% dplyr::group_by(temperature,dataset) %>% dplyr::group_split(.)
    
    df.out10$temperature<-as.factor(df.out10$temperature)
    df.out10<-df.out10 %>% dplyr::group_by(temperature,dataset) %>% dplyr::group_split(.)
    
    df1<-dplyr::bind_rows(df1)
    df1<-df1 %>% dplyr::group_split(temperature,dataset)
    #apply correction factor by temperature to original data
    df3<-purrr::map2(df1,df.out3,function(x,y)x %>% dplyr::mutate(correction3=y$correction3))
    df5<-purrr::map2(df3,df.out5,function(x,y)x %>% dplyr::mutate(correction5=y$correction5))
    df<-purrr::map2(df5,df.out10,function(x,y)x %>% dplyr::mutate(correction10=y$correction10))
    #bind_rows
    df<-dplyr::bind_rows(df)
    
    df3 <- df%>% 
      dplyr::mutate(norm_value3= ifelse(is.na(correction3),value,value * correction3)) %>% dplyr::ungroup(.)
    df5 <- df3%>% 
      dplyr::mutate(norm_value5= ifelse(is.na(correction5),value,value * correction5)) %>% dplyr::ungroup(.)
    df<- df5%>% 
      dplyr::mutate(norm_value10= ifelse(is.na(correction10),value,value * correction10)) %>% dplyr::ungroup(.)
    
    
    df <- df %>% 
      dplyr::rename("uniqueID"="Accession", "C"="temperature","I3"="norm_value3","I5"="norm_value5","I10"="norm_value10")
    df<-df %>% dplyr::ungroup(.) %>% dplyr::select(-value,-correction3,-correction5,-correction10)
    
    df<-dplyr::bind_rows(df)
    df$sample_name<-str_replace(df$sample_name,"S","\u03A6")
    df<-df %>% distinct(.)
    return(df)
  }else{
    
    temperatures <- sort(temperatures)
    if(any(names(df)=="uniqueID")){
      df<-df %>% dplyr::rename("Accession"="uniqueID")
    }
    if(any(names(df)=="I")){
      df<-df %>% dplyr::rename("value"="I")
    }
    df$Accession<-as.factor(df$Accession)
    df$sample<-as.factor(df$sample)
    df.jointP <- suppressWarnings(df %>%
                                    dplyr::group_split(Accession,sample) %>% 
                                    purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                                 dplyr::mutate(.,T7 = try(mean(value[temperature == temperatures[7]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                               T9 = try(mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]],na.rm=TRUE)),
                                                               T10 = try(mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]],na.rm=TRUE)))))
    df.jointP<- dplyr::bind_rows(df.jointP)
    if(isTRUE(filters)){
      df.jointP<-df.jointP %>% dplyr::filter(T7 >= 0.4 & T7 <= 0.6)
      df.jointP<-df.jointP %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
      if(any(names(df.jointP)=="T10")){
        df.jointP<- df.jointP %>% dplyr::filter(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
      }
    }
    if(nrow(df)==0){
      return(warning("Please disable filters, all data was filtered out."))
    }
    
    ## split[[i]] by sample group and filter
    l.bytype <- split.data.frame(df.jointP, df.jointP$sample)
    
    ## determine which sample (F1 through FN) contains the greatest number of curves and use this for normalization
    n.filter <- lapply(l.bytype, nrow)
    df.normP <- l.bytype[[which.max(n.filter)]]
    norm.accessions <- df.normP$Accession
    
    ## calculate median for each sample group
    
    df.mynormset <- df %>% base::subset(Accession %in% norm.accessions)
    
    df.median <- df %>%
      dplyr::group_by(sample,temperature,dataset) %>%
      dplyr::mutate(value = median(value,na.rm=TRUE))
    
    
    ## fit curves to the median data for each replicate (F1 through FN)
    df.fit <- df.median %>%
      dplyr::group_by(sample,dataset) %>% 
      dplyr::mutate(fit = try(list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                           start = c(Pl=0, a = 550, b = 10),
                                           data = list(x=temperature,y=value),
                                           na.action = na.exclude,
                                           algorithm = "port",
                                           lower = c(0.0,1e-5,1e-5),
                                           upper = c(1.5,15000,300),
                                           control = nls.control(maxiter = 50)),
                                       silent = TRUE)))) 
    
    df.fit<- df.fit%>% 
      dplyr::mutate(fitted_values = ifelse(!is.logical(fit[[1]]),list(data.frame(fitted_values=predict(fit[[1]]))),NA)) %>% 
      dplyr::select(sample,dataset,fitted_values,temperature) %>% ungroup(.)
    
    ## calculate the fitted values
    check<-data.frame(fitted_values=unique(df.fit$fitted_values[df.fit$dataset=="treated"][[1]]),dataset="treated",temperature=unique(df.fit$temperature))
    check3<-data.frame(fitted_values=unique(df.fit$fitted_values[df.fit$dataset=="vehicle"][[1]]),dataset="vehicle",temperature=unique(df.fit$temperature))
    check3<-rbind(check,check3)
    check_<-df.fit %>% dplyr::select(-sample,-fitted_values) %>% unique
    check<-check_ %>% dplyr::right_join(check3,by=c("temperature","dataset"))
    
    col_n<-dplyr::intersect(names(df),names(check))
    
    check<-data.table(check)
    
    df<-data.table(df)
    
    setkeyv(df,cols=col_n)
    setkeyv(check,cols=col_n)
    
    #right_join
    test <- data.frame(merge(df,check, all.y=TRUE))
    
    
    ## calculate ratios between the fitted curves and the median values
    df.out <- test %>% as.data.frame() %>% 
      dplyr::mutate(correction = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
      dplyr::select('sample','temperature','correction','dataset')
    
    
    ## apply normalization factor to data
    
    df<-df %>% dplyr::right_join(df.out,by=intersect(names(df),names(df.out)))
    df <- df %>% 
      dplyr::mutate(norm_value = ifelse(is.na(correction),value,value * correction)) %>% dplyr::ungroup(.)
    df <- df %>% 
      dplyr::mutate(norm_value = ifelse(is.na(correction),value,value * correction)) %>% dplyr::ungroup(.)
    df <- df %>% 
      dplyr::rename("uniqueID"="Accession", "C"="temperature","I"="norm_value")
    df<-df %>% dplyr::ungroup(.) %>% dplyr::select(-value,-correction)
    if(isTRUE(filters)){
      df<-df %>% dplyr::select(-T10,-T7,-T9)
    }
    
    df<-df %>% distinct(.)
    return(df)
  }
}

#
#' fit curves to CETSA data
#'
#' fit curves to normalized or non-normalized CETSA data
#'
#' @param df. data frame containing CETSA data
#' @param normalized_data  Boolean.  If true then use normalized data otherwise use non-normalized
#' @param n_cores.  Number of cores for parallel processing
#' @param separator.  Character used to separate sample name from replicate.  If NULL then
#'     replicates are treated as individual samples.
#'
#' @return list
#'
#' @import nls2
#' @import dplyr
#' @importFrom tidyr separate
#'
#' @importFrom tibble rowid_to_column
#'
#' @export

cetsa_fit<-function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    #c(Pl=0, a = 550, b = 10
    fine_start <- expand.grid(p=c(0,0.5),k=seq(0,1000,by=100),m=seq(5,45,by=10))
    new_start <- nls2::nls2(y ~ fit.cetsa(p, k, m, t),
                            data = myData,
                            start = fine_start,
                            algorithm = "grid-search",#note: check other ones
                            control = nls.control(warnOnly=T,maxiter=5000))
    nls2::nls2(y ~ fit.cetsa(p, k, m, t),
               data = myData,
               start = new_start,
               control = nls.control(warnOnly=F),
               algorithm = "grid-search",
               lower = c(0, 1, 10),
               upper = c(0.4, 100000, 100))
  }, error = function(err) {
    return(NA)
  })
  return(result)
}


#' calculate CETSA statistics
#'
#' calculate CETSA statistics - top and bottom plateau gradient and goodness of fit
#'
#' @param df.  data frame containing curve-fitted CETSA data
#' @param plateau_temps.  Vector containing plateau values (T1A, T1B, T2A, T2B)
#'
#' @return data frame containing curves with stats
#'
#' @import dplyr
#' @importFrom tibble as.tibble
#'
#' @export
stats_cetsa <- function(df, plateau_temps = c(37,40,64,67)) {
  df <- df%>%
    mutate(gof = gof(fit))
  
  responses <- apply(df, 1, function(x) {
    tryCatch({
      predict(x['fit'][[1]], list(t = plateau_temps))
    }, error = function(e) {
      rep(NA, length(plateau_temps))
    })
  })
  df.responses <- as.data.frame(as.matrix(t(responses)))
  df.responses$ref <- as.numeric(row.names(df.responses))
  df.responses$slopeStart <- apply(df.responses, 1, function(x) abs(x[2] - x[1]) / (plateau_temps[2] - plateau_temps[1]))
  df.responses$slopeEnd <- apply(df.responses, 1, function(x) abs(x[4] - x[3]) / (plateau_temps[4] - plateau_temps[3]))
  
  df <- cbind(df, slope_start=df.responses$slopeStart, slope_end=df.responses$slopeEnd)
  return(as.tibble(df))
}

#' tag CETSA curves
#'
#' Apply filters and tag curves
#'
#' @param df.  data frame containing curve-fitted CETSA data
#' @param tag_fit.  Boolean.  Tag if poor curve fit
#' @param tag_range.  Boolean.  Tag if Tm is outside experimental temperature range
#' @param tag_gof.  Boolean.  Tag if poor goodness-of-fit
#' @param tag_plateau.  Boolean.  Tag if plateau(s) not flat
#' @param plateau_filter.  Filter value for tag_plateau (default = 0.01)
#' @param gof_filter.  Filter value for tag_gof (default = 0.085)
#'
#' @return data frame containing tagged results
#'
#' @import dplyr
#' @importFrom purrr map
#'
#' @export
tag_cetsa <- function(df, tag_fit = TRUE, tag_range = TRUE, tag_gof = TRUE, tag_plateau = TRUE, plateau_filter = 0.01, gof_filter = 0.085) {
  
  df <- df %>%dplyr::ungroup()
  
  if (tag_fit) {
    df <- df %>%
      mutate(tag_fit = if_else(is.na(fit), 'bad fit', NA_character_))
  }
  
  if (tag_range) {
    df <- df %>%
      mutate(temps = map(df$fit, c('data', 't'))) %>%
      rowwise() %>%
      mutate(t_min = min(temps, na.rm=T)) %>%
      mutate(t_max = max(temps, na.rm=T)) %>%
      ungroup() %>%
      mutate(tag_range = if_else(is.na(fit), 'no data',
                                 if_else(Tm < t_min | Tm > t_max, 'outside experimental range', NA_character_))) %>%
      select (-temps, -t_min, -t_max)
  }
  
  if (tag_gof) {
    df <- df %>%
      mutate(tag_gof = if_else(gof >= gof_filter, 'poor gof', NA_character_))
  }
  
  if (tag_plateau) {
    df <- df %>%
      rowwise() %>%
      mutate(tag_plateau = if_else(min(slope_start, slope_end) >= plateau_filter, 'poor plateau - low and high',
                                   if_else(slope_start >= plateau_filter, 'poor plateau - low',
                                           if_else(slope_end >= plateau_filter, 'poor plateau - high', NA_character_))))
  }
  
  ## add tag_summary column
  use_cols <- names(df)[which(startsWith(names(df), 'tag_'))]
  if (length(use_cols) > 0) {
    df$tag_summary <- apply(df, 1, function(x) {
      cols <- sapply(use_cols, function(y) x[y])
      cols <- cols[!is.na(cols)]
      paste0(cols, collapse = '; ')
    })
    
    df$tagged <- nchar(df$tag_summary) > 0
  }
  
  return(df)
}


#' extract curve fit parameters from a CETSA data frame
#'
#' extract curve fit parameters from a CETSA data frame
#'
#' @param df.  Data frame containing curve fit column
#' @param fit.  Curve fit column name (default = "fit")
#' @param remove_fit.  Boolean.  If true then remove the fit column
#'
#' @import dplyr
#'
#' @export
extract_parameters <- function(df, fit_col = 'fit', remove_fit = TRUE) {
  fit_col_ref <- which(names(df) == fit_col)
  if (length(fit_col_ref) == 0) {
    message('fit column not found')
    return(NULL)
  }
  
  # get variable names
  fit_types <- sapply(df[[fit_col]], class)
  first_match <- which(!fit_types == 'logical')[1]
  param_names <- names(df[[fit_col]][[first_match]]$m$getPars())
  
  df <- df %>%
    rowwise() %>%
    mutate(params = curve_parameters(!!as.name(fit_col)))
  
  # rearrange column order
  if (remove_fit) {
    df <- df[, c(1:(fit_col_ref-1), ncol(df), (fit_col_ref+1):(ncol(df)-1))]
  } else {
    df <- df[, c(1:fit_col_ref, ncol(df), (fit_col_ref+1):(ncol(df)-1))]
  }
  
  return(df)
}


#' convert to wide format
#'
#' convert to wide format
#'
#' @param df. Data frame containing Accession, sample and params columns
#'
#' @return a wide data frame
#'
#' @import dplyr
#' @importFrom tidyr spread
#'
#' @export
long_to_wide_cetsa <- function(df) {
  df %>%
    dplyr::select(Accession, sample, params) %>%
    spread(sample, params)
}


#' curve fitting equation
#'
#' CETSA curve fitting equation
#' Curve equation \eqn{y = \frac{(1 - p)}{(1 + e^{(-k(\frac{1}{t} - \frac{1}{m}))})} + p}
#'
#' @param p.  Curve parameter
#' @param k.  Curve parameter
#' @param m.  Curve parameter
#' @param t.  Curve variable
#'
#' @return position on curve at t
#'
fit.cetsa <- function(p, k, m, t) {
  (1 - p)/(1 + exp(-k*(1/t - 1/m))) + p
}


#' curve fitting function
#'
#' CETSA curve fitting function
#' Fit a curve to a set of data
#'
#' @param d.  data frame containing temperature and value or norm_value columns
#' @param norm.  Boolean.  If true then use norm_value column otherwise use value column
#'
#' @return nls2 curve model
#'
#' @import nls2
#'
## fit curve using mstherm equation
cetsa_fit <- function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    #c(Pl=0, a = 550, b = 10
    fine_start <- expand.grid(p=c(0,0.01),k=seq(500,600),m=seq(5,15,by=10))
    new_start <- nls2::nls2(y ~ fit.cetsa(p, k, m, t),
                            data = myData,
                            start = fine_start,
                            algorithm = "grid-search",#note: check other ones
                            control = nls.control(warnOnly=T,maxiter=5000))
    nls2::nls2(y ~ fit.cetsa(p, k, m, t),
               data = myData,
               start = new_start,
               control = nls.control(warnOnly=F),
               algorithm = "grid-search",
               lower = c(0, 1, 10),
               upper = c(0.4, 100000, 100))
  }, error = function(err) {
    return(NA)
  })
  return(result)
}


#' return Tm from CETSA curve
#'
#' return Tm from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return Tm.
#'
Tm <- function(f) {
  if (length(f) == 1) return(NA)
  pars<-f$fit[[1]]$m$getPars()
  Tm<-pars['a']/(pars['b'] - log(1-pars['Pl'])/(1/2 - pars['Pl']-1))
  return(round(Tm,1)[['a']]) 
}


#' return gof from CETSA curve
#'
#' return gof from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return gof.
#'
gof <- function(f) {
  if (length(f) == 1) return(NA)
  sigma(f)
}


#' return parameters from CETSA curve
#'
#' return parameters from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return json formatted string of parameters
#'
#' @importFrom jsonlite toJSON
#'
curve_parameters <- function(f) {
  if (length(f) == 1) return(NA)
  as.character(toJSON(c(f$data, params=list(as.list(f$m$getPars())))))
}


#' plot a series of CETSA curves
#'
#' plot a series of CETSA curves
#'
#' @param r. Subset of rows from CETSA table
#'
#' @import RColorBrewer
#' @importFrom jsonlite fromJSON
#'
#' @export
plot_cetsa <- function(r) {
  ## plot a series of curves
  cols <- brewer.pal(length(r)-1, 'Dark2')
  legend_inc <- c()
  firstplot <- TRUE
  for (c in 2:length(r)) {
    if (!is.na(r[[c]])) {
      params <- fromJSON(r[[c]])
      
      df <- data.frame(t = params$t, y = params$y) %>%
        rlist::list.group(t) %>%
        dplyr::summarise(av = mean(y, na.rm = T), sd = ifelse(dplyr::n() == 1, 0, sd(y, na.rm = T)))
      
      ## plot experimental points
      if (firstplot) {
        plot(x = df$t, y = df$av, col = cols[c-1], main = r[[1]], xlab="Temperature", ylab="Normalized Response", cex.axis = 0.8, cex.main = 0.8, cex = 0.5, ylim = c(0, 1.5))
        firstplot <- FALSE
      } else {
        points(x = params$t, y = params$y, col = cols[c-1], cex = 0.5)
      }
      
      ## error bars
      if(!all(df$sd == 0)) {
        arrows(df$t, df$av-df$sd, df$t, df$av+df$sd, col = cols[c-1], length=0.05, angle=90, code=3)
      }
      
      ## plot fitted curve
      t_range <- seq(range(params$t)[1], range(params$t)[2], by = 0.2)
      y_pred <- fit.cetsa(p=params$params$p, k=params$params$k, m=params$params$m, t_range)
      lines(t_range, y_pred, col = cols[c-1])
      
      ## add line through Tm
      abline(v = params$params$m, col = cols[c-1], lty = 2, lwd = 0.5)
      
      ## add to legend
      legend_inc <- c(legend_inc, c)
    }
  }
  ## include legend
  if (!firstplot) {
    legend("bottomleft",
           legend = sapply(names(r)[legend_inc], function(x) {paste0(x, ': ', round(params$params$m, 2), '°C')}),
           pch = rep(1, length(legend_inc)),
           col = cols[legend_inc - 1],
           cex = 0.4)
  }
}
#find patterns in data
find_pat = function(pat, x) 
{
  ff = function(.pat, .x, acc = if(length(.pat)) seq_along(.x) else integer(0L)) {
    if(!length(.pat)) return(acc)
    
    if(is.na(.pat[[1L]])) 
      Recall(.pat[-1L], .x, acc[which(is.na(.x[acc]))] + 1L)
    else 
      Recall(.pat[-1L], .x, acc[which(.pat[[1L]] == .x[acc])] + 1L)
  }
  
  return(ff(pat, x) - length(pat))
}  
###Upset plots#################
#To generate information on missing values
################################
upPSM_SN<-function(df_){
  
  df_ <-df_ %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum_File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
  
  df_$dataset<-ifelse(df_$CC==0,"vehicle","treated")
  
  
  df_$sample_name<-paste0(ifelse(str_detect(df_$Spectrum_File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(df_$Spectrum_File,"carrier")==TRUE,"C",NA)),'_',
                          ifelse(str_detect(df_$Spectrum_File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(df_$Spectrum_File,"r_FAIMS")==TRUE,"F",NA)),'_',
                          ifelse(str_detect(df_$Spectrum_File,"S_eFT")==TRUE,"E",ifelse(str_detect(df_$Spectrum_File,"S_Phi")==TRUE,"S",NA)))
  df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                            "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                            "I_Interference"="Isolation_Interference_[%]")
  
  
  #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
  df_1<-dplyr::bind_rows(df_)
  df_1$uniqueID<-as.factor(df_1$uniqueID)
  df_1$dataset<-as.factor(df_1$dataset)
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1<-df_1%>% 
    dplyr::group_split(uniqueID,sample_id,Annotated_Sequence)
  df_1<-purrr::map(df_1,function(x) x %>% dplyr::mutate(missing_pct=(sum(100*is.na(x$I))/length(x$I))) %>% head(.,1)) 
  df_1<-dplyr::bind_rows(df_1)
  
  rank<-df_1 %>% dplyr::filter(C=="126") %>% dplyr::mutate(rank=dplyr::ntile(I,3)) %>% 
    dplyr::select("uniqueID","dataset","sample_name","rank","Annotated_Sequence")
  df_1<-df_1 %>% dplyr::right_join(rank,by=c("uniqueID","dataset","sample_name","Annotated_Sequence"))
  df_1$SNrank<-dplyr::ntile(df_1$S_N,3)
  minSNrank<-min(df_1[df_1$SNrank==3,'S_N'],na.rm=TRUE)
  maxSNrankL<-max(df_1[df_1$SNrank==1,'S_N'],na.rm=TRUE)
  low<-paste0("low"," < ",as.character(maxSNrankL))
  medium<-"medium"
  high<-paste0("high"," > ",as.character(minSNrank))
  df_1<-df_1 %>% dplyr::mutate(rank=ifelse(df_1$rank==1,"low",ifelse(df_1$rank==2,"medium","high")),
                               SNrank=ifelse(df_1$SNrank==3,"high",ifelse(df_1$SNrank==2,"medium","low")))
  
  df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
  df_1$SNrank<-factor(df_1$SNrank,levels=c("high","medium","low"))
  df_1$missing_pct<-round(df_1$missing_pct,1)
  df_1$DeltaM<-round(df_1$DeltaM,1)
  df_1<-df_1%>% 
    dplyr::group_split(sample_id)
  df_1<-purrr::map(df_1,function(x)
    pivot_wider(x,names_from=c(Charge),values_from=Charge,values_fill=NA) %>%
      dplyr::select(-uniqueID,-Spectrum_File,-Annotated_Sequence,-IonInjTime,-MissedCleavages,-S_N,-sample_id,-dataset,-C,-CC,-DeltaM,-PEP,-Modifications,-I_Interference,-missing_pct,-XCorr,-I,-Protein_value))
  
  df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(SNrank,sample_name))
  df_1<-purrr::map(df_1,function(x) 1+x[,sapply(x,class)=="numeric"])
  df<-lapply(df_1,function(x) colnames(x))
  df<-lapply(df,function(x) str_replace(x,x,
                                        paste0("Charge: ","+",x, sep = " ")))
  
  
  
  
  df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
  df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
  df_1<-purrr::map(df_1,function(x)x[!is.na(x$SNrank),])
  
  
  rating_scale = scale_fill_manual(name="Ranked S/N",
                                   values=c("low" ='#fee6ce', "medium" ='#fdae6b', "high"  = '#e6550d'))
  
  show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
  UPD<-purrr::map(df_1,function(x) upset_data(x,names(x),
                                              min_degree=1,
                                              n_intersections=N,
                                              with_sizes %>% dplyr::select(order, intersection) %>%
                                                distinct() %>% arrange(order) %>% pull("intersection"))
                  
  )
  check<-list()
  check<- purrr::map(UPD,function(x)upset(x,names(x),
                                          min_degree=1,
                                          sort_intersections="ascending",
                                          set_sizes=FALSE,
                                          guides='collect',
                                          n_intersections=N,
                                          height_ratio = 0.7,
                                          stripes='white',
                                          base_annotations=list(
                                            '# of PSMs'=intersection_size(
                                              counts=TRUE,
                                            )
                                          ),
                                          annotations =list(
                                            "Bin %"=list(
                                              aes=aes(x=intersection, fill=x$SNrank),
                                              
                                              geom=list(
                                                geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                
                                                geom_text(
                                                  aes(
                                                    label=!!aes_percentage(relative_to='intersection'),
                                                    color=ifelse(!is.na(x$SNrank), 'show', 'hide')
                                                  ),
                                                  stat='count',
                                                  position=position_fill(vjust = .5),
                                                  
                                                ),
                                                
                                                scale_y_continuous(labels=scales::percent_format()),
                                                show_hide_scale,
                                                rating_scale
                                                
                                              )
                                            )
                                          )
                                          
  )+ggtitle(str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
  check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                set_sizes=FALSE,
                guides='collect',
                n_intersections=N,
                
                stripes='white',
                base_annotations=list(
                  '# of PSMs'=intersection_size(
                    counts=TRUE,
                  )
                ),
                annotations =list(
                  "Ranked Intensity  "=list(
                    aes=aes(x=intersection, fill=df_1[[1]]$SNrank),
                    
                    geom=list(
                      geom_bar(stat='count', position='fill', na.rm=TRUE),
                      themes=theme(legend.position="bottom", legend.box = "horizontal"),
                      geom_text(
                        aes(
                          label=!!aes_percentage(relative_to='intersection'),
                          color=ifelse(!is.na(df_1[[1]]$SNrank), 'show', 'hide')
                        ),
                        stat='count',
                        position=position_fill(vjust = .5),
                        
                      ),
                      
                      scale_y_continuous(labels=scales::percent_format()),
                      show_hide_scale,
                      rating_scale
                      
                    )
                  )
                )
                
  )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
  y<-get_legend(check1$patches$plots[[1]])
  data<-unlist(lapply(check,function(x) x$labels$title))
  check<-check[order(data)]
  P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
  
  return(print(P))
  
}
#listUP<-upPSM_SN(df_raw)
#df_ is the input from df_norm which includes missing percentages
#returns segmented data for vehicle and treated
upMV <- function(df_,condition,N,plot_multiple=FALSE,PSMs=FALSE){
  
  if(isTRUE(plot_multiple)){
    
    if(isTRUE(PSMs)){
      #df_<- df_ %>% dplyr::right_join(df.samples,by="sample_id")
      df_ <-df_ %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum_File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
      
      df_$dataset<-ifelse(df_$CC==0,"vehicle","treated")
      
      
      df_$sample_name<-paste0(ifelse(str_detect(df_$Spectrum_File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(df_$Spectrum_File,"carrier")==TRUE,"C",NA)),'_',
                              ifelse(str_detect(df_$Spectrum_File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(df_$Spectrum_File,"r_FAIMS")==TRUE,"F",NA)),'_',
                              ifelse(str_detect(df_$Spectrum_File,"S_eFT")==TRUE,"E",ifelse(str_detect(df_$Spectrum_File,"S_Phi")==TRUE,"S",NA)))
      df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                "I_Interference"="Isolation_Interference_[%]")
      
      
      #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
      df_1<-dplyr::bind_rows(df_)
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$dataset<-as.factor(df_1$dataset)
      df_1$sample_name<-as.factor(df_1$sample_name)
      df_1<-df_1%>% 
        dplyr::group_split(uniqueID,sample_id,Annotated_Sequence)
      df_1<-purrr::map(df_1,function(x) x %>% dplyr::mutate(missing_pct=(sum(100*is.na(x$I))/length(x$I))) %>% head(.,1)) 
      df_1<-dplyr::bind_rows(df_1)
      
      rank<-df_1 %>% dplyr::filter(C=="126") %>% dplyr::mutate(rank=dplyr::ntile(I,3)) %>% 
        dplyr::select("uniqueID","dataset","sample_name","rank","Annotated_Sequence")
      df_1<-df_1 %>% dplyr::right_join(rank,by=c("uniqueID","dataset","sample_name","Annotated_Sequence"))
      df_1$SNrank<-dplyr::ntile(df_1$S_N,3)
      minSNrank<-min(df_1[df_1$SNrank==3,'S_N'],na.rm=TRUE)
      maxSNrankL<-max(df_1[df_1$SNrank==1,'S_N'],na.rm=TRUE)
      low<-paste0("low"," < ",as.character(maxSNrankL))
      medium<-"medium"
      high<-paste0("high"," > ",as.character(minSNrank))
      df_1<-df_1 %>% dplyr::mutate(rank=ifelse(df_1$rank==1,"low",ifelse(df_1$rank==2,"medium","high")),
                                   SNrank=ifelse(df_1$SNrank==3,"high",ifelse(df_1$SNrank==2,"medium","low")))
      
      df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
      df_1$SNrank<-factor(df_1$SNrank,levels=c("high","medium","low"))
      df_1$missing_pct<-round(df_1$missing_pct,1)
      df_1$DeltaM<-round(df_1$DeltaM,1)
      df_1<-df_1%>% 
        dplyr::group_split(sample_id)
      df_1<-purrr::map(df_1,function(x)
        pivot_wider(x,names_from=c(Charge),values_from=Charge,values_fill=NA) %>%
          dplyr::select(-uniqueID,-Spectrum_File,-Annotated_Sequence,-IonInjTime,-MissedCleavages,-S_N,-sample_id,-dataset,-C,-CC,-DeltaM,-PEP,-Modifications,-I_Interference,-missing_pct,-XCorr,-I,-Protein_value))
      
      df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(SNrank,sample_name))
      df_1<-purrr::map(df_1,function(x) 1+x[,sapply(x,class)=="numeric"])
      df<-lapply(df_1,function(x) colnames(x))
      df<-lapply(df,function(x) str_replace(x,x,
                                            paste0("Charge: ","+",x, sep = " ")))
      
      
      
      
      df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
      df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
      df_1<-purrr::map(df_1,function(x)x[!is.na(x$SNrank),])
      
      
      rating_scale = scale_fill_manual(name="Ranked S/N",
                                       values=c("low" ='#fee6ce', "medium" ='#fdae6b', "high"  = '#e6550d'))
      
      show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
      
      check<-list()
      check<- purrr::map(df_1,function(x)upset(x,names(x),
                                               min_degree=1,
                                               set_sizes=FALSE,
                                               guides='collect',
                                               n_intersections=N,
                                               height_ratio = 0.7,
                                               stripes='white',
                                               base_annotations=list(
                                                 '# of PSMs'=intersection_size(
                                                   counts=TRUE,
                                                 )
                                               ),
                                               annotations =list(
                                                 "Bin %"=list(
                                                   aes=aes(x=intersection, fill=x$SNrank),
                                                   
                                                   geom=list(
                                                     geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                     
                                                     geom_text(
                                                       aes(
                                                         label=!!aes_percentage(relative_to='intersection'),
                                                         color=ifelse(!is.na(x$SNrank), 'show', 'hide')
                                                       ),
                                                       stat='count',
                                                       position=position_fill(vjust = .5),
                                                       
                                                     ),
                                                     
                                                     scale_y_continuous(labels=scales::percent_format()),
                                                     show_hide_scale,
                                                     rating_scale
                                                     
                                                   )
                                                 )
                                               )
                                               
      )+ggtitle(str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
      check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                    set_sizes=FALSE,
                    guides='collect',
                    n_intersections=N,
                    
                    stripes='white',
                    base_annotations=list(
                      '# of PSMs'=intersection_size(
                        counts=TRUE,
                      )
                    ),
                    annotations =list(
                      "Ranked Intensity  "=list(
                        aes=aes(x=intersection, fill=df_1[[1]]$SNrank),
                        
                        geom=list(
                          geom_bar(stat='count', position='fill', na.rm=TRUE),
                          themes=theme(legend.position="bottom", legend.box = "horizontal"),
                          geom_text(
                            aes(
                              label=!!aes_percentage(relative_to='intersection'),
                              color=ifelse(!is.na(df_1[[1]]$SNrank), 'show', 'hide')
                            ),
                            stat='count',
                            position=position_fill(vjust = .5),
                            
                          ),
                          
                          scale_y_continuous(labels=scales::percent_format()),
                          show_hide_scale,
                          rating_scale
                          
                        )
                      )
                    )
                    
      )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
      y<-get_legend(check1$patches$plots[[1]])
      data<-unlist(lapply(check,function(x) x$labels$title))
      check<-check[order(data)]
      P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
      
      return(print(P))
      
    }else{
      df_1<-dplyr::bind_rows(df_)
      df_1<-df_1%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature")
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$dataset<-as.factor(df_1$dataset)
      df_1$sample_name<-as.factor(df_1$sample_name)
      df_1$rank<-ifelse(df_1$rank==3,"high",ifelse(df_1$rank==2,"medium","low"))
      df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
      
      df_1<-df_1%>% 
        dplyr::group_split(uniqueID,sample_id)
      df_1<-lapply(df_1,function(x) x[1,]) 
      df_1<-dplyr::bind_rows(df_1)
      
      df_1$missing_pct<-as.numeric(df_1$missing_pct)
      df_1$missing_pct<-round(df_1$missing_pct,0)
      
      
      df_1<-df_1 %>% dplyr::group_split(sample_name)
      df_1<-purrr::map(df_1,function(x)
        pivot_wider(x,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA) %>%
          dplyr::select(-C,-I,-CC,-missing,-dataset,-sample_id,-temp_ref,-id))
      
      df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(rank,sample_name))
      df_1<-purrr::map(df_1,function(x) 1+x[,lapply(x,class)=="numeric"])
      df<-lapply(df_1,function(x) colnames(x))
      df<-lapply(df,function(x) str_replace(x,x,
                                            paste0("missing ",x,"%", sep = " ")))
      
      
      
      
      df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
      df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
      df_1<-purrr::map(df_1,function(x)x[!is.na(x$rank),])
    }
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    check<-list()
    check<- purrr::map(df_1,function(x)upset(x,names(x),
                                             min_degree=1,
                                             set_sizes=FALSE,
                                             guides='collect',
                                             n_intersections=N,
                                             height_ratio = 0.7,
                                             stripes='white',
                                             base_annotations=list(
                                               '# of Replicates'=intersection_size(
                                                 counts=TRUE,
                                               )
                                             ),
                                             annotations =list(
                                               "Bin %"=list(
                                                 aes=aes(x=intersection, fill=x$rank),
                                                 
                                                 geom=list(
                                                   geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                   
                                                   geom_text(
                                                     aes(
                                                       label=!!aes_percentage(relative_to='intersection'),
                                                       color=ifelse(!is.na(rank), 'show', 'hide')
                                                     ),
                                                     stat='count',
                                                     position=position_fill(vjust = .5),
                                                     
                                                   ),
                                                   
                                                   scale_y_continuous(labels=scales::percent_format()),
                                                   show_hide_scale,
                                                   rating_scale
                                                   
                                                 )
                                               )
                                             )
                                             
    )+ggtitle(str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
    check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                  set_sizes=FALSE,
                  guides='collect',
                  n_intersections=N,
                  
                  stripes='white',
                  base_annotations=list(
                    '# of Replicates'=intersection_size(
                      counts=TRUE,
                    )
                  ),
                  annotations =list(
                    "Ranked Intensity  "=list(
                      aes=aes(x=intersection, fill=df_1[[1]]$rank),
                      
                      geom=list(
                        geom_bar(stat='count', position='fill', na.rm=TRUE),
                        themes=theme(legend.position="bottom", legend.box = "horizontal"),
                        geom_text(
                          aes(
                            label=!!aes_percentage(relative_to='intersection'),
                            color=ifelse(!is.na(rank), 'show', 'hide')
                          ),
                          stat='count',
                          position=position_fill(vjust = .5),
                          
                        ),
                        
                        scale_y_continuous(labels=scales::percent_format()),
                        show_hide_scale,
                        rating_scale
                        
                      )
                    )
                  )
                  
    )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
    y<-get_legend(check1$patches$plots[[1]])
    
    P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
    print(P)
  }else{
    ########Plot separately
    
    df_<-df_ %>% subset(sample_name== condition) 
    
    df_<-df_ %>% dplyr::rename("uniqueID"="Accession","C"="temperature","I"="value")
    
    df_$rank<-ifelse(df_$rank==1,"high",ifelse(df_$rank==2,"medium","low"))
    df_$rank<-factor(df_$rank,levels=c("high","medium","low"))
    df_<-df_ %>% dplyr::filter(sample_name==condition)%>%
      dplyr::group_split(uniqueID) 
    df_<-lapply(df_,function(x) x[1,]) %>% dplyr::bind_rows(.)
    df_$missing_pct<-as.numeric(df_$missing_pct)
    df_$missing_pct<-round(df_$missing_pct,0)
    df_$dataset<-as.factor(df_$dataset)
    
    d2<-dplyr::bind_rows(df_)%>% dplyr::filter(dataset=="vehicle") 
    d2$uniqueID<-as.factor(d2$uniqueID)
    d2$sample_name<-as.factor(d2$sample_name)
    d2$C<-as.factor(d2$C)
    d2$rank<-as.factor(d2$rank)
    
    
    d1<-dplyr::bind_rows(df_)%>% dplyr::filter(dataset=="treated") 
    d1$uniqueID<-as.factor(d1$uniqueID)
    d1$sample_name<-as.factor(d1$sample_name)
    d1$C<-as.factor(d1$C)
    d1$rank<-as.factor(d1$rank)
    
    
    d3<-rbind(d1,d2)
    
    d3<-tidyr::pivot_wider(d3,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d1<-pivot_wider(d1,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d2<-pivot_wider(d2,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    
    d1<-d1%>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    d2<-d2 %>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    d3<-d3 %>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    
    rd1<-d1$rank
    rd2<-d2$rank
    rd3<-d3$rank
    
    d1<-1+d1[,lapply(d1,class)=="numeric"]
    d2<-1+d2[,lapply(d2,class)=="numeric"]
    d3<-1+d3[,lapply(d3,class)=="numeric"]
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    colnames(d1) <- paste("missing ", colnames(d1),"%", sep = " ")
    colnames(d2) <- paste("missing ", colnames(d2),"%", sep = " ")
    colnames(d3) <- paste("missing ", colnames(d3),"%", sep = " ")
    
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    p<-upset(d2,names(d2),
             min_degree=1,
             set_sizes=FALSE,
             n_intersections=N,
             stripes='white',
             base_annotations=list(
               '# of Replicates'=intersection_size(
                 counts=TRUE,
               )
             ),
             annotations =list(
               "Ranked Intensity %"=list(
                 aes=aes(x=intersection, fill=d2$rank),
                 geom=list(
                   geom_bar(stat='count', position='fill', na.rm=TRUE),
                   geom_text(
                     aes(
                       label=!!aes_percentage(relative_to='intersection'),
                       color=ifelse(!is.na(rank), 'show', 'hide')
                     ),
                     stat='count',
                     position=position_fill(vjust = .5)
                   ),
                   scale_y_continuous(labels=scales::percent_format()),
                   
                   show_hide_scale,
                   rating_scale
                 )
               )
             )
             
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))
    
    p2<-upset(d1,names(d1),name="Sample bins",
              min_degree=1,
              n_intersections=N,
              set_sizes=FALSE,
              stripes='white',
              base_annotations=list(
                '# of Replicates'=intersection_size(
                  counts=TRUE,
                )
              ),
              annotations =list(
                "Ranked Intensity %"=list(
                  aes=aes(x=intersection, fill=d1$rank),
                  geom=list(
                    geom_bar(stat='count', position='fill', na.rm=TRUE),
                    
                    geom_text(
                      aes(
                        label=!!aes_percentage(relative_to='intersection'),
                        color=ifelse(!is.na(rank), 'show', 'hide')
                      ),
                      stat='count',
                      position=position_fill(vjust = .5)
                    ),
                    scale_y_continuous(labels=scales::percent_format()),
                    
                    show_hide_scale,
                    rating_scale
                  )
                )
              )
              
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    
    d<-names(d3)[which(!names(d3)=="rank")]
    dr<-which(names(d3)=="rank")
    d3<-d3 %>% dplyr::select(-rank)
    d<-as.numeric(unlist(str_extract_all(d,"[[:digit:]]+")))
    ints<-names(d3)[c(order(d,decreasing=FALSE))] %>% as.list(.)
    d3$rank<-rd3
    
    p3<-ComplexUpset::upset(d3,names(d3), 
                            n_intersections=N,
                            set_sizes=FALSE,
                            stripes='white',
                            base_annotations=list(
                              '# of Replicates'=intersection_size(
                                counts=TRUE,
                              )
                            ),
                            annotations =list(
                              "Ranked Intensity %"=list(
                                aes=aes(x=intersection, fill=d3$rank),
                                geom=list(
                                  geom_bar(stat='count', position='fill', na.rm=TRUE),
                                  
                                  geom_text(
                                    aes(
                                      label=!!aes_percentage(relative_to='intersection'),
                                      color=ifelse(!is.na(rank), 'show', 'hide')
                                    ),
                                    stat='count',
                                    position=position_fill(vjust = .5)
                                  ),
                                  scale_y_continuous(labels=scales::percent_format()),
                                  
                                  show_hide_scale,
                                  rating_scale
                                )
                              )
                            ) 
                            
    )+ggtitle(str_replace(d2$condition[1],"S",paste0("\u03A6")))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    return(list(plot(p),plot(p2),plot(p3)))
  }
}


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
  df_1<-purrr::map(df_1,function(x){x %>% dplyr::select(-LineRegion)})#remove Line Region column from one dataset before merging
  
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
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                         dplyr::group_nest(LineRegion,uniqueID) %>%
                         dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=purrr::map(M1,function(x){deviance(x)}),
                                       Rsq=purrr::map(M1,function(x){summary(x)$r.squared}), 
                                       dataset="vehicle",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1):(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1)])))
    
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}), 
                                         dataset="treated",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
    
    
    
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1):(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1)])))
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                         dplyr::group_nest(LineRegion,uniqueID) %>% 
                         dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx( x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="null",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1:which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1])))
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
      testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-dataset)
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
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                         dplyr::group_nest(LineRegion,uniqueID) %>%
                         dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="vehicle",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1):(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1)])))
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}),
                                         dataset="treated",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1):(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1)])))
    
    
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                         dplyr::group_nest(LineRegion,uniqueID) %>% 
                         dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="null",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values[(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))-1):(which(abs(x$M1[[1]]$fitted.values-0.5)==min(abs(x$M1[[1]]$fitted.values-0.5)))+1)])))
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
      testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-dataset)
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
    #tlresults<-tlresults %>% keep(function(x)  sum(data.frame(x)[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <10)#move data with extremely large RSS values 
    # tlresults<-tlresults %>% keep(function(x) sum(data.frame(x)[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <1.5)
    tlresults<-tlresults %>% keep(function(x) sum(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),na.rm=TRUE) > sum(unlist(x[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),na.rm=TRUE))#remove data with extremely large RSS values 
    tlresults<-tlresults %>% keep(function(x) mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'Tm']),na.rm=TRUE) < mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'Tm']),na.rm=TRUE))
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
  Nsum<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "null")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
  
  #get the summed rss values for vehicle
  Rssv<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>%head(.,1))
  #get the summed rss values for treated
  Rsst<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
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
  names(Nsum)<-c("RSSn","Tmn","dataset","uniqueID")
  Nsum<-Nsum %>% dplyr::filter(uniqueID %in% CID)
  Nsum<-Nsum %>% dplyr::mutate(id=rownames(Nsum))
  
  Nsum$dataset<-as.factor(Nsum$dataset)
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
    Df1<-Dsum %>% dplyr::left_join(tlresults,by=c("uniqueID")) %>% as.data.frame(.) %>% dplyr::rename("dataset"="dataset.y") %>% 
      dplyr::select(-dataset.x)
    Df1<-Df1 %>% dplyr::group_split(uniqueID)
    
  }
  
  df1<-list()
  #get uniqueID and dataset for stable proteins with decreasing RSS differences
  df1<-purrr::map(Df1,function(x) x %>% dplyr::select(uniqueID,dataset) %>% head(.,1))
  df1<-data.frame(dplyr::bind_rows(df1))
  
  #unlist to data.frame
  #order the original data by RSS differences
  #
  DFN<- dplyr::bind_rows(DFN)
  DFN$uniqueID<-as.vector(DFN$uniqueID)
  df1$uniqueID<-as.vector(df1$uniqueID)
  
  
  df2<-df1 %>% dplyr::right_join(DFN,by=c("uniqueID")) %>% dplyr::rename("dataset"="dataset.y") %>% 
    dplyr::select(-dataset.x)
  
  
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
  
  null<-Df1 %>% subset(dataset == "null")
  
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
  
  Pred$Treatment<-null$dataset[1]##################
  Pred<-na.omit(Pred)
  Pred$C<-as.numeric(as.vector(Pred$C))
  Pred$I<-as.numeric(as.vector(Pred$I))
  PLN<-ggplot2::ggplot(Pred, ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1$uniqueID[1],str_replace(DF1$sample_name[1],"S",paste0("\u03A6"))))+
    ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ 
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
    annotate("text", x=60, y=min(Pred$I),label=paste("RSS= ",round(sum(unlist(null$rss)),3)))+
    annotate("text",
             x = Pred[which(round(Pred$fit,1)==0.5)[1],]$C,
             y = 0.45,
             label=paste(Pred[which(round(Pred$fit,1)==0.5)[1],]$C),
             colour="red"
    )+theme(legend.position="bottom")
  
  
  DF_f<-df2 %>%subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated')) %>% subset(dataset=="vehicle")
  
  vehicle<-Df1 %>% subset(dataset == "vehicle")
  
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
  
  Pred1$Treatment<-vehicle$dataset[1]##################
  Pred1<-na.omit(Pred1)
  rownames(Pred1)<-1:nrow(Pred1)
  Pred1$C<-as.numeric(as.vector(Pred1$C))
  Pred1$I<-as.numeric(as.vector(Pred1$I))
  
  
  DF_f1<-data.frame()
  DF_f1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated'))
  if(length(unique(DF_f1$dataset))==1){
    DF_f1<-DF_f1 
  }else{
    DF_f1<-DF_f1 %>% subset(dataset =="treated")
  }
  
  treated<-data.frame()
  treated<-Df1 %>% subset(dataset == "treated")
  
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
  
  Pred2$Treatment<-treated$dataset[1]##################
  Pred2<-na.omit(Pred2)
  rownames(Pred2)<-as.vector(1:nrow(Pred2))
  #Area under the curve using trapezoid rule
  
  P1_AUC <- pracma::trapz(Pred1$I[(which(abs(Pred1$I-0.5)==min(abs(Pred1$I-0.5)))-1):(which(abs(Pred1$I-0.5)==min(abs(Pred1$I-0.5)))+1)])
  P2_AUC <- pracma::trapz(Pred2$I[(which(abs(Pred2$I-0.5)==min(abs(Pred2$I-0.5)))-1):(which(abs(Pred2$I-0.5)==min(abs(Pred2$I-0.5)))+1)])
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
  DF1$dataset<-as.factor(DF1$dataset)
  
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
    Pred1<-lapply(Pred1,function(x) x %>% dplyr::mutate(replicate=row.names(.)))
    Pred2<-Pred2 %>% dplyr::group_split(uniqueID,C)
    Pred2<-lapply(Pred2,function(x) x %>% dplyr::mutate(replicate=row.names(.)))
    
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
               x = 50,
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
               x = 50,
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
      miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
      miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
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
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.35, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "vehicle"),'rss']))+
                                                                                                    sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "treated"),'rss'])),3)),size=3.5)+
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.45, label=  paste("\u0394", "AUC = ",AUCd),size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.55, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"),size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.65, label= paste("missing  ",Pred1$missing_v[1],"%"),colour="#00BFC4",size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.75, label= paste("missing  ",Pred2$missing_t[1],"%"),colour="#F8766D",size=3.5)+
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
               y = 0.55,
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
               y = 0.55,
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
      miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
      miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
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
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.35, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "vehicle"),'rss']))+
                                                                                                    sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "treated"),'rss'])),3)),size=3.5)+
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.45, label=  paste("\u0394", "AUC = ",AUCd),size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.55, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"),size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.65, label= paste("missing  ",Pred1$missing_v[1],"%"),colour="#00BFC4",size=3.5)+ 
        ggplot2::annotate("text", x=min(Pred2$C)+5, y= -0.75, label= paste("missing  ",Pred2$missing_t[1],"%"),colour="#F8766D",size=3.5)+
        annotate("text",
                 x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = 0.55,
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
      miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
      miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
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
        ggplot2::annotate("text", x=min(Pred$C), y=min(Pred2$I)+0.45, label= paste("missing % v",Pred2$missing_v[1]))+ 
        ggplot2::annotate("text", x=min(Pred$C), y=min(Pred2$I)+0.35, label= paste("missing % t",Pred2$missing_t[1]))+
        annotate("text",
                 x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1),
                 y = 0.55,
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
                            dataset=c(Pred$dataset[1],Pred2$dataset[1]),
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


#Spline
spstat<-function(DF,df,df1,norm=FALSE,Ftest=TRUE,show_results=TRUE,filters=TRUE){
  if(any(class(df)=="list")){
    df<-df %>% purrr::keep(function(x) is.data.frame(x))
    df1<-df1 %>% purrr::keep(function(x) is.data.frame(x))
    DF<-DF %>% purrr::keep(function(x) is.data.frame(x))
    
    df<-dplyr::bind_rows(df)
    df1<-dplyr::bind_rows(df1)
    DF<-dplyr::bind_rows(DF)
  }
  #plot spline results
  
  df1$I<-as.numeric(as.vector(df1$I))
  df$I<-as.numeric(as.vector(df$I))
  DF$I<-as.numeric(as.vector(DF$I))
  #mutate to get CV values
  DF<-DF %>% dplyr::group_split(C,uniqueID) 
  DF<- purrr::map(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df<-df %>% dplyr::group_split(C,uniqueID,dataset) 
  df<- purrr::map(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df1<-df1 %>% dplyr::group_split(C,uniqueID,dataset) 
  df1<- purrr::map(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  
  #convert to data frame
  
  df<-dplyr::bind_rows(df)
  df1<-dplyr::bind_rows(df1)
  DF<-dplyr::bind_rows(DF)
  
  #aggregate column 
  #switch from factor to numeric
  #convert factor to numeric columns
  df1$C<-as.numeric(as.vector(df1$C))
  df$C<-as.numeric(as.vector(df$C))
  DF$C<-as.numeric(as.vector(DF$C))
  
  df1$I<-as.numeric(as.vector(df1$I))
  df$I<-as.numeric(as.vector(df$I))
  DF$I<-as.numeric(as.vector(DF$I))
  
  #remove NA vaues
  df <- df[,which(unlist(lapply(df, function(x)!all(is.na(x))))),with=F]
  df1 <- df1[,which(unlist(lapply(df1, function(x)!all(is.na(x))))),with=F]
  DF <- DF[,which(unlist(lapply(DF, function(x)!all(is.na(x))))),with=F]
  #convert back to list
  DF<-DF %>% dplyr::group_split(uniqueID)
  df<-df %>% dplyr::group_split(uniqueID)
  df1<-df1 %>% dplyr::group_split(uniqueID)
  
  if(!isTRUE(norm)){
    
    #alternative spline fit method : Generalized Additive Models
    #fit penalized splines
    # if(any(names(df[[1]])=="time_point")){
    #   df<-purrr::map(df,function(x) dplyr::bind_rows(x) %>% dplyr::group_split(uniqueID,time_point))
    # }
    m <- purrr::map(df,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(x$I ~ s(x$C,k=5), data = x , method = "REML")))))
    m<-m %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m<-purrr::map(m,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                   sum = list(summary(.$M1[[1]])),
                                                   Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                                   rss=deviance(.$M1[[1]]),
                                                   CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
                                                   AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
                                                   rsq=summary(x$M1[[1]])$r.sq,
                                                   n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
                                                   sample_name=.$sample_name))
    m1 <- purrr::map(df1,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data = x , method = "REML")))))
    m1<-m1 %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m1<-purrr::map(m1,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                     sum = list(summary(.$M1[[1]])),
                                                     Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                                     rss=deviance(.$M1[[1]]),
                                                     CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
                                                     AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
                                                     rsq=summary(x$M1[[1]])$r.sq,
                                                     n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
                                                     sample_name=.$sample_name))
    #m1<-lapply(df1,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    
    
    mn<- purrr::map(DF,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data =x, method = "REML",fit=TRUE)))))
    mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    mn<-purrr::map(mn,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                     sum = list(summary(.$M1[[1]])),
                                                     Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                                     rss=deviance(.$M1[[1]]),
                                                     CV_pct=ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
                                                     AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
                                                     rsq=summary(.$M1[[1]])$r.sq,
                                                     n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
                                                     sample_name=.$sample_name,
                                                     RSS0=deviance(.$M1[[1]])
    ))
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
    m<-m %>% dplyr::group_split(uniqueID)
    m1<-m1%>% dplyr::group_split(uniqueID)
    mn<-mn%>% dplyr::group_split(uniqueID)
    
    #calculate RSS
    m<-purrr::map2(m,m1,function(x,y) x %>% dplyr::mutate(RSSv=deviance(x$M1[[1]]),
                                                          RSSt=deviance(y$M1[[1]]),
                                                          RSS1=deviance(x$M1[[1]])+deviance(y$M1[[1]])
    ))
    m1<-purrr::map2(m,m1,function(x,y) y %>% dplyr::mutate(RSSv=deviance(x$M1[[1]]),
                                                           RSSt=deviance(y$M1[[1]]),
                                                           RSS1=deviance(x$M1[[1]])+deviance(y$M1[[1]])
    ))
    
    check<-purrr::map2(m1,mn,function(x,y) x %>% dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
    check<-dplyr::bind_rows(check)
    
    #4031 proteins have rssDiff> 0 
    
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(m)
    mean1_1<-dplyr::bind_rows(m1)
    mean3<-dplyr::bind_rows(mn)
    #equal column names
    mean1<-mean1  %>% dplyr::mutate(rss=RSS1)%>% dplyr::select(-RSSv,-RSSt,-RSS1)
    mean1_1<-mean1_1  %>% dplyr::mutate(rss=RSS1)%>% dplyr::select(-RSSv,-RSSt,-RSS1)
    mean3<-mean3  %>% dplyr::mutate(rss=RSS0)%>% dplyr::select(-RSS0)
    
    #filter out proteins with neg RSSdiff
    if(isTRUE(filters)){
      check<-check %>% dplyr::filter(rssDiff>0,missing_pct<=20,rsq>0.8)
      check<-unique(check$uniqueID)
      mean1<-mean1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
      mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
      mean3<-mean3 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
      
      check<-intersect(mean1$uniqueID,mean1_1$uniqueID)
      
      mean1<-mean1 %>% dplyr::filter(uniqueID %in% check)
      mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check)
      mean3<-mean3 %>% dplyr::filter(uniqueID %in% check)
      
    }
    
    
    #split into lists by uniqueID
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
    
    #Cliff
    results<-dplyr::bind_rows(mean1,mean1_1,mean3)
    
    if(isTRUE(Ftest)){
      
      #Calculate rss0 and rss1 null vs alt
      rss0<-purrr::map(mean3,function(x)data.frame(RSS0=x$rss,
                                                   se0=summary.gam(x$M1[[1]])$se[1],#standard error
                                                   pN0=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                                                   rsq0=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
                                                   np0=summary.gam(x$M1[[1]])$np[1],#number of parameters
                                                   rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                                                   m=summary.gam(x$M1[[1]])$m[1]))
      #of smooth terms in the model
      
      #mean1<-purrr::map2(mean1,mean1_1,function(x,y)x %>% purrr::keep(x$uniqueID[1] %in% y$uniqueID[1]))
      #generate a grid of new values for C
      
      
      #make predictions to generate point confidence intervals
      
      
      rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS=x$rss,
                                                              se=summary.gam(x$M1[[1]])$se[1],#standard error
                                                              pN1=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                                                              rsq=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
                                                              
                                                              rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                                                              m=summary.gam(x$M1[[1]])$m[1],
                                                              Tm = y$Tm[[1]]-x$Tm[[1]],
                                                              se1=summary.gam(y$M1[[1]])$se[1],#standard error
                                                              pN2=summary.gam(y$M1[[1]])$edf[1],#effective degrees of freedom treated
                                                              rsq1=summary.gam(y$M1[[1]])$r.sq[1],#r-squared
                                                              
                                                              rdf1=summary.gam(y$M1[[1]])$residual.df[1],#residual degrees of freedom
                                                              m1=summary.gam(y$M1[[1]])$m[1],
                                                              pA=sum(summary.gam(x$M1[[1]])$edf[1],summary.gam(y$M1[[1]])$edf[1],na.rm=TRUE)))
      
      # rss1<-purrr::map2(rss1,mean1,function(x,y) x %>% dplyr::mutate(fit_v=list(ifelse(class(try(mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",
      #                                                                                                              se.fit=TRUE,newdata.guaranteed = TRUE)))=='try-error',NA,
      #                                                                                  mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",se.fit=TRUE,newdata.guaranteed = TRUE)))))
      
      #bind predicted values to main data frame
      mean3<-purrr::map2(mean3,rss0,function(x,y)cbind(x,y))
      mean1<-purrr::map2(mean1,rss1,function(x,y)cbind(x,y))
      #mean1_1<-purrr::map2(mean1_1,rss1,function(x,y)cbind(x,y))
      #params for null and alternative models
      f0<-lapply(mean3,function(x) data.frame(np=length(x$M1[[1]]$fitted.values)))#number of measurements
      f1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(nA=sum(nrow(x),nrow(y))))#number of measurements alternative
      #bind data frames
      mean3<-purrr::map2(mean3,f0,function(x,y) cbind(x,y))
      mean1<-purrr::map2(mean1,f1,function(x,y) cbind(x,y))
      mean1_1<-purrr::map2(mean1_1,f1,function(x,y) cbind(x,y))
      #calculate parameters 
      pN<-purrr::map(mean3,function(x) x %>% dplyr::select(pN0,np))
      pA<-purrr::map(mean1,function(x)x %>% dplyr::select(pA,nA))
      #degrees of freedom before
      d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN0))
      d2<-purrr::map(pA,function(x)data.frame(d2=x$nA-x$pA))
      
      #delta RSS
      rssDiff<-purrr::map2(rss0,rss1,function(x,y) data.frame(dRSS=x$RSS0-y$RSS))
      
      d2<-lapply(d2,function(x) x$d2[1])
      d1<-lapply(d1,function(x) x$d1[1])
      #RSS1 numerator
      Fvals<-purrr::map2(rssDiff,d1,function(x,y) data.frame(fNum=x$dRSS/y[1]))
      #Rss denominator
      Fd<-purrr::map2(rss1,d2,function(x,y) data.frame(fDen=x$RSS/y[1]))
      
      Fvals<-purrr::map2(Fvals,Fd,function(x,y) data.frame(fStatistic=x$fNum/y$fDen))
      Fvals<-purrr::map2(Fvals,d1,function(x,y) x %>% dplyr::mutate(df1=y[1]))
      Fvals<-purrr::map2(Fvals,d2,function(x,y) x %>% dplyr::mutate(df2=y[1]))
      Fvals<-purrr::map2(Fvals,rssDiff,function(x,y) cbind(x,y))
      #add p-vals
      Fvals<-purrr::map(Fvals,function(x) x %>% dplyr::mutate(pValue = 1 - pf(fStatistic, df1 = x$df1, df2 = x$df2),
                                                              pAdj = p.adjust(pValue,"BH")))
      Fvals<-purrr::map2(Fvals,mean1,function(x,y) x %>% dplyr::mutate(uniqueID=y$uniqueID[1]))
      
      mean1<-purrr::map(mean1,function(x) data.frame(x)%>% dplyr::select(-id,-temp_ref,-C,-I,-M1,-sum) )
      
      mean1<-purrr::map(mean1,function(x) x %>% distinct(.) )
      
      #convert results to list
      testResults<-purrr::map2(mean1,Fvals,function(x,y)x%>% dplyr::right_join(y,by=c("uniqueID")))
      testResults<-purrr::map(testResults,function(x) x[1,])
      testResults<-dplyr::bind_rows(testResults)
      mean1<-dplyr::bind_rows(mean1)
      Unscaled<-ggplot(testResults)+
        geom_density(aes(x=fStatistic),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=fStatistic,y= df(fStatistic,df1=df1,df2=df2)),color="darkred",size = 1.5) +
        theme_bw() +
        ggplot2::xlab("F-values")+
        ggplot2::xlim(0,0.05)
      print(Unscaled)
      #scale variables
      M<-median(testResults$dRSS,na.rm=TRUE)
      V<-mad(testResults$dRSS,na.rm=TRUE)
      #alternative scaling factor sig0-sq
      altScale<-0.5*V/M
      #filter out negative delta rss
      testResults<-testResults %>% dplyr::filter(dRSS>0)
      #effective degrees of freedom
      ed1<-MASS::fitdistr(x=testResults$dRSS, densfun = "chi-squared", start = list(df=2))[["estimate"]]
      ed2<-MASS::fitdistr(x=testResults$rss, densfun = "chi-squared", start = list(df=2))[["estimate"]]
      #scale data
      testScaled <-testResults %>%
        dplyr::mutate(rssDiff = .$dRSS/altScale,
                      rss1 =.$RSS/altScale,
                      d1=ed1,
                      d2=ed2)
      #
      #new F-test
      testScaled<-testScaled %>% dplyr::mutate(Fvals=(dRSS/rss1)*(d2/d1))
      Fvals<-testScaled$Fvals
      d1<-testScaled$d1
      d2<-testScaled$d2
      
      #scaled values
      TestScaled<-ggplot(testScaled)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() + 
        ggplot2::xlab("F-values")+
        ggplot2::xlim(0,0.05)
      print(TestScaled)
      #Define checked as filtered protein IDs
      check<-testScaled$uniqueID
      test<-testScaled %>% dplyr::filter(.$pValue<0.01)
      test$d1<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test$d2<-MASS::fitdistr(x=test$RSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      
      testS<-ggplot(test)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
      print(testS)
      
      
      mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
      mean3<-dplyr::bind_rows(mean3)
      mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      if(isTRUE(show_results)){
        return(testScaled)
      }else{
        return(list(results,testS,Unscaled))
      }
    }
    return(results)
  }
  return(results)
}
spf<-function(spresults,DFN,filters = TRUE){
  if(!any("missing" %in% names(spresults))){
    spresults<-spresults %>% dplyr::mutate(missing=NA,missing_pct=NA,rank=NA)
    DFN<-DFN%>% dplyr::mutate(missing=NA,missing_pct=NA,rank=NA)
  }
  spresults<-spresults %>% dplyr::group_split(uniqueID)
  
  if(!isTRUE(filters))
  {
    sl<-purrr::map(seq_len(length(spresults)),function(x) as.numeric({paste(x)})) 
    sp<-purrr::map2(spresults,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))  
    sp<-dplyr::bind_rows(sp)  
    df1<-data.frame(uniqueID = unique(sp$uniqueID))  
    df2<-dplyr::bind_rows(DFN)  
    df2$C<-as.numeric(as.vector(df2$C)) 
    df2$I<-as.numeric(as.vector(df2$I))  
    if(any(names(sp) %in% c("missing.x","LineRegion","N"))){
      names<-intersect(names(df2),names(sp))
      df2$id<-as.numeric(df2$id)
      df2<-sp %>% left_join(df2, by = names) 
      
    }else{
      names<-intersect(names(df2),names(sp))
      df2$id<-as.numeric(df2$id)
      df2<-sp %>% left_join(df2, by = names) 
    }
    Df1<-spresults
  }else{
    #Apply filters 
    #keep the positive AUC differences
    
    spresults<-spresults %>% keep(function(x) mean(x$AUC[x$dataset=="treated"],na.rm=TRUE)>mean(x$AUC[!x$dataset=="vehicle"],na.rm=TRUE))
    
    spresults<-spresults %>% keep(function(x) max(x$lambda)<1)
    
    if (is.null(nrow(spresults))){
      return(warning("all proteins filtered out by AUC and lambda value"))
    }
    #get Tm and RSS differences
    sp<-purrr::map(spresults, function(x) x %>% dplyr::mutate(Tmd= x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'Tm'][[1]] - x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'Tm'][[1]],
                                                              RSSd = sum(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']) - sum(x[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),
                                                              AUCd = x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'AUC'])[[1]]- x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'AUC'][[1]])
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
      names<-intersect(names(df2),names(sp))
      df2<-sp %>% left_join(df2, by = names) 
    }else{
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
spCI<-function(i,df1,df2,Df1,df.temps,overlay=TRUE,alpha,residuals=FALSE,simulations=FALSE,CI=TRUE,Peptide=FALSE,CARRIER=TRUE){
  if(isTRUE(CARRIER)){
    df.temps<-length(unique(df.temps$temperature))-1
  }else{
    df.temps<-length(unique(df.temps$temperature))
  }
  null<-data.frame()
  i<-i
  
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-as.character(df2$uniqueID)
  
  
  #get original data
  ###########################################
  df1<-df1$uniqueID[i]
  DF1<-df2[which(df2$uniqueID %in% df1),]
  Df1<-data.frame(Df1[[i]])
  null<-Df1[which(Df1$uniqueID %in% df1 & Df1$dataset %in% "null"),]
  ###########################################
  DF_f<-df2 %>% subset(uniqueID %in% df1 & dataset %in% "vehicle")
  vehicle<-Df1 %>% subset(uniqueID %in% df1 & dataset %in% "vehicle")
  ###########################################
  DF_f1<-df2%>% subset(uniqueID %in% df1 & dataset %in% "treated")
  treated<-Df1 %>% subset(uniqueID == df1 & dataset == "treated")
  
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
  BSvarN<-df2 %>% subset(uniqueID == df1 ) 
  BSvar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")
  BSvar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")
  BSVarN<-df2 %>% subset(uniqueID == df1 ) %>%dplyr::group_by(C)# %>%  dplyr::mutate(I=mean(I))
  BSVar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")%>%dplyr::group_by(C)# %>% dplyr::mutate(I=mean(I))
  BSVar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")%>%dplyr::group_by(C)# %>% dplyr::mutate(I=mean(I))
  
  
  BSVarN<-BSVarN %>% dplyr::mutate(dataset="null") 
  BSVar<-BSVar %>% dplyr::mutate(dataset="vehicle")
  BSVar1<-BSVar1 %>% dplyr::mutate(dataset="treated")
  BSVarN$dataset<-as.factor(BSVarN$dataset)
  BSVar$dataset<-as.factor(BSVar$dataset)
  BSVar1$dataset<-as.factor(BSVar1$dataset)
  
  BSVar<-BSVar[!is.na(BSVar$I),]
  BSVar1<-BSVar1[!is.na(BSVar1$I),]
  BSVarN<-BSVarN[!is.na(BSVarN$I),]
  
  fit <-  stats::smooth.spline(x = BSVar$C, y=BSVar$I,cv=F)
  fit1<-  stats::smooth.spline(x = BSVar1$C, y=BSVar1$I,cv=F)
  fitN<-  stats::smooth.spline(x = BSVarN$C, y=BSVarN$I,cv=F)
  
  #####try GAM
  #fit penalized splines
  m <- mgcv::gam(I ~ s(C,k=5), data = BSVar , method = "ML")
  m1<-  mgcv::gam(I ~ s(C,k=5), data =BSVar1, method = "ML")
  mn<-  mgcv::gam(I ~ s(C,k=5), data = BSVarN, method = "ML")
  
  #####try GAM
  
  #Plot boostrapped  residuals with 95%CI
  #PLP<-plot(m, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
  #generate random values from a multivariate normal distribution
  
  #get some parmeters
  Vb <- vcov(m)
  newd <- with(BSVar, data.frame(C = seq(min(C), max(C), length = 30)))%>% as.data.frame(.)
  pred <- predict(m, newd, se.fit = TRUE)%>% as.data.frame(.)#get confidence intervals
  se.fit <- pred$se.fit
  #get some parmeters
  Vb1<- vcov(m1) 
  newd1<- with(BSVar1,data.frame(C = seq(min(C), max(C), length = 30)))%>% as.data.frame(.)
  pred1<- predict(m1,newd1,se.fit = TRUE) %>% as.data.frame(.)
  se.fit1<- pred1$se.fit
  #generate std
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
  pred$dataset<-"vehicle"
  pred$dataset<-as.factor(pred$dataset)
  pred$CI<-"vehicle"
  pred$CI<-as.factor(pred$CI)
  
  plot<-ggplot(pred,mapping= ggplot2::aes(x = C,y=fit ))+
    geom_point(BSVar, mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=factor(CC)))+
    geom_ribbon(aes(ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle("")+ylim(-0.4,1.6)+xlim(37,68)+theme(legend.position="bottom")
  
  pred1<- transform(cbind(data.frame(pred1),newd1),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit))
  pred1$dataset<-"treated"
  pred1$dataset<-as.factor(pred1$dataset)
  pred1$CI<-"treated"
  pred1$CI<-as.factor(pred1$CI)
  pred1$AUC<-pracma::trapz(pred1$fit[(which(abs(pred1$fit-0.5)==min(abs(pred1$fit-0.5)))-1):(which(abs(pred1$fit-0.5)==min(abs(pred1$fit-0.5)))+1)])-pracma::trapz(pred$fit[(which(abs(pred$fit-0.5)==min(abs(pred$fit-0.5)))-1):(which(abs(pred$fit-0.5)==min(abs(pred$fit-0.5)))+1)])
  pred1$AUC<-abs(round(pred1$AUC[1],3))
  
  pred1$RSS<- deviance(m1)+deviance(m)
  
  pred1$RSS<- round(pred1$RSS,3)
  #Residuals
  
  pred1$Tm<-round(treated$Tm[1]-vehicle$Tm[1],1)
  #missing values
  miss_v<-data.frame(NA)
  miss_t<-data.frame(NA)
  #max replicates
  roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  miss_v<-DF1%>% dplyr::filter(dataset=="vehicle") %>% unique(.)
  miss_t<-DF1%>% dplyr::filter(dataset=="treated") %>% unique(.)
  
  Pred<-data.frame(m$fitted.values,m$residuals)
  names(Pred)<- c("fit","rn")
  Pred$dataset<-as.factor("vehicle")
  BSVar$dataset<-as.factor("vehicle")
  Pred1<-data.frame(m1$fitted.values,m1$residuals)
  names(Pred1)<-c("fit","rn")
  Pred1$dataset<-as.factor("treated")
  Preds<-rbind(Pred,Pred1)
  BSVar1$dataset<-as.factor("treated")
  #get fitted value data
  fitted.values<-data.frame(C=BSVar$M1[[1]]$model$`x$C`,fit=predict.gam(BSVar$M1[[1]],se.fit=TRUE))
  names(fitted.values)<-c("C","fit","se.fit")
  fitted.values1<-data.frame(C=BSVar1$M1[[1]]$model$C,fit=predict.gam(BSVar1$M1[[1]],se.fit=TRUE))
  names(fitted.values1)<-c("C","fit","se.fit")
  
  num<-roundUpNice(length(unique(miss_v$C))*length(unique(miss_v$rss)))
  
  #append Tm values on predicted data
  pred1$Tm<-round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)
  if(isTRUE(residuals)){
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=dataset)) +ggplot2::geom_point()+ 
      ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1]," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    print(PLrs)
  }
  if(isTRUE(CI)){
    BSVar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")%>%dplyr::group_by(C) # %>% dplyr::mutate(I=mean(I))
    BSVar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")%>%dplyr::group_by(C)# %>% dplyr::mutate(I=mean(I))
    
    BSVar<-BSVar[!is.na(BSVar$I),]
    BSVar1<-BSVar1[!is.na(BSVar1$I),]
    
    BSVar<-BSVar %>% distinct(.)
    BSVar1<-BSVar1 %>% distinct(.)
    
    
    #fit penalized splines
    m <- mgcv::gam(I ~ s(C,k=5), data = BSVar , method = "ML")
    m1<-  mgcv::gam(I ~ s(C,k=5), data =BSVar1, method = "ML")
    mn<-  mgcv::gam(I ~ s(C,k=5), data = BSVarN, method = "ML")
    
    #####try GAM
    #get some parmeters
    Vb <- vcov(m)
    newd <- with(BSVar, data.frame(C = seq(min(C), max(C), length = 10)))%>% as.data.frame(.)
    BSVar <- BSVar %>% dplyr::mutate(fit=list(predict(m, newd, se.fit = TRUE)))
    
    #get some parmeters
    Vb1<- vcov(m1) 
    newd1<- with(BSVar1,data.frame(C = seq(min(C), max(C), length = 10)))%>% as.data.frame(.)
    BSVar1 <- BSVar1 %>% dplyr::mutate(fit=list(predict(m1, newd1, se.fit = TRUE)))
    if (any(names(BSVar)=="sample.x")){
      BSVar<-BSVar %>% dplyr::rename("sample"="sample.x")
      
    }
    if (any(names(BSVar1)=="sample.x")){
      BSVar1<-BSVar1 %>% dplyr::rename("sample"="sample.x")
      
    }
    
    #append missing value data
    if(length(unique(miss_v$C))==df.temps & length(miss_v$C)>df.temps){
      getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      
      BSVar$missing_v<-rep(getmode(miss_v$missing_pct),nrow(BSVar))
      BSVar$missing_t<-rep(getmode(miss_t$missing_pct),nrow(BSVar))
    }else{
      BSVar$missing_v<-rep((100*(num-df.temps)/num),nrow(BSVar))
      BSVar$missing_t<-rep((100*(num-df.temps)/num),nrow(BSVar))
    }
    p<-data.frame(BSVar$fit[[1]])
    p1<-data.frame(BSVar1$fit[[1]])
    
    fit_v<-p %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                               uprP=fit+(1.96*se.fit),
                               C= seq(min(BSVar$C), max(BSVar$C), length = 10),
                               dataset=BSVar$dataset[1],
                               CI=BSVar$dataset[1])
    fit_t<-p1 %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                                uprP=fit+(1.96*se.fit),
                                C= seq(min(BSVar1$C), max(BSVar1$C), length = 10),
                                dataset=BSVar1$dataset[1],
                                CI=BSVar1$dataset[1])
    
    id<-data.frame(sample=as.factor(unique(BSVar$sample)),replicate=as.factor(seq(unique(BSVar$sample))))
    id1<-data.frame(sample=as.factor(unique(BSVar1$sample)),replicate=as.factor(seq(unique(BSVar1$sample))))
    
    BSVar<-BSVar %>%  dplyr::right_join(id,by="sample")
    BSVar1<-BSVar1%>% dplyr::right_join(id1,by="sample")
    if(isTRUE(Peptide)){
      BSVar$PeptideGroup<-BSVar$replicate
      BSVar1$PeptideGroup<-BSVar1$replicate
      if(any(names(BSVar)=="rank_l")){
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
      }else{
        BSVar$Stroke<-0
        BSVar1$Stroke<-0
      }
      
      plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
        ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=PeptideGroup))+
        geom_point(data=BSVar[BSVar$Stroke==1,],
                   pch=21, fill=NA, size=4, colour="black", stroke=1)+
        ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::annotate("text", x=43, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
        annotate("text",
                 x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
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
        ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=PeptideGroup))+
        geom_point(data=BSVar1[BSVar1$Stroke==1,],
                   pch=21, fill=NA, size=4, colour="black", stroke=1)+
        ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=dataset), alpha = 0.2 ) +
        ggplot2::labs(y = "Relative Solubility",
                      x = "Temperature (\u00B0C)")+
        annotate("text",
                 x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="red",
                 size=3.5
        )+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
        ylim(-0.60,1.6)+xlim(37,68)+
        theme(legend.position="bottom")
      
      return(plot)
    }else{
      
      plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
        ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=replicate))+
        ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::annotate("text", x=43, y=-0.03, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.13, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.23, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+ 
        ggplot2::annotate("text", x=43, y=-0.33, label= paste("missing",round(BSVar$missing_v[1],0),"%"),size=3.5,colour="#00BFC4")+ 
        ggplot2::annotate("text", x=43, y=-0.43, label= paste("missing",round(BSVar$missing_t[1],0),"%"),size=3.5,colour="#F8766D")+
        annotate("text",
                 x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="blue",
                 size=3.5
        )+
        annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = 0.5, yend = 0.5,
                 colour = "blue",linetype=2)+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "blue",linetype=2)+theme(legend.position="bottom")
      
      
      plot<-plot1+
        ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=replicate))+
        ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=dataset), alpha = 0.2 ) +
        ggplot2::labs(y = "Relative Solubility",
                      x = "Temperature (\u00B0C)")+
        annotate("text",
                 x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="red",
                 size=3.5
        )+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
        ylim(-0.5,1.6)+xlim(37,68)+
        theme(legend.position="bottom")
      return(plot)
    }
  }else{
    if(isTRUE(Peptide)){
      BSVar$PeptideGroup<-BSVar$replicate
      BSVar1$PeptideGroup<-BSVar1$replicate
      if(any(names(BSVar)=="rank_l")){
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
      }else{
        BSVar$Stroke<-0
        BSVar1$Stroke<-0
      }
      
      plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
        ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=PeptideGroup))+
        geom_point(data=BSVar[BSVar$Stroke==1,],
                   pch=21, fill=NA, size=4, colour="black", stroke=1)+
        ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::annotate("text", x=43, y=-0.35, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
        annotate("text",
                 x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
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
        ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=PeptideGroup))+
        geom_point(data=BSVar1[BSVar1$Stroke==1,],
                   pch=21, fill=NA, size=4, colour="black", stroke=1)+
        ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::labs(y = "Relative Solubility",
                      x = "Temperature (\u00B0C)")+
        coord_cartesian(xlim = c(37,67))+annotate("text",
                                                  x = 2+round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                                                  y = -0.15,
                                                  label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                                                  colour="red",
                                                  size=3.5
        )+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
        ylim(-0.60,1.5)+xlim(37,68)+
        theme(legend.position="bottom")
      return(plot)
    }else{
      
      plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
        ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,shape=replicate))+
        ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
        ggplot2::annotate("text", x=43, y=-0.03, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.13, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
        ggplot2::annotate("text", x=43, y=-0.23, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+ 
        ggplot2::annotate("text", x=43, y=-0.33, label= paste("missing",round(BSVar$missing_v[1],0),"%"),size=3.5,colour="#00BFC4")+ 
        ggplot2::annotate("text", x=43, y=-0.43, label= paste("missing",round(BSVar$missing_t[1],0),"%"),size=3.5,colour="#F8766D")+
        annotate("text",
                 x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="blue",
                 size=3.5
        )+
        annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 y = 0.5, yend = 0.5,
                 colour = "blue",linetype=2)+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                 xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "blue",linetype=2)+theme(legend.position="bottom")
      
      
      plot<-plot1+
        ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=replicate))+
        ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
        ggplot2::labs(y = "Relative Solubility",
                      x = "Temperature (\u00B0C)")+
        coord_cartesian(xlim = c(37,67))+
        annotate("text",
                 x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                 y = -0.15,
                 label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                 colour="red",
                 size=3.5
        )+
        annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                 colour = "red",linetype=2)+
        annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                 colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
        ylim(-0.60,1.5)+xlim(37,68)+
        theme(legend.position="bottom")
      
      return(plot) 
    }
    
  }
  
}



##################################################
#Sigmoidal function with confidence intervals
###################################################


fitSingleSigmoid <- function(x , y, start =c(Pl=0, a = 550, b = 10))
{
  try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
          start = start,
          data = list(x=x,y=y),
          na.action = na.exclude,
          algorithm = "port",
          lower = c(0.0,1e-5,1e-5),
          upper = c(1.5,15000,250),
          control = nls.control(maxiter = 50)),
      silent = TRUE)
}

repeatFits <- function(x,y, start= c(Pl = 0, a = 550, b=10),
                       seed = NULL, alwaysPermute = FALSE, maxAttempts = 100){
  i <-0
  doFit <- TRUE
  doVaryPars <-alwaysPermute
  
  if(!is.null(seed)){
    set.seed(seed)
  }  
  while (doFit){
    startTmp <-start * (1+doVaryPars*runif(1, -0.5,0.5))
    m <- fitSingleSigmoid(x=x, y=y, start = startTmp)
    i <- i +1
    doFit <- inherits(m,"try-error") & i < maxAttempts
    doVaryPars <- TRUE
  }
  return(m)
}

computeRSS <- function(x,y,start = c(Pl = 0, a = 550, b=10),seed = NULL,
                       alwaysPermute = FALSE,maxAttempts = 50){
  #start model fitting
  fit <- repeatFits(x=x, y=y,start=start,seed=seed,alwaysPermute=alwaysPermute,maxAttempts = maxAttempts)
  if(!inherits(fit,"try-error")){
    #if model fit converged, calculate RSS and parameters
    resid <-residuals(fit)
    rss <- sum(resid^2,na.rm = TRUE)
    fittedValues <- sum(!is.na(resid))
    params <- coefficients(fit)
    #Compute Tm for comparison (not needed for HT)
    tm <- params[["a"]]/(params[["b"]] - log(1-params[["Pl"]])/(1/2 - params[["Pl"]]-1))
  } else {
    #if model did not converge, return default vals
    rss <- NA
    fittedValues <-0
    params <- c(Pl=NA,a=NA,b=NA)
    tm <- NA
  }
  out <- tibble(rss=rss, fittedValues = fittedValues, tm= tm,
                a = params[["a"]],b = params[["b"]], Pl = params[["Pl"]])
  return(out)
  
  
}
# #GoF by RSS comparison for each protein
# #RSS difference
# 
computeRSSdiff <- function(x,y,treatment,maxAttempts = 50, repeatsIfNeg = 100){
  rssDiff <- -1
  repeats <- 0
  alwaysPermute <- FALSE
  
  start0 = start1 <- c(Pl=0,a=550,b=10)
  
  while((is.na(rssDiff) | rssDiff<0) & repeats <= repeatsIfNeg){
    
    nullResults <- computeRSS(x=x, y=y, start = start0, seed=repeats,
                              maxAttempts = maxAttempts,
                              alwaysPermute = alwaysPermute)
    
    altResults <- tibble(x,y,treatment) %>%
      group_by(treatment) %>%
      purrr::map_dfr({
        fit = computeRSS(x=.$x, y = .$y,start = start1, seed=repeats,
                         maxAttempts = maxAttempts,
                         alwaysPermute = alwaysPermute)
        
      }) 
    rss0 <- nullResults$rss
    rss1 <-sum(altResults$rss)#combined alternative RSS values
    rssDiff <- rss0-rss1#difference between null and combined alternative RSS values
    
    if(is.na(rssDiff) | rssDiff <0){
      repeats <- repeats +1
      alwaysPermute <- TRUE
      start1 <- c(Pl = nullResults[["Pl"]],a = nullResults[["a"]],b=nullResults[["b"]])
      
    }
  }
  
  n0 <-nullResults$fittedValues
  n1 <-sum(altResults$fittedValues)
  
  tm <- altResults %>%
    mutate(key = paste0("tm_",treatment)) %>%
    dplyr::select(key,tm) %>% spread(key,tm)
  
  out <- tibble(rss0,rss1,rssDiff,n0,n1,repeats) %>%
    cbind(tm)
  return(out)
}
sigCI <- function(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...)
{
  method <- match.arg(method)
  
  format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
          "%")
  
  ## Taken from confint.nls
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- seq_along(pnames)
  if (is.numeric(parm)) 
    parm <- pnames[parm]
  
  ## Taken from confint.default and modified slightly to use t-distribution
  asCI <- function(object, parm, level)
  {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    #        pct <- stats:::format.perc(a, 3)
    pct <- format.perc(a, 3)
    fac <- qt(a, df.residual(object))
    
    parmInd <- match(parm, pnames)
    ci <- array(NA, dim = c(length(parmInd), 2), dimnames = list(parm, pct))
    ses <- sqrt(diag(vcov(object)))[parmInd]
    ci[] <- cf[parmInd] + ses %o% fac
    ci
  }
  
  ## Taken from confint.nls
  asProf <- function(object, parm, level)
  {
    
    utils::flush.console()
    object <- profile(object, which = parm, alphamax = (1 - level)/4)
    confint(object, parm = parm, level = level, ...)    
  }
  
  switch(method, asymptotic = asCI(object, parm, level), profile = asProf(object, parm, level))
}

sigC<-function(df_,Protein){
  Protein<-as.character(Protein)
  df_$C<-as.numeric(as.character(df_$C))
  DFN<-df_%>% dplyr::filter(uniqueID %in% Protein)
  df_1<-df_%>% dplyr::filter(uniqueID%in%Protein,dataset=="treated")
  df_<-df_%>% dplyr::filter(uniqueID%in%Protein,dataset=="vehicle")
  
  
  
  df_<-df_ %>% dplyr::rename("I"="I3")
  df_1<-df_1 %>% dplyr::rename("I"="I3")
  DFN<-DFN %>% dplyr::rename("I"="I3")
  
  nlm2<-df_  %>%   
    dplyr::mutate(fit=list(try(fitSingleSigmoid(C,I))))
  nlm2$Tm<-NA
  nlm2$Tm<-Tm(nlm2)
  if(class(nlm2$fit[[1]])=='try-error'){
    warning("the function could not fit the data")
  }
  #remove proteins with a plateau value not =  zero 
  CT<-nlm2 %>%
    dplyr::filter( !coef(fit[[1]])[[1]]==0)#find values where Pl = 0 
  if(nrow(CT)==0){
    
    dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95))
    nlm2<- nlm2  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
    
  }else{
    #get confidence intervals for sigmoidal function
    dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95)) #predicted values
    #obtain sigmoidal equation parameters
    nlm2<- CT  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
    
  }
  
  
  #ready confidence intervals for plot
  result<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                        HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                        nV=length(predict(.$fit[[1]])))# lower CI
  result$AUC<-round(pracma::trapz(result$I[(which(abs(result$I-0.5)==min(abs(result$I-0.5)))-1):(which(abs(result$I-0.5)==min(abs(result$I-0.5)))+1)]),2)
  
  result <-result %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]]))
  
  
  nlm1<-list()
  CT<-list()
  dfc<-list()
  
  #sigmoidal fit for treated
  nlm2<-df_1 %>%  
    dplyr::mutate(fit=list(try(fitSingleSigmoid(C,I))))
  nlm2$Tm<-Tm(nlm2)
  
  if(class(nlm2$fit[[1]])=='try-error'){
    return(warning("the sigmoidal function could not fit the treated data"))
  }
  #remove proteins with a plateau value not =  zero 
  CT<-nlm2 %>%
    dplyr::filter( !coef(fit[[1]])[[1]]==0)#find values where Pl = 0 
  if(nrow(CT)==0){
    warning("error with the plateau")
    dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95))
    nlm2<- nlm2  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
    
  }else{
    dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95)) #predicted values
    #obtain sigmoidal equation parameters
    nlm2<- CT  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
    
  }
  
  #get confidence intervals for sigmoidal function
  #ready confidence intervals for plot
  result1<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                         HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                         nV=length(predict(.$fit[[1]])))# lower CI
  result1$AUC<-round(pracma::trapz(result1$I[(which(abs(result1$I-0.5)==min(abs(result1$I-0.5)))-1):(which(abs(result1$I-0.5)==min(abs(result1$I-0.5)))+1)]),2)
  
  result1<-result1 %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]]))
  
  #remove fit column
  Pred<- result %>% dplyr::select(-fit)
  Pred1<- result1 %>% dplyr::select(-fit)
  
  return(rbind(Pred,Pred1))
}

sigfit<-function(SigF,Peptide=FALSE){
  
  #
  if(isTRUE(Peptide)){
    Pred<-SigF %>%
      subset(dataset=="vehicle") %>% 
      dplyr::select(uniqueID,dataset ,C,I,CC,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,sample_name,sample,missing_pct)
    Pred1<-SigF%>%
      subset(dataset=="treated") %>% 
      dplyr::select(uniqueID,dataset ,C,I,CC,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,sample_name,sample,missing_pct)
    Pred$LOW<-Pred$LOW[[1]]
    Pred$HI<-Pred$HI[[1]]
    
    Pred1$LOW<-Pred1$LOW[[1]]
    Pred1$HI<-Pred1$HI[[1]]
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    Pred1$dAUC<-abs(as.double(round(Pred1$AUC[1]-Pred$AUC[1],2)))
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=dataset))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
      ggplot2::annotate("text", x=min(Pred1$C)+3, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+3, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+3, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      annotate("text",
               x = round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = round(Pred$Tm[1],1), xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)+ylim(-0.6,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
    P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 )+
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = round(min(Pred$C),1), xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2)+theme(legend.position="bottom")
    
    
    print(P1)
  }else{
    Pred<-SigF %>%
      subset(dataset=="vehicle") %>% 
      dplyr::select(uniqueID,dataset,C,I,CC,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,sample_name,sample,missing_pct,rank)
    Pred1<-SigF[[i]]%>%
      subset(dataset=="treated") %>% 
      dplyr::select(uniqueID,dataset,C,I,CC,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,sample_name,sample,missing_pct,rank)
    Pred$LOW<-Pred$LOW[[1]]
    Pred$HI<-Pred$HI[[1]]
    
    Pred1$LOW<-Pred1$LOW[[1]]
    Pred1$HI<-Pred1$HI[[1]]
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    Pred1$dAUC<-as.double(round(Pred1$AUC[1]-Pred$AUC[1],2))
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=dataset))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=-0.65, label= paste("missing vehicle",round(Pred1$missing_pct[1],0),"%"))+ 
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=-0.75, label= paste("missing treated",round(Pred1$missing_pct[1],0),"%"))+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = round(Pred$Tm[1],1), xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 )+
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = round(min(Pred$C),1), xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2) +ylim(-0.8,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
    print(P1)
    
    
  }
  
}

#plot global histograms for variability 
his_sp<-function(Df1,df.temps,MD=FALSE){
  if(any(names(Df1)=="C")){
    Df1<-dplyr::bind_rows(Df1) 
    Df1<-Df1%>% dplyr::rename("temperature"="C")
    
  }
  if(any(names(dplyr::bind_rows(Df1))=="sample_name.x")){
    test<-dplyr::bind_rows(Df1) %>% dplyr::rename("sample_name"=ifelse(any(names(.)=="sample_name.x"),"sample_name.x","sample_name.y"))
    test<-test%>% dplyr::select(Dataset,CV_pct,dataset,C,sample_name) %>% unique(.)
    test<-test %>% dplyr::rename("temperature"="C")
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
    
    
    test<-test %>% dplyr::mutate(carrier=ifelse(str_detect(test$sample_name,"NOcarrier"),"+ Carrier","- Carrier"),
                                 FAIMS = ifelse(str_detect(test$sample_name,"NO_FAIMS"),"+ FAIMS","- FAIMS"),
                                 Phi = ifelse(str_detect(test$sample_name,"PhiSDM"),"+ PhiSDM","- PhiSDM"))
  }else{
    test<-dplyr::bind_rows(Df1)
    test<-test %>% unique(.)
    if (any(names(test)=="C")){
      test<-test %>% dplyr::rename("temperature"="C")
    }
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
  }
  
  
  test<-test %>% dplyr::mutate(carrier=ifelse(str_detect(test$sample_name,"NOcarrier"),"+ Carrier","- Carrier"),
                               FAIMS = ifelse(str_detect(test$sample_name,"NO_FAIMS"),"+ FAIMS","- FAIMS"),
                               Phi = ifelse(str_detect(test$sample_name,"PhiSDM"),"+ PhiSDM","- PhiSDM"))
  
  if(isTRUE(MD)){
    test$Dataset<-as.factor(test$Dataset)
    
    ggplot(test,aes(y=CV_pct,x=Dataset,fill=Dataset))+
      facet_grid(c(~temperature,~carrier),scales="free_x")+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,200)
  }else{  
    
    test$Dataset<-test$dataset
    ggplot(test,aes(y=CV_pct,x=Dataset,fill=Dataset))+
      facet_grid(~temperature)+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,200)
    
    
  }
}
#plot global histograms for variability 
Violin_panels<-function(df_raw,df.temps,MD=TRUE){
  if(any(names(df_raw)=="C")){
    Df1<-dplyr::bind_rows(df_raw) 
    Df1<-Df1%>% dplyr::rename("temperature"="C")
    
  }
  if(any(names(dplyr::bind_rows(Df1))=="sample_name.x")){
    test<-dplyr::bind_rows(Df1) %>% dplyr::rename("sample_name"=ifelse(any(names(.)=="sample_name.x"),"sample_name.x","sample_name.y"))
    test<-test%>% dplyr::select(Dataset,CV_pct,dataset,C,sample_name) %>% unique(.)
    test<-test %>% dplyr::rename("temperature"="C")
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum_File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$dataset<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum_File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum_File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum_File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum_File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum_File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum_File,"S_Phi")==TRUE,"S",NA)))
    test<-test %>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                 "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                 "I_Interference"="Isolation_Interference_[%]")
    
    
  }else{
    test<-df_raw
    test<-test %>% dplyr::left_join(df.temps,by="temp_ref")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum_File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$dataset<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum_File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum_File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum_File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum_File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum_File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum_File,"S_Phi")==TRUE,"S",NA)))
    test<-test %>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                 "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                 "I_Interference"="Isolation_Interference_[%]")
    
    
    
  }
  
  
  if(isTRUE(MD)){
    test$sample_name<-as.factor(test$sample_name)
    
    P<-ggplot(test,aes(y=I_Interference,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("Isolation Interference (%)")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,100)+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=IonInjTime,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("Ion Injection Time (ms)")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,75)+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=PEP,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      scale_fill_brewer(palette = "Dark2") +
      geom_boxplot(width=0.1) +
      ggplot2::ylab("PEP Score")+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      guides(fill=guide_legend(title="Method"))
    P<-ggplot(test,aes(y=DeltaM,x=sample_name,fill=sample_name))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      scale_fill_brewer(palette = "Dark2") +
      ggplot2::ylab(expression(paste(Delta, " M (ppm)")))+
      ggplot2::xlab("Method")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(-15,15)+
      guides(fill=guide_legend(title="Method"))
    
  }
}
#Convert to MSStatsTMT
MSStats_converter<-function(df_raw,solvent){
  editPSMs2<-data.frame()
  editPSMs2<-df_raw %>% dplyr::rename("ProteinName"="Accession",
                                      "PeptideSequence"="Annotated_Sequence",
                                      "Run"="Spectrum_File",
                                      "PSM"="PSMs_Peptide_ID",
                                      "Channel"="temp_ref",
                                      "Intensity"="value")
  #Condition, Bioreplicate and TechRepMixture need to be filled
  
  editPSMs2$Condition<-ifelse(editPSMs2$Channel=="126","Norm",0)
  editPSMs2$BioReplicate<-ifelse(str_detect(editPSMs2$Run,solvent)=="TRUE","vehicle","treated")
  editPSMs2<-editPSMs2 %>% 
    dplyr::mutate(Mixture=paste0(ifelse(str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                 ifelse(str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                 ifelse(str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))
  
  editPSMs2$TechRepMixture<-1
  editPSMs2<-editPSMs2 %>% dplyr::select(ProteinName,PeptideSequence,Charge,PSM,Mixture,TechRepMixture,Run,Channel,Condition,BioReplicate,Intensity)
  
  Annotation<-editPSMs2 %>% dplyr::select(Run,TechRepMixture,Channel,Condition,Mixture,BioReplicate)
  Annotation$Fraction<-1
  Annotation<-Annotation %>% dplyr::group_split(Run)
}



#TPP TR Reader
TPPbenchmark<-function(f,volcano=TRUE){
  f<-"~/Files/Scripts/Files/TPP_results"
  f<-list.files(f)
  #read data
  df_TPP<-lapply(f,function(x) read_excel(x,.name_repair = "unique"))
  #extract experiment names
  f<-str_extract_all(f,c("C_F_E","C_F_S","C_nF_E","C_nF_S","nC_F_E","nC_F_S","nC_nF_E","nC_nF_S"))
  f<-lapply(f,function(x) data.frame(sample_name=as.factor(x)))
  
  #join data
  df_TPP<-purrr::map2(df_TPP,f,function(x,y)cbind(x,y))
  
  
  #select columns of interest
  df_TPP<-dplyr::bind_rows(df_TPP) %>% 
    select(Protein_ID,fulfills_all_4_requirements,
           model_converged_DMSO_1,model_converged_DMSO_2,model_converged_MEKi_1,model_converged_MEKi_2,
           sample_name,
           diff_meltP_MEKi_1_vs_DMSO_1,diff_meltP_MEKi_2_vs_DMSO_2,
           pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>% dplyr::filter(model_converged_DMSO_1=="Yes" & model_converged_DMSO_2 =="Yes" ,model_converged_MEKi_2=="Yes" & model_converged_MEKi_1 =="Yes")
  
  df_TPP$Protein_ID<-as.factor(df_TPP$Protein_ID)
  df_TPP$sample_name<-str_replace(df_TPP$sample_name,"S","\u03A6")
  #get names
  
  df_TPP<-df_TPP %>% dplyr::group_split(sample_name)
  if(isTRUE(volcano)){
    
    df_TPP1<-purrr::map(df_TPP,function(x)x %>% dplyr::select(Protein_ID,sample_name,pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2,diff_meltP_MEKi_1_vs_DMSO_1,diff_meltP_MEKi_2_vs_DMSO_2) %>% 
                          pivot_longer(c(pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2),
                                       names_to = c("hi","dTm"),
                                       names_pattern = c("(.+)pVal_(.+)"),
                          ) %>% dplyr::rename("p_dTm"="value"))
    df_TPP2<-purrr::map(df_TPP1,function(x)x %>% dplyr::select(Protein_ID,sample_name,diff_meltP_MEKi_1_vs_DMSO_1,diff_meltP_MEKi_2_vs_DMSO_2,p_dTm) %>% 
                          pivot_longer(c(diff_meltP_MEKi_1_vs_DMSO_1,diff_meltP_MEKi_2_vs_DMSO_2),
                                       names_to = c("hi","set"),
                                       names_pattern = "(.+)diff_(.+)"
                          ) %>% dplyr::rename("dTm"="value") %>% dplyr::select(-hi,-set))
    df_TPP3<-purrr::map(df_TPP2,function(x) x[!is.na(x$dTm),])
    df_TPP2<-dplyr::bind_rows(df_TPP3) 
    df_TPP2$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_TPP2$diffexpressed[df_TPP2$dTm > 1 & df_TPP2$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_TPP2$diffexpressed[df_TPP2$dTm < -1 & df_TPP2$p_dTm < 0.05] <- "Destabilized"
    df_TPP2$delabel <- NA
    df_TPP2$delabel[df_TPP2$Protein_ID %in%c("P36507","Q02750")] <- as.character(df_TPP2$Protein_ID[df_TPP2$Protein_ID %in%c("P36507","Q02750")])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name,Protein_ID)
    df_TPP2<-purrr::map(df_TPP2,function(x) x[1,])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name)
    df_TPP2<-purrr::map(df_TPP2,function(x)
      ggplot(data=x,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-1, 1), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=c("blue", "black", "red"))+
        labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
        geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
        geom_text_repel(aes(dTm, -log10(p_dTm),label=delabel))+ggtitle(x$sample_name[1])+
        theme(legend.position="bottom", legend.box = "horizontal")+ylim(-0.1,55))
    
    
    return(df_TPP2)
  }
  #for upset plots
  df_1<-purrr::map(df_TPP,function(x)
    pivot_wider(x,names_from=pVal_adj_MEKi_1_vs_DMSO_1,values_from=c(pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2),values_fill=NA))
  
  # df_1<-dplyr::bind_rows(df_1)
  df_1<-dplyr::bind_rows(df_TPP)
  
  df_1$model_converged_DMSO_1[df_1$model_converged_DMSO_1=="Yes"]<-"TRUE"
  df_1$model_converged_DMSO_2[df_1$model_converged_DMSO_2=="Yes"]<-"TRUE"
  df_1$model_converged_MEKi_1[df_1$model_converged_MEKi_1=="Yes"]<-"TRUE"
  df_1$model_converged_MEKi_2[df_1$model_converged_MEKi_2=="Yes"]<-"TRUE"
  df_1$fulfills_all_4_requirements[df_1$fulfills_all_4_requirements=="Yes"]<-"TRUE"
  
  df_1$model_converged_DMSO_1[df_1$model_converged_DMSO_1=="No"]<-"FALSE"
  df_1$model_converged_DMSO_2[df_1$model_converged_DMSO_2=="No"]<-"FALSE"
  df_1$model_converged_MEKi_1[df_1$model_converged_MEKi_1=="No"]<-"FALSE"
  df_1$model_converged_MEKi_2[df_1$model_converged_MEKi_2=="No"]<-"FALSE"
  df_1$fulfills_all_4_requirements[df_1$fulfills_all_4_requirements=="No"]<-"FALSE"
  
  df_1$sample_name<-str_replace(df_1$sample_name,"S","\u03A6")
  df_1<-lapply(df_1,function(x) x %>% dplyr::rename("Passed_4_filters"="fulfills_all_4_requirements"))
  df_1<-dplyr::bind_rows(df_1) 
  df_1$Passed_4_filters<-as.factor(df_1$Passed_4_filters)
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1<-df_1 %>% dplyr::group_split(sample_name)
  
  check<-list()
  check<-purrr::map(df_1,function(x) upset(x,colnames(x)[!colnames(x) %in% c("sample_name","Protein_ID")],
                                           min_degree=2,
                                           set_sizes=FALSE,
                                           guides='collect',
                                           n_intersections=5,
                                           stripes='white',
                                           sort_intersections_by="cardinality",
                                           stat='count',
                                           position=position_fill(vjust = .5)
                                           # show_hide_scale,
                                           # rating_scale
                                           # 
  )+ggtitle(paste0("TPP results: ",x$sample_name[1])))
  
  rating_scale = scale_fill_manual(name="High Quality Melt Curves",
                                   values=c("TRUE" ='#fee6ce', "FALSE" ='#fdae6b', "NA"  = '#e6550d'))
  
  show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
  
  check<-list()
  check<- purrr::map(df_1,function(x)upset(x,names(x),
                                           min_degree=1,
                                           set_sizes=FALSE,
                                           guides='collect',
                                           n_intersections=5,
                                           height_ratio = 0.7,
                                           stripes='white',
                                           base_annotations=list(
                                             '# of Protein Curves'=intersection_size(
                                               counts=TRUE,
                                               aes=aes(fill=sample_name,
                                                       color=sample_name,
                                                       palette="Dark2")
                                             ) 
                                             
                                           ),
                                           annotations =list(
                                             'Bin %'=list(
                                               aes=aes(x=intersection, fill=Passed_4_filters),
                                               geom=list(
                                                 geom_bar(stat='count', position='fill', na.rm=TRUE),
                                                 geom_text(
                                                   aes(
                                                     label=!!aes_percentage(relative_to='intersection'),
                                                     color="black"
                                                   ),
                                                   stat='count',
                                                   position=position_fill(vjust = .5)
                                                 ),
                                                 scale_y_continuous(labels=scales::percent_format()),
                                                 
                                               )
                                             )
                                           )
                                           
                                           
  )+ggtitle(paste0("TPP results: ",x$sample_name[1])))
  #y<-get_legend(check$patches$plots[[1]])
  colors<-as.list(brewer.pal(n = 8, name = "Dark2"))
  test<-purrr::map2(check,colors,function(x,y) x+scale_fill_manual(palette))
  data<-unlist(lapply(check,function(x) x$labels$title))
  check<-check[order(data)]
  P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO")
  
  return(print(P))
  
  
}
#CV<-his_sp(res_sp[[3]],df.temps,MD=FALSE)
UpSet_curves<-function(f,Trilinear=FALSE,Splines=FALSE,Sigmoidal=TRUE,Peptide=FALSE,filter=FALSE){
  
  if(isTRUE(Peptide) & any(names(f)=="rank_l")){
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=as.factor(sample_name),dataset=as.factor(dataset)) %>%
      dplyr::select(-rank,-rank_l,-C,-I,-temp_ref,-CV_pct,-missing_pct) %>%
      dplyr::group_split(uniqueID,dataset,sample_name,sample)
    
    f<-purrr::map(f,function(x) x %>% group_by(sample,dataset) %>% dplyr::summarise(uniqueID=uniqueID,
                                                                                    dataset=dataset,
                                                                                    sample_name=sample_name,
                                                                                    M1=M1,
                                                                                    sample=sample,
                                                                                    Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                                                                    rss=mean(rss,na.rm=TRUE),
                                                                                    rsq=mean(rsq,na.rm=TRUE),
                                                                                    AUC=mean(AUC,na.rm=TRUE)) %>% 
                    ungroup(.) %>% distinct(.))
  }else if(!isTRUE(Peptide)){
    f<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=sample_name,dataset=as.factor(dataset)) %>% 
      dplyr::select(-rank,-C,-I,-temp_ref,-CV_pct,-missing_pct) %>%
      dplyr::group_split(uniqueID,dataset,sample_name,sample)
    
    f<-purrr::map(f,function(x) x %>%
                    group_by(sample,dataset) %>%
                    dplyr::summarise(uniqueID=uniqueID,
                                     dataset=dataset,
                                     sample_name=sample_name,
                                     M1=M1,
                                     sample=sample,
                                     Coverage=Coverage,
                                     MW_kDa=MW_kDa,
                                     Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                     rss=mean(rss,na.rm=TRUE),
                                     rsq=mean(rsq,na.rm=TRUE),
                                     AUC=mean(AUC,na.rm=TRUE)) %>%
                    ungroup(.) %>% distinct(.)) 
  }else if (isTRUE(Trilinear)){#if this is a trilinear result
    #f<-f %>% dplyr::group_split(uniqueID,dataset)
    # f<-lapply(f,function(x) dplyr::bind_rows(x))
    # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=sample_name,dataset=as.factor(dataset)) %>%
      dplyr::select(-rsq,-CI,-data) %>%
      dplyr::group_split(uniqueID,dataset,sample_name)
    
    
    f<-purrr::map(f,function(x) x %>% dplyr::summarise(uniqueID=uniqueID,
                                                       dataset=dataset,
                                                       M1=M1,
                                                       Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                                       rss=sum(rss,na.rm=TRUE),
                                                       rsq=mean(rsq,na.rm=TRUE),
                                                       AUC=mean(AUC,na.rm=TRUE),
                                                       rssDiff=rssDiff,
                                                       sample_name=x$sample_name,
                                                       Fvals=Fvals,
                                                       rssDiff=rssDiff,
                                                       pV=pV,
                                                       pAdj=pAdj) %>% 
                    ungroup(.) %>% distinct(.))
  }
  
  
  if(isTRUE(Trilinear)){
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="lm"))
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$dataset=="treated"][1]>x$Tm[x$dataset=="vehicle"][1],"Stabilized","Destabilized")),
                                                    dTm=x$Tm[x$dataset=="treated"][1]-x$Tm[x$dataset=="vehicle"][1]))
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    
    df_1<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=as.character(sample_name)) %>% dplyr::group_split(uniqueID,sample_name)
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,model_converged,stabilized,rsq_greater_than_0.8) %>% 
      pivot_wider(names_from=sample_name,values_from=c(model_converged)) %>% distinct(.)
    
  }else if(isTRUE(Splines)){
    f<-dplyr::bind_rows(f) %>% dplyr::select(-sample) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="gam"))
    
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$dataset=="treated"][1]>x$Tm[x$dataset=="vehicle"][1],"Stabilized","Destabilized")),
                                                    dTm=x$Tm[x$dataset=="treated"][1]-x$Tm[x$dataset=="vehicle"][1]))
    
    f<-purrr::map(f,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(any(class(M1[[1]])=="gam"),1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    
    
    df_TPP<-dplyr::bind_rows(f) %>% 
      select(uniqueID,sample_name,Tm,rss,AUC,rsq,rsq_greater_than_0.8,model_converged,stabilized) %>% 
      distinct(.)
    
    df_1<-df_TPP %>% dplyr::select(uniqueID,sample_name,Tm,model_converged,rsq_greater_than_0.8,stabilized) %>% distinct(.)
    
    df_1<-dplyr::bind_rows(df_1)%>% dplyr::mutate(sample_name=as.character(sample_name)) %>% dplyr::group_split(uniqueID,sample_name) 
    
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,model_converged,stabilized) %>% 
      pivot_wider(names_from=sample_name,values_from=c(model_converged)) %>% distinct(.)
    
    
  }
  
  df_<-df_ %>% dplyr::mutate(uniqueID=as.character(uniqueID),stabilized=as.character(stabilized))
  df_ <- df_ %>%
    mutate_if(sapply(df_, is.factor), as.numeric)
  
  #keep the first row out of redundant peptide group data
  
  #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$dataset=="treated"]>x$Tm[x$dataset=="vehicle"],1,0)))
  # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,dataset) %>% 
  #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
  #                           RSQ=mean(Rsq,na.rm=TRUE)))
  # f<-dplyr::bind_rows(f) 
  # f<-f %>% dplyr::mutate(rss_a=ifelse(dataset=="vehicle"|dataset=="treated",RSS,NA),
  #                        rss_n=ifelse(dataset=="null",RSS,NA))
  # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
  
  #select columns of interest
  
  
  #,rsq_greater_than_0.8,stabilized
  rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                   values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
  
  check<-list()
  check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","uniqueID","dataset","sample","Tm","dTm","rsq_greater_than_0.8","stabilized")],
               min_degree=2,
               set_sizes=FALSE,
               guides='collect',
               n_intersections=10,
               stripes='white',
               sort_intersections_by="cardinality",
               base_annotations=list(
                 '# of Protein Events'=intersection_size(
                   counts=TRUE
                 )
               ),
               # base_annotations=list(
               #   'Number of protein events'=upset_annotate(
               #     '..count..',
               #     list(
               #       geom_bar(aes(fill=df_$stabilized)),
               #       scale_fill_manual(values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
               #     )
               #   )
               # )#,
               
               # annotations =list(
               #   'Stabilized Percentage'=list(
               #     aes=aes(x=intersection, fill=as.factor(df_$stabilized)),
               #     geom=list(
               #       geom_bar(stat='count', position='fill', na.rm=TRUE),
               #       geom_text(
               #         aes(
               #           label=!!aes_percentage(relative_to='group'),
               #           group=df_$stabilized
               #         ),
               #         stat='count',
               #         position=position_fill(vjust = .5)
               #       ),
               #       scale_y_continuous(labels=scales::percent_format()),
               #       rating_scale
               #     )
               #   )
               # ,
               # 'Tm'=(
               #   # note that aes(x=intersection) is supplied by default and can be skipped
               #   ggplot(mapping=aes(y=log10(df_$dTm)))+
               #     # checkout ggbeeswarm::geom_quasirandom for better results!
               #     geom_jitter(aes(color=log2(df_$dTm), na.rm=TRUE))+
               #     geom_violin(alpha=0.5, na.rm=TRUE)
               # )
               #),
               width_ratio=0.1,
               height_ratio=0.8,
               themes=upset_default_themes(text=element_text(size=15,colour="black"))
  )+ggtitle(paste0("Number of fitted ", ifelse(isTRUE(Splines),"Spline","Trilinear") ," curves (",ifelse(isTRUE(Peptide),"Peptide","Protein"),"-level ", ifelse(isTRUE(filter),"filtered","unfiltered"),"): "))
  
  
  print(check)
  
  
}
volcano_data<-function(f,Trilinear=FALSE,Splines=FALSE,Sigmoidal=TRUE,Peptide=FALSE,filter=FALSE){
  
  if(isTRUE(Peptide) & any(names(f)=="rank_l")){
    f<-f %>% 
      dplyr::group_split(uniqueID,dataset,sample)
    
    addTm<-function(x){
      f<-dplyr::bind_rows(x) %>%
        dplyr::group_by(uniqueID,dataset,sample) %>% # group individual curve data to calculate individual Tm
        dplyr::mutate(sample_name=as.factor(sample_name),
                      dataset=as.factor(dataset),
                      Tm=try(with(.,stats::approx(I,C, xout=min(I,na.rm=TRUE)+(0.5*(max(I, na.rm=TRUE)-min(I, na.rm=TRUE))))$y))) %>% 
        dplyr::select(-rank,-C,-I,-temp_ref,-CV_pct,-missing_pct) %>% dplyr::ungroup(.)
    }
    f<-mclapply(f,addTm,mc.cores=availableCores())
    #remove data where TM cannot be calculated
    
    f<-f %>% purrr::keep(function(x) !is.null(x))
    f<-dplyr::bind_rows(f)
    f<-f %>% dplyr::group_split(uniqueID)
    f<-f %>% purrr::keep(function(x) any(x$dataset %in% c("vehicle")) & any(x$dataset%in% c("treated")))
    FC<-function(x){ x %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% 
        dplyr::mutate(v_Tm=mean(x$Tm[x$dataset == "vehicle"],na.rm=TRUE),
                      t_Tm=mean(x$Tm[x$dataset == "treated"],na.rm=TRUE),
                      dTm=try(mean(Tm[dataset == "treated"],na.rm=TRUE)-mean(Tm[dataset == "vehicle"],na.rm=TRUE)),
                      FC=log2(v_Tm/t_Tm),
                      variance_equal_vt = var.test(x$Tm ~ x$dataset)$p.value,
                      p_dTm= try(ifelse(all(x$variance_equal_vt < 0.05),
                                        t.test(x$Tm ~ x$dataset, data = x,
                                               var.equal = ifelse(all(x$variance_equal_vt<0.05),FALSE,TRUE))$p.value,NA)))
    }
    
    f<-mclapply(f,FC,mc.cores=availableCores())
    ######
    f<-f %>% dplyr::group_split(uniqueID,sample_name)
    f<-purrr::map(f,function(x) x[1,])
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID)
    f<-purrr::map(f,function(x) x %>% dplyr::ungroup(.) %>% group_by(sample,dataset) %>% 
                    dplyr::mutate(uniqueID=uniqueID,
                                  dataset=dataset,
                                  sample_name=sample_name,
                                  M1=M1,
                                  sample=sample,
                                  Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                  rss=mean(rss,na.rm=TRUE),
                                  rsq=mean(rsq,na.rm=TRUE),
                                  AUC=mean(AUC,na.rm=TRUE)
                    ) %>% 
                    ungroup(.) %>% distinct(.))
    
  }else if(!isTRUE(Peptide)){
    f<-f %>% 
      dplyr::group_split(uniqueID,dataset,sample)
    addTm<-function(x){
      f<-x %>%
        dplyr::group_by(uniqueID,dataset,sample) %>% # group individual curve data to calculate individual Tm
        dplyr::mutate(sample_name=as.factor(sample_name),
                      dataset=as.factor(dataset),
                      Tm=with(.,stats::approx(I,C, xout=min(I,na.rm=TRUE)+(0.5*(max(I, na.rm=TRUE)-min(I, na.rm=TRUE))))$y)) %>% 
        dplyr::select(-rank,-C,-I,-CV_pct,-missing_pct) %>% dplyr::ungroup(.)
    }
    f<-mclapply(f,addTm,mc.cores=availableCores())
    f<-dplyr::bind_rows(f)
    f<-f %>% dplyr::group_split(uniqueID)
    f<-f %>% purrr::keep(function(x) any(x$dataset %in% c("vehicle")) & any(x$dataset%in% c("treated")))
    FC<-function(x){ x %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% 
        dplyr::mutate(v_Tm=mean(x$Tm[x$dataset == "vehicle"],na.rm=TRUE),
                      t_Tm=mean(x$Tm[x$dataset == "treated"],na.rm=TRUE),
                      dTm=try(mean(Tm[dataset == "treated"],na.rm=TRUE)-mean(Tm[dataset == "vehicle"],na.rm=TRUE)),
                      FC=log2(v_Tm/t_Tm),
                      variance_equal_vt = var.test(x$Tm ~ x$dataset)$p.value,
                      p_dTm= try(ifelse(all(x$variance_equal_vt < 0.05),
                                        t.test(x$Tm ~ x$dataset, data = x,
                                               var.equal = ifelse(all(x$variance_equal_vt<0.05),FALSE,TRUE))$p.value[1],NA)))
    }
    
    f<-mclapply(f,FC,mc.cores=availableCores())
    
  }else if (isTRUE(Trilinear)){#if this is a trilinear result
    #f<-f %>% dplyr::group_split(uniqueID,dataset)
    # f<-lapply(f,function(x) dplyr::bind_rows(x))
    # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
    f<-dplyr::bind_rows(f) %>%
      dplyr::group_by(uniqueID,dataset,sample) %>% # group individual curve data to calculate individual Tm
      dplyr::mutate(sample_name=as.factor(sample_name),
                    dataset=as.factor(dataset),
                    Tm=with(.,stats::approx(I,C, xout=min(I,na.rm=TRUE)+(0.5*(max(I, na.rm=TRUE)-min(I, na.rm=TRUE))))$y)) %>% 
      dplyr::select(-rsq,-CI,-data) %>% dplyr::ungroup(.)
    f<-f %>% 
      dplyr::group_split(uniqueID)
    FC<-function(x){ x %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% 
        dplyr::mutate(v_Tm=mean(x$Tm[x$dataset == "vehicle"],na.rm=TRUE),
                      t_Tm=mean(x$Tm[x$dataset == "treated"],na.rm=TRUE),
                      dTm=mean(Tm[dataset == "treated"],na.rm=TRUE)-mean(Tm[dataset == "vehicle"],na.rm=TRUE),
                      FC=log2(v_Tm/t_Tm),
                      variance_equal_vt = var.test(x$Tm ~ x$dataset)$p.value,
                      p_dTm= ifelse(all(x$variance_equal_vt < 0.05),
                                    t.test(x$Tm ~ x$dataset, data = x,
                                           var.equal = ifelse(all(x$variance_equal_vt<0.05),FALSE,TRUE))$p.value,NA))
    }
    
    f<-mclapply(f,FC,mc.cores=availableCores())
  }
  if(isTRUE(Trilinear)){
    f<-f %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="lm"))
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$dataset=="treated"][1]>x$Tm[x$dataset=="vehicle"][1],"Stabilized","Destabilized"))))
    
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    df_<-dplyr::bind_rows(f) %>% 
      select(uniqueID,FC,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,rsq_greater_than_0.8,model_converged,stabilized) %>% 
      distinct(.)
    df_$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 1 & df_$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -1 & df_$p_dTm < 0.05] <- "Destabilized"
    df_$delabel <- NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750")])
    
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-1, 1), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=c("blue", "black", "red"))+
      labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
      geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_)+
      geom_text_repel(aes(dTm, -log10(p_dTm),label=delabel))+ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")
    
  }else if(isTRUE(Splines)){
    f<-f %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
    f<-dplyr::bind_rows(f) %>% dplyr::select(-sample,-variance_equal_vt,-v_Tm,-t_Tm,-k_,-CC) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="gam"))
    
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(dTm>0,"Stabilized",ifelse(dTm<0,"Destabilized","No")))))
    
    f<-purrr::map(f,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(any(class(M1[[1]])=="gam"),1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    
    
    df_<-dplyr::bind_rows(f) %>% 
      select(uniqueID,FC,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,rsq_greater_than_0.8,model_converged,stabilized) %>% 
      distinct(.)
    df_$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 1 & df_$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -1 & df_$p_dTm < 0.05] <- "Destabilized"
    df_$delabel <- NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750")])
    
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-1, 1), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=c("blue", "black", "red"))+
      labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
      geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
      geom_text_repel(aes(dTm, -log10(p_dTm),label=delabel))+ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")
    
    
  }
  
  return(check)
  
}

###############################
memory.limit(175921900000)#set for 16 GB RAM
plan(multicore,workers=availableCores())
##################################

# #prepare a list of proteins
# setwd("~/Cliff prot pep")

#rename to the folder where your PSM file is located
f<- list.files(pattern='*PEPTIDES2.xlsx')
f<- list.files(pattern='*PROTEINS.xlsx')

#New
f<- list.files(pattern='*PSMs.xlsx')
f<- list.files(pattern='*_0.xlsx')
f<- list.files(pattern='*Proteins.xlsx')
#Covid
# f<- list.files(pattern='*Proteins.xlsx')
# f<- list.files(pattern='*PSMs.xlsx')


# df_raw<-df_raw %>% dplyr::left_join(df.samples,by=c("temp_ref",))


# PSMs<-read_excel(f)
# PSMs<-PSMs %>% dplyr::rename("uniqueID"="Protein","sample_name"="Mixture","dataset"="BioReplicate","temp_ref"="Channel","I"="Abundance")
# PSMs<-PSMs %>% dplyr::select(uniqueID,sample_name,dataset,temp_ref,I)
# PSMs<-PSMs %>% dplyr::mutate(sample_name=ifelse(!is.na(str_match(PSMs$sample_name,'[:punct:][:digit:][:digit:][:digit:][:punct:][:digit:]')),str_match(PSMs$sample_name,'[:punct:][:digit:][:digit:][:digit:][:punct:][:digit:]'),PSMs$sample_name)) 
# PSMs<-PSMs %>% dplyr::mutate(sample_name=str_replace(PSMs$sample_name," ","_")) %>%
#   dplyr::left_join(df.temps,by="temp_ref") %>%
#   dplyr::rename("C"="temperature") %>% 
#   dplyr::select(-temp_ref)
# d<-PSMs
# d<-d %>% dplyr::mutate(CC=ifelse(stringr::str_detect(d$sample_name,"DMSO"),0,1),
#                        rank= dplyr::ntile(I,3),
#                        sample=sample_name)
#rename to the folder where your Protein file is located
# f<-"~/Cliff prot pep/Proteins.xlsx"
# f<-"C:/Users/figue/OneDrive - Northeastern University/CETSA R/CP_Exploris_20200811_DMSOvsMEKi_carrier_FAIMS_PhiSDM_PEPTIDES.xlsx"
#df_raw <- read_cetsa("~/Files/Scripts/Files/PSM_validator","~/Files/Scripts/Files/PSM_validator","_Proteins",Peptide=FALSE,Batch=FALSE,CFS=TRUE,solvent="DMSO")     
#df_raw <- read_cetsa("~/Files/Scripts/Files/Covid","~/Files/Scripts/Files/Covid","_Proteins",Peptide=FALSE,CFS=FALSE,Batch=FALSE)                                                              
df_raw <- read_cetsa("~/Files/Scripts/Files/CONSENSUS","~/Files/Scripts/Files/CONSENSUS","_Proteins",Peptide=FALSE,Batch=FALSE)                                                              
#saveRDS(df_raw,"df_raw.RDS")

#filter Peptides
filter_Peptides<-function(df_,S_N,PEP,XCor,Is_Int,Missed_C,Mods,Charg,DeltaMppm,filter_rank=FALSE,shared_proteins=FALSE,CFS=TRUE){
  #remove shared proteins
  if(!isTRUE(shared_proteins)){
    
    df_<-df_[!str_detect(df_$Accession,";"),]
  }
  if(any(str_detect(names(df_),"PEP"))=="FALSE"){
    df_<-df_ %>% dplyr::mutate(Percolator_PEP=0)
  }
  #rename problematic headers
  check<-names(head(df_))
  ch<-str_replace_all(check,paste0("[","[:punct:]","]"),paste0("_"))
  ch<-str_replace_all(ch,paste0("[:punct:]","[:punct:]","[:punct:]"),paste0(""))
  ch<-str_replace_all(ch,paste0("[:punct:]","[:punct:]"),paste0(""))
  #set new names
  names(df_)<-ch
  #filter
  if(any("S_N" %in% names(df))==TRUE){
    df_<-df_ %>% dplyr::filter(Average_Reporter_S_N>S_N,Percolator_PEP<PEP,Charge<Charg,Missed_Cleavages<Missed_C,abs(DeltaMppm_)<DeltaM_ppm)
    df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","S_N"="Average_Reporter_S_N","PEP"="Percolator_PEP",
                              "DeltaM"="DeltaMppm_","IonInjTime"="Ion_Inject_Timems_",
                              "I_Interference"="Isolation_Interference_")
    #rank by the highest intensity channel
    rank<-df_%>% dplyr::filter(temp_ref=="126") %>% 
      dplyr::mutate(rank=dplyr::ntile(.$S_N,3)) %>% dplyr::select(uniqueID,Spectrum_File,Annotated_Sequence,Charge,rank,S_N,PEP,Missed_Cleavages,DeltaM,sample_id,Protein_value)
    #remove na values in uniqueID's
    rank<-rank %>% dplyr::filter(!is.na(uniqueID))
    #convert to data.table
    rank<-data.table(rank)
    df_<-data.table(df_)
    #get unique PSMs
    rank<-unique(rank)
    
    #Use I column to distinguish duplicated rank values/PSM
    col_n<-c("uniqueID","Spectrum_File","Annotated_Sequence","Charge","S_N","PEP","Missed_Cleavages","DeltaM","sample_id","Protein_value")
    #set join columns
    setkeyv(rank,cols=col_n)
    setkeyv(df_,cols=col_n)
    #right_join
    df_ <- data.frame(merge(df_,rank, all.y=TRUE))
    
    
    
  }else{
    
    #df_<-df_ %>% dplyr::rename("Missed_Cleavages"="#_Missed_Cleavages") 
    df_<-df_ %>% dplyr::filter(Percolator_PEP<PEP,Charge<Charg,MissedCleavages<Missed_C,abs(DeltaM)<DeltaMppm)
    df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","PEP"="Percolator_PEP")
    
    
    if(length(XCor)==2){
      df_<-df_ %>% dplyr::filter(XCorr>XCor[1],XCorr<XCor[2])
    }else{
      df_<-df_ %>% dplyr::mutate(XCor_l=ifelse(Charge==2 & XCorr > 1.8,TRUE,ifelse(Charge>2 & XCorr > XCor,TRUE,FALSE)))
    }
    
    #remove the carrier channel 
    rank<-df_ %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(I,3)) %>% dplyr::select(-temp_ref,-I)
    #remove na values in uniqueID's
    rank<-rank %>% dplyr::select(Annotated_Sequence,uniqueID,Modifications,sample_id,sample_name,dataset,rank)
    #convert to data.table
    rank<-data.table(rank)
    df_<-data.table(df_)
    #get unique PSMs
    rank<-unique(rank)
    
    #Use I column to distinguish duplicated rank values/PSM
    col_n<-dplyr::intersect(names(rank),names(df_))
    #set join columns
    setkeyv(rank,cols=col_n)
    setkeyv(df_,cols=col_n)
    #right_join
    df_ <- data.frame(merge(df_,rank, all.y=TRUE))
    
    
    df_<-df_ %>% dplyr::mutate(rank_l=ifelse(rank==3,TRUE,FALSE)) 
    
    df_<-df_ %>% dplyr::group_split(sample_id,uniqueID)
    #ove data that has at least one high ranking intensity value and passed the XCor filters by charge state
    if(isTRUE(filter_rank)){
      df_<-df_ %>% purrr::keep(function(x) any(x$XCor_l==TRUE & x$rank_l==TRUE,na.rm=TRUE))
    }else{
      df_<-df_ %>% purrr::keep(function(x) any(x$XCor_l==TRUE & x$rank_l==FALSE,na.rm=TRUE))
      
    }
    
  }
  df_<-dplyr::bind_rows(df_)
  
  
  
  return(df_)
  
}
df_raw1<-df_raw %>% dplyr::group_split(sample_name)
df_raw1<-furrr::future_map(df_raw1,function(x) try(filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,filter_rank=FALSE,shared_proteins=TRUE,CFS=TRUE)))
df_raw1<-dplyr::bind_rows(df_raw1)%>% dplyr::group_split(sample_name)
df_raw1<-df_raw1 %>% purrr::keep(function(x) !nrow(x)==0)

#Sum abundances
Sum_Ab<-function(x){
  x<-x %>% dplyr::group_by(uniqueID,sample_id,dataset,temp_ref) %>% dplyr::mutate(I=sum(I,na.rm=TRUE)) 
  x<-x %>% distinct(.)
}
df_raw1<-mclapply(df_raw1,Sum_Ab,mc.cores=availableCores())

df_raw1<-df_raw1 %>% purrr::keep(function(x) !length(x)==0)
#new for PSMs
#df_raw<-df_raw %>% dplyr::rename("sample_name"="Spectrum_File")
#annotate protein data with missing values
MID<-df_raw[is.na(df_raw$value),]
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = unique(df_raw$temp_ref),temperature = c(40, 42.1, 43.8, 46.5, 50, 54, 57.3, 60.1, 62, 64), stringsAsFactors = FALSE)

df.t <- function(n,protein_path,sample_mapping_name=NA){
  if(!is.logical(sample_mapping_name)){
    TMT<-read_xlsx(sample_mapping_name) %>% 
      dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample"="Sample","sample_name"="MS_sample_number","dataset"="Time_point")  
    TMT$dataset<-as.factor(TMT$dataset)
  }else{
    
    TMT<-data.frame(NA)
    if(n==10){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67), stringsAsFactors = FALSE)
    }else if (n==11){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,68), stringsAsFactors = FALSE)
    }else if (n == 16){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C','132N','132C','133N','133C','134N'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,69,71,73,75,77,99), stringsAsFactors = FALSE)
    }
  }
  return(TMT)
}
df.temps<-df.t(11)
#Covid
df.temps<-df.t(16,sample_mapping_name="sample_mapping.xlsx")
# df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(67, 64, 60.4, 57.1, 53.8, 50.5, 47.2, 43.9, 40.6, 37.3), stringsAsFactors = FALSE)
#df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4'), sample_name = c('DMSO_1','DMSO_2', '655_1','655_2'), stringsAsFactors = FALSE)

df.s <- function(data_path,n,rep_,bio_,vehicle_name,treated_name,Batch=FALSE,PSM=FALSE){#n is df_raw, rep is tech rep
  
  if(any(names(n)=="sample_name") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,sample_name) %>% unique(.)
    return(n)
  }
  if(any(names(n)=="Spectrum_File") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,Spectrum_File) %>% unique(.)
    return(n)
  }
  if (!isTRUE(PSM)){
    if(isTRUE(Batch)){
      n<-data.frame(sample_name=unique(n$sample_name),sample_id=unique(n$sample_id))
      return(n)
    }
    f<-data_path
    find<-c('[:upper:][[:digit:]]+')
    check<-list()
    check<-purrr::map(seq(f),function(x){
      data.frame(sample_name= str_extract_all(names(read_xlsx(f[x]))[str_detect(names(read_xlsx(f[x])),"Abundance")],find)[[1]])
    })
    
    df.samples<-data.frame(sample_name=data.frame(f),sample_id=dplyr::bind_rows(check) )
    names(df.samples)<-c("sample_name","sample_id")
    return(df.samples)
  }
  
  b<-2*rep_*bio_
  samples<-data.frame(sample_name=c(paste0(rep(as.factor(vehicle_name),rep_*bio_)),paste0(rep(as.factor(treated_name),rep_*bio_))),
                      sample_id=NA)
  samples<-rownames_to_column(samples)
  samples$sample_id<-paste0("F",samples$rowname)
  samples<-samples %>% dplyr::select(-rowname)
  samples$sample_id<-as.factor(samples$sample_id)
  samples$sample_name<-as.factor(samples$sample_name)
  return(samples)
  
}

df.samples<-df.s(f,dplyr::bind_rows(df_raw),3,2,"DMSO","TREATED",Batch=TRUE,PSM=TRUE)
#Peptides
df_raw<-df_raw %>% group_split(sample_name)
df_clean <- furrr::future_map(df_raw,function(x) clean_cetsa(x, temperatures = df.temps, samples = df.samples,Peptide=FALSE,solvent="DMSO",CFS=TRUE,CARRIER=TRUE))#assgns temperature and replicate values

#Covid data
#df_clean<-purrr::map(seq_along(df_clean),function(x) rbind(df_clean[[1]],x))

#df_clean<-dplyr::bind_rows(df_clean) %>% dplyr::group_split(sample_name)

#normalize data
df_norm <- furrr::future_map(df_clean,function(x) normalize_cetsa(x, df.temps$temperature,Peptide=FALSE,filters=FALSE)) #normalizes according to Franken et. al. without R-squared filter


# rm(df_raw,df_clean)

Int_plot<-function(df_norm,Peptide=FALSE){
  df_norm$sample_name<-f$sample_name<-str_replace(df_norm$sample_name,"S","\u03A6")
  if(!isTRUE(Peptide)){
    list<-ggplot2::ggplot(df_norm,mapping=aes(x=C,y=I))+
      geom_jitter(position=position_jitter(2),alpha=0.5)+
      geom_boxplot(mapping=aes(color=as.factor(C)))+xlab('Temperature (\u00B0C)')+
      ylab("Normalized intensity protein")+
      ggtitle(df_norm$sample_name[1])+
      ylim(-0.1,10)+
      theme(legend.position="bottom")+ labs(colour = "Temperature (\u00B0C)")
  }else{
    list<-ggplot2::ggplot(df_norm,mapping=aes(x=C,y=I3))+
      geom_jitter(position=position_jitter(2),alpha=0.5)+
      geom_boxplot(mapping=aes(color=as.factor(C)))+xlab('Temperature (\u00B0C)')+
      ylab("Normalized intensity (top 25%)")+
      ggtitle(df_norm$sample_name[1])+
      ylim(-0.1,5)+
      theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")
  }
}
plot_I<-purrr::map(df_norm1,function(x) Int_plot(x,Peptide=FALSE))
check<-ggplot2::ggplot_build(plot_I[[2]])
y<-get_legend(check$plot)
data<-unlist(lapply(plot_I,function(x) x$labels$title))
plot_I<-plot_I[order(data)]
P<-ggarrange(plotlist=plot_I,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
pdf("Intensity_values_Peptide_filt.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P
dev.off()
##Generate upset plots for missing value data###
df_<-df_clean %>% dplyr::rename("sample_id"="sample")%>% dplyr::select(-missing_pct,-value,-missing,-rank)

df_<-data.frame()

df_<-df_raw%>% dplyr::right_join(df_,by=c("Accession","sample_id"))


##
listUP<-upMV(df_,"C_F_S",5,plot_multiple=TRUE,PSMs=FALSE)
pdf("UpsetMV.pdf",pointsize= 14,paper= "a4r", width = 0.001, height = 0.001,)
listUP

dev.off()

df_norm1<-df_norm

#df_norm<-purrr::map(df_norm1,function(x) x[str_detect(x$uniqueID,c("P36507","Q02750")),])

# df_norm<-dplyr::bind_rows(df_norm) %>% dplyr::group_split(uniqueID)
# df_norm<-df_norm %>% purrr::keep(function(x) nrow(x)>1)
df_norm<-purrr::map(df_norm1,function(x)x %>% dplyr::filter(uniqueID %in% c("P36507","Q02750")))


PlotTrilinear<-function(df_norm,target,df.temps,Ft,filt,Peptide=FALSE,show_results=FALSE){
  df_norm$CC<-ifelse(df_norm$dataset=="vehicle",0,1)
  if(isTRUE(Peptide) & any(names(df_norm)=="XCorr")){
    #remove columns not needed for curve fitting 
    df_norm<-df_norm %>% dplyr::select(-XCorr,-temp_ref,-Modifications,-MissedCleavages,-DeltaM,-"Annotated_Sequence",-tidyr::contains("PEP"),-Charge)
    if(any(names(df_norm)=="I_Interference")){
      df_norm<-df_norm %>% dplyr::select(-I_Interference,-IonInjTime,-S_N,-Spectrum_File)
    }
    
  }
  if(any(names(df_norm)=="I3")){
    df_norm<-df_norm %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
  }
  ##SCRIPT STARTS HERE
  DF<-df_norm %>% dplyr::group_split(uniqueID) #split null dataset only by protein ID
  d_<-df_norm %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,dataset) #split vehicle dataset
  d_1<-df_norm %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,dataset) #split treated dataset
  if(length(d_1)==0){
    d_1<-d_
  }
  if(length(d_)==0){
    warning(paste0("No vehicle curves found for ",df_norm$sample_name[1]))
  }
  #convert to data frame for uniqueID presence
  DF<-dplyr::bind_rows(DF)
  d_<-dplyr::bind_rows(d_)
  d_1<-dplyr::bind_rows(d_1)#1 
  #make sure uniqueIDs are present for treated and vehicle
  CID<-intersect(d_1$uniqueID,d_$uniqueID)
  CID<-intersect(CID,DF$uniqueID)
  DF<-DF %>% subset(uniqueID %in% CID)
  DF$LineRegion<-1
  d_<-d_%>% subset(uniqueID %in% CID)
  d_$LineRegion<-1
  d_1<-d_1%>% subset(uniqueID %in% CID)
  d_1$LineRegion<-1
  #split dataset into equal-sized lists
  DF<-DF %>%  dplyr::group_split(uniqueID) 
  d_<-d_ %>% dplyr::group_split(uniqueID) 
  d_1<-d_1 %>% dplyr::group_split(uniqueID) 
  
  #preallocate list
  results<-vector(mode = "list", length(d_))
  results_t<-vector(mode = "list",length(d_1))
  results_n<-vector(mode = "list",length(DF))
  
  
  results<-suppressWarnings(DLR(d_))#First guess at line regions
  results_t<-suppressWarnings(DLR(d_1))
  results_n<-suppressWarnings(DLR(DF))
  
  
  
  #reassign shared points between line regions
  df_<-suppressWarnings(purrr::map2(results,d_,function(x,y) CP(x,y,PSM=Peptide)))
  df_1<-suppressWarnings(purrr::map2(results_t,d_1,function(x,y) CP(x,y,PSM=Peptide)))
  DFN<-suppressWarnings(purrr::map2(results_n,DF,function(x,y) CP(x,y,PSM=Peptide)))
  
  
  #remove results to save space 
  rm(results,results_t,results_n,d_,d_1,DF)#10
  df_<-dplyr::bind_rows(df_) %>% dplyr::group_split(uniqueID)
  df_1<-dplyr::bind_rows(df_1)%>% dplyr::group_split(uniqueID)
  DFN<-dplyr::bind_rows(DFN)%>% dplyr::group_split(uniqueID)
  
  df_<-lapply(df_,function(x)x[order(x$C),])
  df_1<-lapply(df_1,function(x)x[order(x$C),])
  DFN<-lapply(DFN,function(x)x[order(x$C),])
  #get # identifier
  df_<-purrr::map2(df_,seq(df_),function(x,y)x %>% dplyr::mutate(N=y))
  df_1<-purrr::map2(df_1,seq(df_1),function(x,y)x %>% dplyr::mutate(N=y))
  DFN<-purrr::map2(DFN,seq(DFN),function(x,y)x %>% dplyr::mutate(N=y))
  
  #prealloate variables
  tlresults<-list()
  tlresults_PI<-list()
  #split data 
  df_<-dplyr::bind_rows(df_)
  df_1<-dplyr::bind_rows(df_1)
  DFN<-dplyr::bind_rows(DFN)
  #confidence intervals
  
  tlresults<-tlstat(DFN,df_,df_1,norm=FALSE,Filters=filt,Ftest=Ft,show_results=show_results)
  if(isTRUE(show_results)){
    
    return(tlresults)
  }
  #return filtered lists
  res<-tlf(tlresults,DFN,APfilt=FALSE,PF=FALSE)
  
  i=which(res[[1]]$uniqueID %in% target)
  plotTL1<-tlCI(i,res[[1]],res[[2]],res[[3]],overlay=TRUE,residuals=FALSE,df.temps=df.temps,PSMs=Peptide)
  
  return(plotTL1)
}

plot<-purrr::map(df_norm,function(x) try(PlotTrilinear(x,"P36507",df.temps,Ft=FALSE,filt=FALSE,Peptide=FALSE,show_results=FALSE)))

check<-ggplot2::ggplot_build(plot[[1]])
y<-get_legend(check$plot)
data<-order(unlist(lapply(plot,function(x) x$labels$title)))
plot<-plot[data]
P1<-ggarrange(plotlist=plot,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plot2<-purrr::map(df_norm,function(x) try(PlotTrilinear(x,"Q02750",df.temps,Ft=FALSE,filt=FALSE,Peptide=TRUE)))
check<-ggplot2::ggplot_build(plot2[[1]])
y<-get_legend(check$plot)
data<-order(unlist(lapply(plot,function(x) x$labels$title)))
plot2<-plot2[data]
P2<-ggarrange(plotlist=plot2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

#saveIDs filtered

pdf("Target_curves_trilinear_unfilt_unshared.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P1
P2
dev.off()

#plot Number of curves
Check<-UpSet_curves(plot,Trilinear=TRUE,Splines=FALSE,Sigmoidal=FALSE,Peptide=TRUE,filter=TRUE)
pdf("Number_of_curves_upset_trilinear_peptide.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
Check
dev.off()
#df1 <- only IDs in order desc(stability)
#df2<-original data in order  
#Df1 <- ordered spline results 
###############################
plot_Splines<-function(x,Protein,df.temps,MD=FALSE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE){
  Filters=Filters
  fT=fT
  MD=MD
  
  if(isTRUE(MD)){
    
    if(isTRUE(Peptide) & any(names(x)=="S_N")){
      
      x<-x %>% dplyr::select(-I_Interference,-IonInjTime,-S_N,-Spectrum_File)
      x<-x %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
      
    }else if(isTRUE(Peptide) & any(names(x)=="XCor_l")){
      x<-x %>% dplyr::select(-"Annotated_Sequence",-XCor_l,-Charge,-PEP,-DeltaM,-Modifications,-XCorr,-MissedCleavages) %>% unique(.)
      x<-x %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
    }else if(isTRUE(Peptide)){
      x<-x %>% dplyr::select(-"Annotated_Sequence",-Charge,-tidyr::contains("PEP"),-DeltaM,-Modifications,-XCorr,-MissedCleavages) %>% unique(.)
      x<-x %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
    }else if (any(names(x)=="I3")){
      x<-x %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
      
      
    }
    DFN<-x 
    df_<-x %>% dplyr::filter(dataset=="vehicle")
    df_1<-x %>% dplyr::filter(dataset=="treated")
    
    if(nrow(df_1)==0){
      df_1<-df_
    }
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    spresults<-spstat(DFN,df_,df_1,Ftest=fT,norm=FALSE,show_results=TRUE,filters=Filters)
    if(isTRUE(show_results)){
      return(spresults)
    }
    if(any(class(spresults)=="list")){
      res_sp<-spf(spresults[[1]],DFN,filters=FALSE)
    }else{
      res_sp<-spf(spresults,DFN,filters=FALSE)
    }
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide)
    
    return(Pred1)
  }else{
    DFN<-x %>% dplyr::filter(uniqueID %in% as.character(Protein))
    df_<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),dataset=="vehicle")
    df_1<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),dataset=="treated")
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    spresults<-spstat(DFN,df_,df_1,Ftest=fT,norm=FALSE,show_results=TRUE,filters=Filters)
    if(isTRUE(show_results)){
      return(spresults)
    }
    
    if(class(spresults)=="list"){
      res_sp<-spf(spresults[[1]],DFN,filters=Filters)
    }else{
      res_sp<-spf(spresults,DFN,filters=Filters)
    }
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05)
    return(Pred1)
  }
  
}

plotS <- furrr::future_map(df_norm,function(x) try(plot_Splines(x,"P36507; Q02750",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE)))
check<-ggplot2::ggplot_build(plotS[[1]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS,function(x) x$labels$title))
plotS<-plotS[order(data)]
P3<-ggarrange(plotlist=plotS,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
#For covid dataset
# check<-dplyr::bind_rows(df_norm) %>% dplyr::group_split(time_point)
# plotS2 <- purrr::map(check,function(x) try(plot_Splines(x,"P0DTC2",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE)))
# 
plotS2 <- purrr::map(df_norm1,function(x) try(plot_Splines(x,"Q02750",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=TRUE,Peptide=FALSE)))
check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS2,function(x) x$labels$title))
plotS2<-plotS2[order(data)]
P2<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS <- furrr::future_map(df_norm,function(x) try(plot_Splines(x,"P36507",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE)))
check<-ggplot2::ggplot_build(plotS[[1]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS,function(x) x$labels$title))
plotS<-plotS[order(data)]
P1<-ggarrange(plotlist=plotS,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)


pdf("Protein_Target_curves_MD_PD_HR_Settings_unfilt_unshared.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P1
P2

dev.off()
#plot volcano 
check<-purrr::map(plotS2,function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE)))
check1<-ggplot2::ggplot_build(check[[1]])
y<-get_legend(check1$plot)
P1<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("volcano_splines_Protein_panels.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P1
dev.off()

check<-TPPbenchmark(f,volcano=TRUE)
check1<-ggplot2::ggplot_build(check[[1]])
y<-get_legend(check1$plot)
P1<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("volcano_TPP_Protein_panels.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P1
dev.off()
#ot Number of curves
Check<-UpSet_curves(plotS2,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE)
pdf("Number_of_curves_upset_splines_Protein.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
Check
dev.off()
###################################################                                                                                                                                                                                                                                                                                                                            ##################################################
#Sigmoidal function with confidence intervals
###################################################

df_<-df_norm1

plot_Sigmoidal<-function(df_,Protein,Peptide){
  #remove duplicated intensity values
  df_<-df_ %>% distinct(.)
  PlSig<-try(sigC(df_,Protein))
  
  ID<-unique(PlSig$uniqueID)
  
  sig<- try(sigfit(PlSig,Peptide=Peptide))
  return(sig)
}
plotS2<-ggplot()
plotS2<-purrr::map(df_,function(x) try(plot_Sigmoidal(x,"P36507",Peptide=FALSE)))
check<-ggplot2::ggplot_build(plotS2[[3]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS2,function(x) x$labels$title))
plotS2<-plotS2[order(data)]
P2<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2<-ggplot()
plotS2<-purrr::map(df_,function(x) try(plot_Sigmoidal(x,"Q02750",Peptide=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[3]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS2,function(x) x$labels$title))
plotS2<-plotS2[order(data)]
P1<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("Sigmoidal_curves_MD_PD_HR_Settings_unfilt.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P2
P1

dev.off()

