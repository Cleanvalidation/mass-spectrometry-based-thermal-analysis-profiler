
# For 1D experiments, the abundance headers should be of the form “Ref_[Temperature]_[V/T Replicate]”,
# where [Temperature] is replaced with respective temperature values in decimal format, and [V/T Replicate] 
# is replaced with a string referring to vehicle/treatment  and replicate number.
#associate TMT channel with Temperature
df.t <- function(n,protein_path,sample_mapping_name=NA){
  if(!is.logical(sample_mapping_name)){
    TMT<-read_xlsx(sample_mapping_name) %>% 
      dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample"="Sample","sample_name"="MS_sample_number","dataset"="Time_point")  
    TMT$dataset<-as.factor(TMT$dataset)
  }else{
    
    TMT<-data.frame(NA)
    if(n==10){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), Temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67), stringsAsFactors = FALSE)
    }else if (n==11){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), Temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,68), stringsAsFactors = FALSE)
    }else if (n == 16){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C','132N','132C','133N','133C','134N'), Temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,69,71,73,75,77,99), stringsAsFactors = FALSE)
    }
  }
  return(TMT)
}
df.temps<-df.t(11)
#associate sample_names with dataset
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
  if(.Platform$OS.type=="windows"){
    Proteins<-furrr::future_map(f,function(x) read_xl(x))
    PSMs<-furrr::future_map(h,function(x) read_xl(x))
  }else{
    Proteins<-mclapply(f,read_xl,mc.cores = availableCores())
    PSMs<-mclapply(h,read_xl,mc.cores = availableCores())
  }
  if (isTRUE(Peptide)){ #if this is a PSM of peptide group file
    if(length(list.files(peptide_path,pattern="PeptideGroups.xlsx"))>0){
      f<-list.files(peptide_path,pattern="PeptideGroups.xlsx")
    }
    #first get Peptide file read
    if(.Platform$OS.type=="windows"){
      df.raw<-furrr::future_map(f,read_xl)
    }else{
      df.raw<-mclapply(f, read_xl, mc.cores = availableCores())
    }
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
        return(warning(paste0("Please check columns from Protein Group PD output:",setdiff(c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n"),check))))
      }
      df2<-data.table(df.raw) 
      df3<-melt(data = df2, 
                id.vars = c("Accession","Annotated_Sequence","Percolator_PEP","Charge","MissedCleavages","DeltaM","Modifications","XCorr","n"),
                measure.vars = names(df2)[grepl( "Abundance" , names( df2 ) )],
                variable.name = "id",
                value.name = "value")
      df3<-df3[, id:=as.character(id)]
      df3<-df3 %>% dplyr::group_split(n)
      
      cols_PD<-function(x){x %>% dplyr::mutate(temp_ref = unlist(str_extract(x$id,find)),
                                               sample_id= unlist(str_extract(x$id,"F[[:digit:]]+"))) %>% 
          dplyr::select(-id)
      }
      if(.Platform$OS.type=="windows"){
        df3<-furrr::future_map(df3,cols_PD)
      }else{
        df3<-mclapply(df3,cols_PD,mc.cores=availableCores())
      }
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
        df2<-df2 %>% dplyr::group_by(Accession,sample_id)
        df2<-df2 %>%  dplyr::mutate(missing_pct=(100*sum(as.numeric(is.na(value)))/length(value)))%>% ungroup(.)
        
        
      }
      
      
      check_PD<-function(x) {
        if(any("Grouped" %in% names(x))==TRUE){
          warning(paste0("please check sample file ",x$sample_name," for Abundance columns")) 
        }
      }
      if(.Platform$OS.type=="windows"){
        ch<-furrr::future_map(df2,check_PD)
      }else{
        ch<-mclapply(df2,check_PD,mc.cores=availableCores())
        
      }
    }
    
    #return(df2)
  }else{
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
#reads in the file complete with sample name & dataset
Protein_File <- read_cetsa("~/CONSENSUS","~/CONSENSUS","_Proteins",Peptide=FALSE,Batch=FALSE,CFS=TRUE,solvent="DMSO")     

#find abundance data
names<-names(Protein_File)[str_detect(names(Protein_File),"Abundance")]

paste0("Ref_[",Temperature,"]_[",dataset," ",Replicate,"]")