pacman::p_load(minpack.lm,
               rlist,
               data.table,
               knitr,
               ggthemes,
               gridExtra,
               grid,
               readxl,
               nls2,
               stats,
               ggplot2,
               pkgcond,
               rlist,
               pracma,
               tidyverse,
               splines,
               mgcv,
               furrr,
               fs,
               ComplexUpset,
               patchwork,
               caret,
               ggpubr,
               furrr,
               parallel,
               janitor,
               ggrepel,
               ggpubr,
               pROC,
               STRINGdb,
               UniProt.ws,
               'org.Hs.eg.db',
               zebrafish.db,
               VGAM,
               RColorBrewer,
               TPP,
               DEqMS,
               limma,
               VGAM,
               RColorBrewer,
               flow,
               tidylog,
               scam,
               BiocParallel,
               tidylog,
               clipr,
               report,
               equatiomatic,
               plyr,
               performance)
#install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("minpack.lm",
                       "rlist",
                       "data.table",
                       "knitr",
                       "ggthemes",
                       "AnnotationDbi",
                       "BiocGenerics",
                       "caret",
                       "ggpubr",
                       "ggrepel",
                       "STRINGdb",
                       "Uniprot.ws",
                       "org.Hs.Eg.db",
                       "zebrafish.db",
                       "TPP",
                       "DEqMS",
                       "equatiomatic"))

library(tidyverse)
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
read_cetsa <- function(protein_path,peptide_path,Prot_Pattern,Peptide=FALSE,Frac=TRUE,CFS=TRUE,solvent="DMSO",CARRIER=TRUE,rank=TRUE,sub=NA,temperatures=temps,baseline="min",NORM="QUANTILE",keep_shared_proteins=FALSE){
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
  
  read_xl<-function(x) {
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
    if (any(names(df2)=="Master.Protein.Accessions")&!any(names(df2)=="Accession")){
      df2<-df2 %>% dplyr::mutate("Accession"="Master.Protein.Accessions")
    }else if(any(stringr::str_detect(names(df2),"Master"))&!any(names(df2)=="Accession")){
      df2<-df2 %>% dplyr::mutate("Accession"=names(df2)[stringr::str_detect(names(df2),"Master")])
    }else if(any(stringr::str_detect(names(df2),"Accession"))){
      df2<-df2 %>% dplyr::mutate("Accession"=names(df2)[stringr::str_detect(names(df2),"Accession")][1])
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
                                   treatment=ifelse(stringr::str_detect(Spectrum.File,solvent),"vehicle","treated"),
                                   CC=ifelse(stringr::str_detect(Spectrum.File,solvent),0,1))
        if(any(stringr::str_detect(df2$File.ID,"."))){
          df2<-df2 %>% dplyr::mutate(Fraction=stringr::str_remove(.$File.ID,"[[:upper:]]+[[:digit:]]+."),
                                     sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
        }else{#if theres no fractionation
          df2<-df2 %>% dplyr::mutate(sample_id = stringr::str_extract(.$File.ID,"F[[:digit:]]+"))
        }
        
      }else{#if this experiment does not have Carrier FAIMS and PhiSDM involved shorten he file names
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
  
  
  read_PD<-function(x){#rename data for peptide groups
    
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
      dplyr::mutate(treatment=ifelse(stringr::str_detect(.$id,solvent),"vehicle","treated"),
                    CC=ifelse(stringr::str_detect(.$id,solvent),0,1),
                    sample_id = stringr::str_extract(.$id,"F[[:digit:]]+"),
                    temp_ref = unlist(stringr::str_extract(.$id,"[:digit:][:digit:][:digit:][N|C]|126|131")),
                    value = as.numeric(value))
    
    return(df2)
  }
  
  #read_PSMs and proteins
  if (.Platform$OS.type=="windows"){
    Proteins<-parallel::mclapply(f,read_xl)
    PSMs<-parallel::mclapply(h,read_xl)
  }else{
    Proteins<-parallel::mclapply(f,read_xl,mc.cores = future::availableCores())
    PSMs<-parallel::mclapply(h,read_xl,mc.cores = future::availableCores())
  }
  
  if(any(stringr::str_detect(Peptide,c("PG","PSMs")))){
    
    if(!is.na(sub)&any(class(sub)=="numeric")){#if there's only a subset of PSMs
      n<-as.numeric(sub)
      PSMs<-dplyr::bind_rows(PSMs) %>%
        dplyr::group_by(Protein.Accessions) %>%
        group_split()
      #subset a set number of PSMs
      PSMs<-PSMs[1:n]
    }else if(is.na(sub)){
      PSMs<-PSMs
    }else if(sub=="Filter"&keep_shared_proteins==FALSE){
      PSMs<-furrr::future_map(PSMs,function(x)filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,80,filter_rank=FALSE,keep_shared_proteins=FALSE,CFS=CFS,Frac=Frac))
      
    }else{
      PSMs<-furrr::future_map(PSMs,function(x)filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,80,filter_rank=FALSE,keep_shared_proteins=TRUE,CFS=CFS,Frac=Frac))
    }
    
    if (Peptide=="PSMs"){ #if this is a PSM file
      #split PSM subset into lists by way of sample_id name
      PSMs<-dplyr::bind_rows(PSMs) %>%
        dplyr::group_by(sample_name) %>% 
        dplyr::group_split()
      if(any(stringr::str_detect(names(PSMs[[1]]),"File_ID"))){
        PSMs<-dplyr::bind_rows(PSMs) %>%
          dplyr::rename("File.ID"="File_ID")%>%
          dplyr::group_by(sample_name) %>% 
          dplyr::group_split()
      }
      if(any(stringr::str_detect(names(PSMs[[1]]),"dataset"))){
        PSMs<-dplyr::bind_rows(PSMs) %>%
          dplyr::rename("treatment"="dataset")%>%
          dplyr::group_by(sample_name) %>% 
          dplyr::group_split()
      }
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
            message("Scaling factor calculated per protein,annotated sequence and treatment")
          }
        )
        )
        x<-dplyr::bind_rows(x)
        return(x)
      }
      )
      #if any data has the bioreplicate number, truncate the sample_id name contents
      if(!isTRUE(CFS)&any(stringr::str_detect(unique(dplyr::bind_rows(df2)$sample_name),"_[[:digit:]]+"))){
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
      
      #make sure the data is ready to be processed by sample_id name
      if(any(names(df2[[1]])=="sample_name")){
        df2<-dplyr::bind_rows(df2) %>% 
          dplyr::group_by(sample_name) %>%
          dplyr::group_split(.)
      }
      if(any(stringr::str_detect(names(df2[[1]]),"Accession"))){
        df2<-dplyr::bind_rows(df2) %>% dplyr::rename("uniqueID"="Accession")
      }
      if(isTRUE(Frac)){
        PSMs<-dplyr::bind_rows(df2) %>% dplyr::select(uniqueID,Spectrum.File,sample_name,sample_id,Fraction)
      }else{#if this isnt fractionated
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
      df2<-df2 %>% 
        tidylog::pivot_longer(cols=colnames(df2)[stringr::str_detect(colnames(df2),"[:digit:][:digit:][:digit:][N|C]|126|131")],
                              names_to = "id",
                              values_to ="value") %>% 
        dplyr::mutate(treatment= ifelse(stringr::str_detect(.$Spectrum.File,solvent),"vehicle","treated"),
                      CC= ifelse(stringr::str_detect(.$Spectrum.File,solvent),0,1),
                      sample_id = stringr::str_extract(.$id,"[[:upper:]]+[[:digit:]]+"),
                      temp_ref = stringr::str_extract(.$id,"[:digit:][:digit:][:digit:][N|C]|126|131"),
                      value = as.numeric(value),
                      sample_name = ifelse(length(unique(.$sample_id))==4,
                                           unique(stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))[2],
                                           stringr::str_extract(stringr::str_to_lower(.$Spectrum.File),"[[:lower:]]+_[[:digit:]]+"))) 
      if(!any(stringr::str_detect(names(dplyr::bind_rows(df2)),"temperature"))){
        df2<-df2 %>% dplyr::right_join(temperatures)%>% 
          dplyr::group_by(sample_name) %>% 
          dplyr::group_split()
      }
      return(df2)
    }else if(Peptide=="PG"){
      if(length(list.files(peptide_path,pattern="PeptideGroups.xlsx"))>0){
        f<-list.files(peptide_path,pattern="PeptideGroups.xlsx")
      }
      #first get Peptide group file read
      if (.Platform$OS.type=="windows"){
        df.raw<-parallel::mclapply(f, read_xl)
      }else{
        df.raw<-parallel::mclapply(f, read_xl, mc.cores = future::availableCores())
      }
      #get row number for peptide groups in case of batch files present
      df.raw1<-purrr::map2(df.raw,seq(df.raw),function(x,y) x %>% dplyr::mutate(n=y))
      #names<-purrr::map2(df.raw,f,function(x,y) ifelse(length(names(x))<45,warning(paste0("Please check the columns on file names",y)),print("All files have all necessary columns")))
      df.raw<-dplyr::bind_rows(df.raw1)
      
      PG<-read_PD(df.raw)#read in peptide groups
      
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
    missing_label<-function(x) {
      x<-x %>% dplyr::group_by(Accession,sample_id,sample_name,treatment) %>%
        distinct(.) %>% dplyr::mutate(missing=is.na(value),
                                      missing_pct=sum(is.na(value))/nrow(.))
      return(x)
    }
    if (.Platform$OS.type=="windows"){
      Proteins<-parallel::mclapply(Proteins,missing_label)
    }else{
      Proteins<-parallel::mclapply(Proteins,missing_label,mc.cores=future::availableCores())
    }
    if(isTRUE(Frac)){#if this is a fractionated treatment
      
      if(baseline=="min"){
        rank_label<-function(x){#if baseline is min and this is a fractionated treatment
          x<-x%>% dplyr::filter(temperature==min(temperature,na.rm=TRUE)) %>%
            dplyr::mutate(rank=dplyr::ntile(value,3))%>%
            dplyr::select(sample_id,sample_name,Accession,rank,id,Fraction,Spectrum.File)%>%
            dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>% dplyr::select(-id) %>% distinct(.)
          return(x)
        }
      }else{#if baseline is max and this is a fractionated treatment
        rank_label<-function(x){
          x<-x%>% dplyr::filter(temp_ref==unique(x$temperature)[length(unique(x$temperature))]) %>%
            dplyr::mutate(rank=dplyr::ntile(value,3))%>%
            dplyr::select(sample_id,sample_name,Accession,rank,id,Fraction,Spectrum.File)%>%
            dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>% dplyr::select(-id) %>% distinct(.)
          return(x)
        }
      }
    }else{#if this is unfractionated
      if(baseline=="min"){#if baseline is min and this is an unfractionated treatment
        rank_label<-function(x){
          x<-x%>% dplyr::filter(temperature==min(temperature,na.rm=TRUE)) %>%
            dplyr::mutate(rank=dplyr::ntile(value,3))%>%
            dplyr::select(sample_id,sample_name,Accession,rank,id,Spectrum.File)%>%
            dplyr::filter(!is.na(rank),!is.na(id)) %>% dplyr::select(-id) %>% distinct(.)
          return(x)
        }
      }else{#if baseline is max and this is an unfractionated treatment
        rank_label<-function(x){
          x<-x%>% dplyr::filter(temp_ref==unique(x$temperature)[length(unique(x$temperature))]) %>%
            dplyr::mutate(rank=dplyr::ntile(value,3))%>%
            dplyr::select(sample_id,sample_name,Accession,rank,id,Spectrum.File)%>%
            dplyr::filter(!is.na(rank),!is.na(id),rank==min(.$rank,na.rm=TRUE)) %>%
            dplyr::select(-id) %>% distinct(.)
          return(x)
        }
      }
    }
    if(isTRUE(rank)){
      if (.Platform$OS.type=="windows"){
        df_raw_D_R<-parallel::mclapply(Proteins,rank_label)
      }else{
        df_raw_D_R<-parallel::mclapply(Proteins,rank_label,mc.cores=future::availableCores())
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

#'Choose PSMs
#'
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
  df2<-dplyr::bind_rows(df2)
  #group data by sample_id name
  df2<-df2 %>%
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
  # df2<-dplyr::bind_rows(df2) %>%
  #   dplyr::group_split(sample_name)
  # ggplot(df2[[1]],mapping=aes(x=temp_ref,y=value))+geom_boxplot()+ylim(0,5)
  # #
  
  if(any(stringr::str_detect(names(df2[[1]]),"Accession"))){
    df2<-purrr::map(df2,function(x) x %>% dplyr::rename("uniqueID"="Accession"))
  }
  return(df2)
}
#' Clean data
#'
#' Clean CETSA data
#'
#' @param df.  Data frame returned by read_cetsa
#' @param temperatures.  Data frame of temperatures related to TMT tags - columns =
#'     temp_ref, temperature
#' @param samples.  Optional data frame of sample_id names - columns = sample_id, sample_name
#' @param separator.  Character used to separate parts of sample_id name.  Sample name is constructed as
#'     sampleroot_A_B where sampleroot is the name of the sample_id, A is the number of the biological
#'     replicate and B is the number of the technical replicate.
#'
#' @return a data frame of clean data
#'
#' @import dplyr
#'
#' @export
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
    error = function(cond){return(NA)}),mc.cores=future::availableCores())
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
#' Normalize by PSM
#' 

normalize_psm <- function(df, temperatures,filters=FALSE,CARRIER=TRUE) {
  if(isTRUE(CARRIER)){
    df<-df %>% dplyr::filter(!temperature=="68")
    temperatures<-temperatures[temperatures$temperature<68,]
  }
  if(any(names(df)=="value")){
    df<-df %>% dplyr::rename("I"="value")
  }
  
  if(any(names(df)=="uniqueID")){
    df<-df %>% dplyr::rename("Accession"="uniqueID","value"="I")
  }
  
  df$Accession<-as.factor(df$Accession)
  df<-data.frame(df)
  #rank by intensity at the lowest temperature channel
  df_filt<-df %>% dplyr::filter(temp_ref==temp_ref[temperature=min(df$temperature,na.rm=TRUE)]) %>% dplyr::mutate(rank=dplyr::ntile(dplyr::desc(I),3)) %>%
    dplyr::select(Accession,Annotated_Sequence,I,sample_id,treatment,sample_name,rank)
  df_filt$rank<-as.factor(df_filt$rank)
  
  #select top3 top5 and top10 peptides according to rank
  df_filt3 <- df_filt %>%dplyr::filter(!is.na(I)) %>%  
    dplyr::arrange(dplyr::desc(I)) %>% dplyr::group_by(rank)
  df_filt3<-df_filt3[sample(c(1:6000),6000,replace=TRUE),]
  
  df_filt5 <- df_filt %>%dplyr::filter(!is.na(I)) %>%  
    dplyr::arrange(dplyr::desc(I)) %>% dplyr::group_by(rank)
  df_filt5<-df_filt5[sample(c(1:12000),12000,replace=TRUE),]
  
  df_filt10 <- df_filt %>%dplyr::filter(!is.na(I)) %>%  
    dplyr::arrange(dplyr::desc(I)) %>% dplyr::group_by(rank)
  df_filt10<-df_filt10[sample(c(1:18000),18000,replace=TRUE),]
  
  df_filt3<-df_filt3 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
  df_filt5<-df_filt5 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
  df_filt10<-df_filt10 %>% dplyr::ungroup(.) %>%  dplyr::select(-I,-rank)
  #preserve original treatment
  df1<-df
  df<-df %>% group_by(sample_id)
  
  name<-dplyr::intersect(names(df1),names(df_filt3))
  #Only keep curves with the topN values
  df3<-df1 %>% dplyr::right_join(df_filt3,by=name)
  df5<-df1 %>% dplyr::right_join(df_filt5,by=name)
  df10<-df1 %>% dplyr::right_join(df_filt10,by=name)
  #remove missing values 
  df3<-df3[!is.na(df3$I),]
  df5<-df5[!is.na(df5$I),]
  df10<-df10[!is.na(df10$I),]
  
  #Calculate fold changes acros temperatures 7, 9 and 10
  
  
  df3<-dplyr::bind_rows(df3) %>%
    dplyr::group_split(Accession,sample_id)
  df5<-dplyr::bind_rows(df5) %>%
    dplyr::group_split(Accession,sample_id)
  df10<-dplyr::bind_rows(df10) %>%
    dplyr::group_split(Accession,sample_id)
  if(order(df3[[1]]$temperature)[1]==length(df3[[1]]$temperature)){
    #if the temperatures are reversed, relabel 7,9th and 10th temperature channels
    df.jointP3 <- suppressWarnings(purrr::map(df3,function(x) x %>%  
                                                dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE))+2,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)),
                                                              T9 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE))+1,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)),
                                                              T10 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)))))
    
    df.jointP5 <- suppressWarnings(purrr::map(df5,function(x) x %>%  
                                                dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE))+2,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)),
                                                              T9 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE))+1,"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)),
                                                              T10 = try(mean(x[which(x$temperature==max(x$temperature,na.rm=TRUE)),"I"]$I/x[which(x$temperature==min(x$temperature,na.rm=TRUE)),"I"]$I)))))
    df.jointP10 <- suppressWarnings(purrr::map(df10,function(x) x %>%  
                                                 dplyr::mutate(.,T7 = try(mean(x[which(x$temperature==max(x$temperature))+2,"I"]$I/x[which(x$temperature==min(x$temperature)),"I"]$I)),
                                                               T9 = try(mean(x[which(x$temperature==max(x$temperature))+1,"I"]$I/x[which(x$temperature==min(x$temperature)),"I"]$I)),
                                                               T10 = try(mean(x[which(x$temperature==max(x$temperature)),"I"]$I/x[which(x$temperature==min(x$temperature)),"I"]$I)))))
  }else{
    df.jointP3 <- suppressWarnings(purrr::map(df3,function(x) x %>%  
                                                dplyr::mutate(.,T7 = try(mean(x[x$temperature %in% temperatures$temperature[7],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)),
                                                              T9 = try(mean(x[x$temperature %in% temperatures$temperature[9],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)),
                                                              T10 = try(mean(x[x$temperature %in% temperatures$temperature[10],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)))))
    
    df.jointP5 <- suppressWarnings(purrr::map(df5,function(x) x %>%  
                                                dplyr::mutate(.,T7 = try(mean(x[x$temperature %in% temperatures$temperature[7],]$I/x[x$temperature %in% temperatures$temperature[1],]$I,na.rm=TRUE)),
                                                              T9 = try(mean(x[x$temperature %in% temperatures$temperature[9],]$I/x[x$temperature %in% temperatures$temperature[1],]$I,na.rm=TRUE)),
                                                              T10 = try(mean(x[x$temperature %in% temperatures$temperature[10],]$I/x[x$temperature %in% temperatures$temperature[1],]$I,na.rm=TRUE)))))
    
    
    df.jointP10 <- suppressWarnings(purrr::map(df10,function(x) x %>%  
                                                 dplyr::mutate(.,T7 = try(mean(x[x$temperature %in% temperatures$temperature[7],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)),
                                                               T9 = try(mean(x[x$temperature %in% temperatures$temperature[9],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)),
                                                               T10 = try(mean(x[x$temperature %in% temperatures$temperature[10],"I"]$I/x[x$temperature %in% temperatures$temperature[1],"I"]$I,na.rm=TRUE)))))
    
    
  }
  
  
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
    df1 <- suppressWarnings(df1 %>% dplyr::group_split(Accession,sample_id) %>% 
                              purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                           dplyr::mutate(.,T7 = try(mean(x[x$temperature %in% temperatures[7],]$I/x[x$temperature %in% temperatures[1],]$I,na.rm=TRUE)),
                                                         T9 = try(mean(x[x$temperature %in% temperatures[9],]$I/x[x$temperature %in% temperatures[1],]$I,na.rm=TRUE)),
                                                         T10 = try(mean(x[x$temperature %in% temperatures[10],]$I/x[x$temperature %in% temperatures[1],]$I,na.rm=TRUE)))))
    
    
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
  l.bytype3 <- split.data.frame(df.jointP3, df.jointP3$sample_id)
  l.bytype5 <- split.data.frame(df.jointP5, df.jointP5$sample_id)
  l.bytype10 <- split.data.frame(df.jointP10, df.jointP10$sample_id)
  
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
  
  
  df.median3 <- df.jointP3 %>%
    dplyr::group_by(temperature,treatment) %>%
    dplyr::mutate(value = median(I,na.rm=TRUE))
  
  df.median5 <- df.jointP5 %>%
    dplyr::group_by(temperature,treatment) %>%
    dplyr::mutate(value = median(I,na.rm=TRUE))
  
  df.median10 <- df.jointP10 %>%
    dplyr::group_by(temperature,treatment) %>%
    dplyr::mutate(value = median(I,na.rm=TRUE))
  
  
  #fit curves to median 
  ## fit curves to the median data for each sample (F1 through FN)
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
    ),silent=TRUE)))
  df.fit3<-df.fit3 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
  DF.fit3<-df.fit3
  df.fit3<-df.fit3 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
  if(!isTRUE(nrow(df.fit3)==0)){
    warning("Initial subset of 3000 PSMs has failed to fit.  No correction factors will be processed.")
    df.fit3<-dplyr::bind_rows(DF.fit3)
  }else{
    df.fit3<-df.fit3%>% 
      dplyr::mutate(fitted_values3 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
      dplyr::select(sample_id,fitted_values3,temperature,treatment) %>% ungroup(.)
  }
  
  
  
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
    ),silent=TRUE)))
  df.fit5<-df.fit5 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
  df.fit5<-df.fit5 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
  df.fit5<-dplyr::bind_rows(df.fit5)
  df.fit5<- df.fit5%>% 
    dplyr::mutate(fitted_values5 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
    dplyr::select(sample_id,fitted_values5,temperature,treatment) %>% ungroup(.)
  
  df.fit10 <- df.median10 %>%
    dplyr::group_by(sample_id,treatment) %>% 
    dplyr::mutate(fit = list(try(nls(formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
                                     start = c(Pl=0, a = 550, b = 10),
                                     data = list(x=temperature,y=value),
                                     na.action = na.exclude,
                                     algorithm = "port",
                                     lower = c(0.0,1e-5,1e-5),
                                     upper = c(1.5,15000,300),
                                     control = nls.control(maxiter = 50)
    ),silent=TRUE)))
  df.fit10<-df.fit10 %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
  df.fit10<-df.fit10 %>% purrr::keep(function(x) class(x$fit[[1]])=='nls')
  df.fit10<-dplyr::bind_rows(df.fit10)
  
  df.fit10<- df.fit10%>% 
    dplyr::mutate(fitted_values10 = try(list(data.frame(fitted_values=predict(fit[[1]]))),silent=TRUE)) %>% 
    dplyr::select(sample_id,fitted_values10,temperature,treatment) %>% ungroup(.)
  ## calculate the fitted values
  d3<-df.fit3 %>% dplyr::group_split(sample_id,treatment,temperature)
  d5<-df.fit5 %>% dplyr::group_split(sample_id,treatment,temperature)
  d10<-df.fit10 %>% dplyr::group_split(sample_id,treatment,temperature)
  
  
  #unnest fitted values from list and name value column and keep fitted values and temps
  check3 <-purrr::map(d3,function(x) x [1,])
  check5 <-purrr::map(d5,function(x) x [1,])
  check10 <-purrr::map(d10,function(x) x [1,])
  
  check3<-purrr::map(check3,function(x) x %>% tidyr::unnest(cols=fitted_values3) %>% unique(.))
  check5<-purrr::map(check5,function(x) x %>% tidyr::unnest(cols=fitted_values5) %>% unique(.)) 
  check10<-purrr::map(check10,function(x) x %>% tidyr::unnest(cols=fitted_values10) %>% unique(.)) #%>% dplyr::mutate(temperature=temperatures))
  # 
  #bind_rows
  check3<-dplyr::bind_rows(check3)
  check5<-dplyr::bind_rows(check5)
  check10<-dplyr::bind_rows(check10)
  
  name3<-dplyr::intersect(names(check3),names(df1))
  
  test3<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check3,name3)
  test5<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check5,name3)
  test10<-df1 %>% dplyr::group_by(temperature) %>% dplyr::right_join(check10,name3)
  
  ## calculate ratios between the fitted curves and the median values
  df.out3 <- test3 %>%
    dplyr::mutate(correction3 = try(fitted_values/I,silent=TRUE)) %>%
    dplyr::select('sample_id','temperature','I','fitted_values','correction3')
  #df.out3<-df.out3 %>% dplyr::select(-fitted_values,-I,-sample)
  
  df.out5 <- test5 %>%
    dplyr::mutate(correction5 = try(fitted_values / I,silent=TRUE)) %>%
    dplyr::select('sample_id','temperature','I','fitted_values','correction5')
  #df.out5<-df.out5 %>% dplyr::select(-fitted_values,-I,-sample)
  
  df.out10 <- test10 %>%
    dplyr::mutate(correction10 = try(fitted_values /I,silent=TRUE)) %>%
    dplyr::select('sample_id','temperature','I','fitted_values','correction10')
  # df.out10<-df.out10 %>% dplyr::select(-fitted_values,-I,-sample)
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
  
  df1<-dplyr::bind_rows(df1)
  df1<-df1 %>% dplyr::group_split(temperature)
  #apply correction factor by temperature to original data
  df3<-purrr::map2(df1,df.out3,function(x,y)tryCatch(x %>% dplyr::mutate(correction3=y$correction3[1])))
  df5<-purrr::map2(df3,df.out5,function(x,y)tryCatch(x %>% dplyr::mutate(correction5=y$correction5[1])))
  df<-purrr::map2(df5,df.out10,function(x,y)tryCatch(x %>% dplyr::mutate(correction10=y$correction10[1])))
  #bind_rows
  df<-dplyr::bind_rows(df)
  
  # df3 <- df%>% 
  #   dplyr::mutate(norm_value3= ifelse(is.na(correction3),I,I* correction3)) %>% dplyr::ungroup(.)
  # df5 <- df3%>% 
  #   dplyr::mutate(norm_value5= ifelse(is.na(correction5),I,I* correction5)) %>% dplyr::ungroup(.)
  # df<- df5%>% 
  #   dplyr::mutate(norm_value10= ifelse(is.na(correction10),I,I* correction10)) %>% dplyr::ungroup(.)
  # 
  # df<-df %>% dplyr::mutate(norm_value3=ifelse(is.na(norm_value3),NA,norm_value3),
  #                          norm_value5=ifelse(is.na(norm_value5),NA,norm_value5),
  #                          norm_value10=ifelse(is.na(norm_value10),NA,norm_value10)
  # )
  
  df <- df %>% 
    dplyr::rename("uniqueID"="Accession", "C"="temperature","C3"="correction3","C5"="correction5","C10"="correction10")
  df<-df %>% dplyr::ungroup(.) %>% dplyr::select(-I)
  
  if(isTRUE(filters)){
    if (any(names(df)==T7)){
      df<-df %>% dplyr::select(-T10,-T7,-T9)
    }
  }
  return(df)
  
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
    dplyr::select(Accession, sample_id, params) %>%
    spread(sample_id, params)
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
    fine_start <- expand.grid(p=c(0,0.1),k=seq(0,100),m=seq(5,45,by=10))
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


cetsa_fit2 <- function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    #c(Pl=0, a = 550, b = 10
    fine_start <- expand.grid(p=c(0,0.2),k=seq(0,1000),m=seq(5,100,by=5))
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
  fgh2 <- deriv(value ~ (1 - p)/(1 + exp(-k*(1/t - 1/m))) + p, c("p","k","m"), function(p,k,m){} ) 
  #create a grid
  t<-seq.int(37,67,0.05)
  beta2.est <- coef(new_start)
  f.new <- fgh2(beta2.est[1],beta2.est[2],beta2.est[3])
  
  g.new <- base::attr(f.new,"gradient")
  g.new1 <- t(as.vector(attr(f.new,"gradient")))
  V.beta2 <- vcov(new_start)
  GS<-rowSums((g.new%*%V.beta2)*g.new)
  #95% CI
  alpha <- 0.05
  deltaf <- sqrt(GS)*qt(1-alpha/2,summary(new_start)$df[2])
  df.delta <- data.frame(temperature=t, value=f.new, lwr.conf=f.new-deltaf, upr.conf=f.new+deltaf)
  sigma2.est <- summary(new_start)$sigma
  deltay <- sqrt(GS + sigma2.est^2)*qt(1-alpha/2,summary(new_start)$df[2])
  df.delta[c("lwr.pred","upr.pred")] <- cbind(f.new - deltay,f.new + deltay)
  #pl<-ggplot(d)+geom_point(mapping=aes(x=temperature,y=value))
  # 
  # pl + geom_ribbon(data=df.delta, aes(x=t, ymin=lwr.pred, ymax=upr.pred), alpha=0.1, fill="blue") +
  #   geom_ribbon(data=df.delta, aes(x=t, ymin=lwr.conf, ymax=upr.conf), alpha=0.2, fill="#339900") +
  #   geom_line(data=df.delta, aes(x=t, y=f.new), colour="#339900", size=1)
  return(list(result,df.delta))
}


#' return Tm from CETSA curve
#'
#' return Tm from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return Tm.
#'
Tmv <- function(f) {
  if (length(f) == 1) return(NA)
  pars<-f$fit[[1]]$m$getPars()
  Tm<-pars['a']/(pars['b'] - log(1-pars['Pl'])/(1/2 - pars['Pl']-1))
  return(round(Tm,1)[['a']]) 
}
#calculateTm for splines
calcTm<-function(x) {tryCatch(expr={
  
  y<-with(x, stats::approx(x$fitted_values[[1]][1:length(x$C)],x$C,
                           xout=min(x$fitted_values[[1]],
                                    na.rm=TRUE)+
                             (0.5*(max(x$fitted_values[[1]],
                                       na.rm=TRUE)-
                                     min(x$fitted_values[[1]], na.rm=TRUE))))$y)
},error=function(cond){
  warning("Tm value not able to be calculated")
  return(NA)
})
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
  
  df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                            "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                            "I_Interference"="Isolation_Interference_[%]")
  
  
  #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
  df_1<-dplyr::bind_rows(df_)
  df_1$uniqueID<-as.factor(df_1$uniqueID)
  df_1$treatment<-as.factor(df_1$treatment)
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1<-df_1%>% 
    dplyr::group_split(uniqueID,sample_id,Annotated_Sequence)
  df_1<-purrr::map(df_1,function(x) x %>% dplyr::mutate(missing_pct=(100*sum(is.na(x$I))/length(x$I))) %>% head(.,1)) 
  df_1<-dplyr::bind_rows(df_1)
  
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
      dplyr::select(-uniqueID,-Spectrum.File,-Annotated_Sequence,-IonInjTime,-MissedCleavages,-S_N,-sample_id,-treatment,-C,-CC,-DeltaM,-PEP,-Modifications,-I_Interference,-missing_pct,-XCorr,-I,-Protein_value))
  
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
upMV <- function(df_,condition,plot_multiple=FALSE,PSMs=FALSE,Frac=TRUE){
  
  if(isTRUE(plot_multiple)){
    df_<-dplyr::bind_rows(df_)
    if(isTRUE(PSMs)){
      #df_<- df_ %>% dplyr::right_join(df.samples,by="sample_id")
      df_ <-df_ %>% dplyr::mutate(CC=ifelse(stringr::str_detect(treatment,"vehicle")==TRUE,0,1))#concentration values are defined in uM
      
      if(any(stringr::str_detect(names(df_),"Accession"))){
        df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                  "IonInjTime"="Ion_Inject_Time_[ms]","I_Interference"="Isolation_Interference_[%]")
      }else{
        df_<-df_%>% dplyr::rename("I"="value","C"="temperature","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                  "IonInjTime"="Ion_Inject_Time_[ms]","I_Interference"="Isolation_Interference_[%]")
      }
      
      #saveRDS(df_,"df_raw_PSMs_Cliff.rds")
      df_1<-dplyr::bind_rows(df_)
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$treatment<-as.factor(df_1$treatment)
      df_1$sample_name<-as.factor(df_1$sample_name)
      if(isTRUE(Frac)){
        df_1<-dplyr::bind_rows(df_1)%>% 
          dplyr::group_split(uniqueID,sample_id,Annotated_Sequence,Fraction,treatment,sample_name)
      }else{
        df_1<-dplyr::bind_rows(df_1)%>% 
          dplyr::group_split(uniqueID,sample_id,Annotated_Sequence,treatment,sample_name)
      }
      
      df_1<-dplyr::bind_rows(df_1)
      if(any(is.na(df_1$rank))){
        if(any(stringr::str_detect(names(df_1),"Fraction"))){
          rank<-df_1 %>% dplyr::select(-rank) %>% dplyr::filter(C==min(.$C,na.rm=TRUE)) %>% dplyr::mutate(rank=dplyr::ntile(dplyr::desc(.$I),3),
                                                                                                          SNrank=dplyr::ntile(dplyr::desc(.$S_N),3)) %>% 
            dplyr::select("uniqueID","treatment","sample_name","sample_id","rank","Fraction","Annotated_Sequence","missing_pct","SNrank","Charge","S_N",-temp_ref)
          maxSNrank<-min(rank[rank$SNrank==1,'S_N'],na.rm=TRUE)
          minSNrank<-max(rank[rank$SNrank==3,'S_N'],na.rm=TRUE)
          low<-paste0("low"," < ",as.character(minSNrank))
          medium<-"medium"
          high<-paste0("high"," > ",as.character(maxSNrank))
          rank<-rank %>% dplyr::mutate(rank=ifelse(rank==1,"high",ifelse(rank==2,"medium","low")),
                                       SNrank=ifelse(SNrank==3,"low",ifelse(SNrank==2,"medium","high")))
          rank$SNrank<-factor(rank$SNrank)
          rank$rank<-factor(rank$rank)
          df_1<-dplyr::bind_rows(rank)
          
        }else{
          rank<-df_1 %>% dplyr::select(-rank) %>% dplyr::filter(C==min(.$C,na.rm=TRUE)) %>% dplyr::mutate(rank=dplyr::ntile(dplyr::desc(.$I),3),
                                                                                                          SNrank=dplyr::ntile(dplyr::desc(.$S_N),3)) %>% 
            dplyr::select("uniqueID","treatment","sample_name","sample_id","rank","Annotated_Sequence","missing_pct","SNrank","Charge","S_N",-temp_ref)
          maxSNrank<-min(rank[rank$SNrank==1,'S_N'],na.rm=TRUE)
          minSNrank<-max(rank[rank$SNrank==3,'S_N'],na.rm=TRUE)
          low<-paste0("low"," < ",as.character(minSNrank))
          medium<-"medium"
          high<-paste0("high"," > ",as.character(maxSNrank))
          rank<-rank %>% dplyr::mutate(rank=ifelse(rank==1,"high",ifelse(rank==2,"medium","low")),
                                       SNrank=ifelse(SNrank==3,"low",ifelse(SNrank==2,"medium","high")))
          rank$SNrank<-factor(rank$SNrank)
          rank$rank<-factor(rank$rank)
          df_1<-dplyr::bind_rows(rank)
        }
        #name<-dplyr::intersect(names(df_1),names(rank))
        #df_1<-df_1 %>% dplyr::right_join(rank,by=name)
        
      }
      #df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
      df_1$missing_pct<-100-round(df_1$missing_pct,1)
      df_1<-dplyr::bind_rows(df_1)%>% 
        dplyr::group_split(sample_name)
      #separate data by charge states
      Df_1<-purrr::map(df_1,function(x)x %>% dplyr::select(uniqueID,Annotated_Sequence,treatment,sample_id,Charge,rank,sample_name,SNrank) %>% distinct(.))
      
      hi<-dplyr::bind_rows(Df_1) %>% pivot_wider(names_from="Charge",values_from="Charge")
      #df_1<-purrr::map(hi,function(x)x %>% dplyr::rename_if(is.numeric))
      df_1<-hi %>% rename_at(names(.)[which(sapply(hi,class)=="numeric")],~paste0("Charge: ", .))
      df_1<-dplyr::bind_rows(df_1) %>% dplyr::group_split(sample_name)
      if(any(names(df_1[[1]])=="SNrank")){
        df_1<-purrr::map(df_1,function(x)x[!is.na(x$SNrank),])
      }
      
      rating_scale = scale_fill_manual(name="Ranked S/N",
                                       values=c("high" = '#e6550d', "medium" ='#fdae6b', "low"  = '#fee6ce'))
      
      show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
      
      check<-list()
      check<- purrr::map(df_1,function(x)upset(x,names(x)[which(stringr::str_detect(names(x),"Charge"))],
                                               min_degree=1,
                                               set_sizes=FALSE,
                                               guides='collect',
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
                                               
      )+ggtitle(stringr::str_replace(x$sample_name[1],"S",paste0("\u03A6"))))
      check1<-upset(df_1[[1]],names(df_1[[1]])[which(stringr::str_detect(names(df_1[[1]]),"Charge"))],
                    min_degree=1,
                    set_sizes=FALSE,
                    guides='collect',
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
      
    }else{#If this is a protein file
      df_1<-dplyr::bind_rows(df_)
      df_1<-df_1%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature")
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1$treatment<-as.factor(df_1$treatment)
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
          dplyr::select(-C,-I,-missing,-treatment,-sample_id))
      
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
    check<- purrr::map(df_1,function(x)upset(x,names(x)[!names(x)%in%c("rank","sample_name")],
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
    
    P<-ggarrange(plotlist=check,ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
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
    df_$treatment<-as.factor(df_$treatment)
    
    d2<-dplyr::bind_rows(df_)%>% dplyr::filter(treatment=="vehicle") 
    d2$uniqueID<-as.factor(d2$uniqueID)
    d2$sample_name<-as.factor(d2$sample_name)
    d2$C<-as.factor(d2$C)
    d2$rank<-as.factor(d2$rank)
    
    
    d1<-dplyr::bind_rows(df_)%>% dplyr::filter(treatment=="treated") 
    d1$uniqueID<-as.factor(d1$uniqueID)
    d1$sample_name<-as.factor(d1$sample_name)
    d1$C<-as.factor(d1$C)
    d1$rank<-as.factor(d1$rank)
    
    
    d3<-rbind(d1,d2)
    
    d3<-tidyr::pivot_wider(d3,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d1<-pivot_wider(d1,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d2<-pivot_wider(d2,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    
    d1<-d1%>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    d2<-d2 %>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    d3<-d3 %>% dplyr::select(-uniqueID,-C,-I,-missing,-treatment,-sample_id,-sample_name)
    
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


#Spline
# spstat<-function(DF,df,df1,Ftest=TRUE,show_results=TRUE,filters=TRUE,scaled_dof=FALSE){
#   if(any(class(df)=="list")){
#     df<-df %>% purrr::keep(function(x) is.data.frame(x))
#     df1<-df1 %>% purrr::keep(function(x) is.data.frame(x))
#     DF<-DF %>% purrr::keep(function(x) is.data.frame(x))
#     
#     df<-dplyr::bind_rows(df)
#     df1<-dplyr::bind_rows(df1)
#     DF<-dplyr::bind_rows(DF)
#   }
#   #plot spline results
#   
#   df1$I<-as.numeric(as.vector(df1$I))
#   df$I<-as.numeric(as.vector(df$I))
#   DF$I<-as.numeric(as.vector(DF$I))
#   #mutate to get CV values
#   DF<-DF %>% dplyr::group_split(C,uniqueID) 
#   DF<- purrr::map(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
#   df<-df %>% dplyr::group_split(C,uniqueID,treatment) 
#   df<- purrr::map(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
#   df1<-df1 %>% dplyr::group_split(C,uniqueID,treatment) 
#   df1<- purrr::map(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
#   
#   #convert to data frame
#   
#   df<-dplyr::bind_rows(df)
#   df1<-dplyr::bind_rows(df1)
#   DF<-dplyr::bind_rows(DF)
#   
#   #aggregate column 
#   #switch from factor to numeric
#   #convert factor to numeric columns
#   df1$C<-as.numeric(as.vector(df1$C))
#   df$C<-as.numeric(as.vector(df$C))
#   DF$C<-as.numeric(as.vector(DF$C))
#   
#   df1$I<-as.numeric(as.vector(df1$I))
#   df$I<-as.numeric(as.vector(df$I))
#   DF$I<-as.numeric(as.vector(DF$I))
#   
#   #remove NA vaues
#   df <- df[,which(unlist(lapply(df, function(x)!all(is.na(x))))),with=F]
#   df1 <- df1[,which(unlist(lapply(df1, function(x)!all(is.na(x))))),with=F]
#   DF <- DF[,which(unlist(lapply(DF, function(x)!all(is.na(x))))),with=F]
#   #convert back to list
#   DF<-DF %>% dplyr::group_split(uniqueID)
#   df<-df %>% dplyr::group_split(uniqueID)
#   df1<-df1 %>% dplyr::group_split(uniqueID)
#   
#   
#   #alternative spline fit method : Generalized Additive Models
#   #fit penalized splines
#   # if(any(names(df[[1]])=="time_point")){
#   #   df<-purrr::map(df,function(x) dplyr::bind_rows(x) %>% dplyr::group_split(uniqueID,time_point))
#   # }
#   m <- purrr::map(df,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(x$I ~ s(x$C,k=5), data = x , method = "REML")))))
#   m<-m %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
#   #check significance and refit data with more k 
#   m<-purrr::map(m,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
#                                                  sum = list(summary(.$M1[[1]])),
#                                                  Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
#                                                  rss=deviance(.$M1[[1]]),
#                                                  CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
#                                                  AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
#                                                  rsq=summary(x$M1[[1]])$r.sq,
#                                                  n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
#                                                  sample_name=.$sample_name))
#   m1 <- purrr::map(df1,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data = x , method = "REML")))))
#   m1<-m1 %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
#   #check significance and refit data with more k 
#   m1<-purrr::map(m1,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
#                                                    sum = list(summary(.$M1[[1]])),
#                                                    Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
#                                                    rss=deviance(.$M1[[1]]),
#                                                    CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
#                                                    AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
#                                                    rsq=summary(x$M1[[1]])$r.sq,
#                                                    n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
#                                                    sample_name=.$sample_name))
#   #m1<-lapply(df1,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
#   
#   
#   mn<- purrr::map(DF,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data =x, method = "REML",fit=TRUE)))))
#   mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
#   #check significance and refit data with more k 
#   mn<-purrr::map(mn,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
#                                                    sum = list(summary(.$M1[[1]])),
#                                                    Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=TRUE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
#                                                    rss=deviance(.$M1[[1]]),
#                                                    CV_pct=ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
#                                                    AUC = pracma::trapz(.$M1[[1]]$fit[(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fit-0.5)))+1)]),
#                                                    rsq=summary(.$M1[[1]])$r.sq,
#                                                    n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
#                                                    sample_name=.$sample_name,
#                                                    RSS0=deviance(.$M1[[1]])
#   ))
#   lm<-length(m)
#   lm1<-length(m1)
#   
#   if(lm<lm1){
#     CID<-unique(dplyr::bind_rows(m)$uniqueID)
#   }else{
#     CID<-unique(dplyr::bind_rows(m1)$uniqueID)
#   }
#   
#   m<-dplyr::bind_rows(m)
#   m1<-dplyr::bind_rows(m1)
#   mn<-dplyr::bind_rows(mn)
#   
#   
#   #filter
#   m <-m %>% dplyr::filter(uniqueID %in% CID)
#   m1<-m1%>% dplyr::filter(uniqueID %in% CID)
#   mn<-mn%>% dplyr::filter(uniqueID %in% CID)
#   
#   CID<-dplyr::intersect(unique(m$uniqueID),unique(m1$uniqueID))
#   
#   m <-m %>% dplyr::filter(uniqueID %in% CID)
#   m1<-m1%>% dplyr::filter(uniqueID %in% CID)
#   mn<-mn%>% dplyr::filter(uniqueID %in% CID)
#   #split 
#   m<-m %>% dplyr::group_split(uniqueID)
#   m1<-m1%>% dplyr::group_split(uniqueID)
#   mn<-mn%>% dplyr::group_split(uniqueID)
#   
#   #calculate RSS
#   m<-purrr::map2(m,m1,function(x,y) x %>% dplyr::mutate(RSSv=deviance(x$M1[[1]]),
#                                                         RSSt=deviance(y$M1[[1]]),
#                                                         RSS1=deviance(x$M1[[1]])+deviance(y$M1[[1]])
#   ))
#   m1<-purrr::map2(m,m1,function(x,y) y %>% dplyr::mutate(RSSv=deviance(x$M1[[1]]),
#                                                          RSSt=deviance(y$M1[[1]]),
#                                                          RSS1=deviance(x$M1[[1]])+deviance(y$M1[[1]])
#   ))
#   
#   check<-purrr::map2(m1,mn,function(x,y) x %>% dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
#   check<-dplyr::bind_rows(check)
#   
#   #4031 proteins have rssDiff> 0 
#   
#   #convert to df and split by uniqueID 
#   mean1<-dplyr::bind_rows(m)
#   mean1_1<-dplyr::bind_rows(m1)
#   mean3<-dplyr::bind_rows(mn)
#   #equal column names
#   mean1<-mean1  %>% dplyr::mutate(rss=RSS1)%>% dplyr::select(-RSSv,-RSSt,-RSS1)
#   mean1_1<-mean1_1  %>% dplyr::mutate(rss=RSS1)%>% dplyr::select(-RSSv,-RSSt,-RSS1)
#   mean3<-mean3  %>% dplyr::mutate(rss=RSS0)%>% dplyr::select(-RSS0)
#   
#   #filter out proteins with neg RSSdiff
#   if(isTRUE(filters)){
#     check<-check %>% dplyr::filter(rssDiff>0,missing_pct<=20,rsq>0.8)
#     check<-unique(check$uniqueID)
#     mean1<-mean1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
#     mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
#     mean3<-mean3 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
#     
#     check<-intersect(mean1$uniqueID,mean1_1$uniqueID)
#     
#     mean1<-mean1 %>% dplyr::filter(uniqueID %in% check)
#     mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check)
#     mean3<-mean3 %>% dplyr::filter(uniqueID %in% check)
#     
#   }
#   
#   
#   #split into lists by uniqueID
#   mean1<-mean1 %>% dplyr::group_split(uniqueID)
#   mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
#   mean3<-mean3 %>% dplyr::group_split(uniqueID)
#   
#   
#   #Cliff
#   results<-dplyr::bind_rows(mean1,mean1_1,mean3)
#   
#   if(isTRUE(Ftest)){
#     
#     #Calculate rss0 and rss1 null vs alt
#     rss0<-purrr::map(mean3,function(x)data.frame(RSS0=x$rss,
#                                                  se0=summary.gam(x$M1[[1]])$se[1],#standard error
#                                                  pN0=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
#                                                  rsq0=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
#                                                  np0=summary.gam(x$M1[[1]])$np[1],#number of parameters
#                                                  rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
#                                                  m=summary.gam(x$M1[[1]])$m[1]))
#     #of smooth terms in the model
#     
#     #mean1<-purrr::map2(mean1,mean1_1,function(x,y)x %>% purrr::keep(x$uniqueID[1] %in% y$uniqueID[1]))
#     #generate a grid of new values for C
#     
#     
#     #make predictions to generate point confidence intervals
#     
#     
#     rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS=x$rss,
#                                                             se=summary.gam(x$M1[[1]])$se[1],#standard error
#                                                             pN1=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
#                                                             rsq=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
#                                                             
#                                                             rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
#                                                             m=summary.gam(x$M1[[1]])$m[1],
#                                                             Tm = y$Tm[[1]]-x$Tm[[1]],
#                                                             se1=summary.gam(y$M1[[1]])$se[1],#standard error
#                                                             pN2=summary.gam(y$M1[[1]])$edf[1],#effective degrees of freedom treated
#                                                             rsq1=summary.gam(y$M1[[1]])$r.sq[1],#r-squared
#                                                             
#                                                             rdf1=summary.gam(y$M1[[1]])$residual.df[1],#residual degrees of freedom
#                                                             m1=summary.gam(y$M1[[1]])$m[1],
#                                                             pA=sum(summary.gam(x$M1[[1]])$edf[1],summary.gam(y$M1[[1]])$edf[1],na.rm=TRUE)))
#     
#     # rss1<-purrr::map2(rss1,mean1,function(x,y) x %>% dplyr::mutate(fit_v=list(ifelse(class(try(mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",
#     #                                                                                                              se.fit=TRUE,newdata.guaranteed = TRUE)))=='try-error',NA,
#     #                                                                                  mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",se.fit=TRUE,newdata.guaranteed = TRUE)))))
#     int<-dplyr::intersect(names(mean3[[1]]),names(rss0[[1]]))
#     int1<-dplyr::intersect(names(mean1[[1]]),names(rss1[[1]]))
#     
#     mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
#     mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
#     
#     #bind predicted values to main data frame
#     mean3<-purrr::map2(mean3,rss0,function(x,y)cbind(x,y))
#     mean1<-purrr::map2(mean1,rss1,function(x,y)cbind(x,y))
#     
#     
#     #mean1_1<-purrr::map2(mean1_1,rss1,function(x,y)cbind(x,y))
#     #params for null and alternative models
#     f0<-lapply(mean3,function(x) data.frame(np=length(x$M1[[1]]$fitted.values)))#number of measurements
#     f1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(nA=sum(nrow(x),nrow(y))))#number of measurements alternative
#     
#     int<-dplyr::intersect(names(mean3[[1]]),names(f0[[1]]))
#     int1<-dplyr::intersect(names(mean1[[1]]),names(f1[[1]]))
#     
#     
#     mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
#     mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
#     
#     #bind data frames
#     mean3<-purrr::map2(mean3,f0,function(x,y) cbind(x,y))
#     mean1<-purrr::map2(mean1,f1,function(x,y) cbind(x,y))
#     mean1_1<-purrr::map2(mean1_1,f1,function(x,y) cbind(x,y))
#     #calculate parameters 
#     pN<-purrr::map(mean3,function(x) x %>% dplyr::select(pN0,np))
#     pA<-purrr::map(mean1,function(x)x %>% dplyr::select(pA,nA))
#     #degrees of freedom before
#     d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN0))
#     d2<-purrr::map(pA,function(x)data.frame(d2=x$nA-x$pA))
#     
#     #delta RSS
#     rssDiff<-purrr::map2(rss0,rss1,function(x,y) data.frame(dRSS=x$RSS0-y$RSS,dTm=y$Tm[1],rss1=y$RSS[1]))
#     
#     d2<-lapply(d2,function(x) x$d2[1])
#     d1<-lapply(d1,function(x) x$d1[1])
#     
#     #RSS1 numerator
#     Fvals<-purrr::map2(rssDiff,d1,function(x,y) data.frame(fNum=x$dRSS/y[1]))
#     #Rss denominator
#     Fd<-purrr::map2(rss1,d2,function(x,y) data.frame(fDen=x$RSS/y[1]))
#     
#     Fvals<-purrr::map2(Fvals,Fd,function(x,y) data.frame(fStatistic=x$fNum/y$fDen))
#     Fvals<-purrr::map2(Fvals,d1,function(x,y) x %>% dplyr::mutate(df1=y[1]))
#     Fvals<-purrr::map2(Fvals,d2,function(x,y) x %>% dplyr::mutate(df2=y[1]))
#     
#     Fvals<-purrr::map2(Fvals,rssDiff,function(x,y) cbind(x,y))
#     Fvals<-purrr::map(Fvals,function(x) x %>% dplyr::mutate(Fvals=(x$dRSS/x$rss1)*(x$df2/x$df1)))
#     #add p-vals
#     Fvals<-purrr::map(Fvals,function(x) x %>% dplyr::mutate(pValue = 1 - pf(fStatistic, df1 = x$df1, df2 = x$df2),
#                                                             pAdj = p.adjust(pValue,"BH")))
#     Fvals<-purrr::map2(Fvals,mean1,function(x,y) x %>% dplyr::mutate(uniqueID=y$uniqueID[1]))
#     
#     
#     mean1<-purrr::map(mean1,function(x) data.frame(x)%>% dplyr::select(-temp_ref,-C,-I,-M1,-sum) )
#     
#     mean1<-purrr::map(mean1,function(x) x %>% distinct(.) )
#     names1<-dplyr::intersect(names(mean1),names(Fvals))
#     #convert results to list
#     testResults<-purrr::map2(mean1,Fvals,function(x,y)x%>% dplyr::right_join(y,by=names1))
#     testResults<-purrr::map(testResults,function(x) x[1,])
#     testResults<-dplyr::bind_rows(testResults)
#     mean1<-dplyr::bind_rows(mean1)
#     Unscaled<-ggplot(testResults)+
#       geom_density(aes(x=fStatistic),fill = "steelblue",alpha = 0.5) +
#       geom_line(aes(x=fStatistic,y= df(fStatistic,df1=df1,df2=df2)),color="darkred",size = 1.5) +
#       theme_bw() +
#       ggplot2::xlab("F-values")+
#       ggplot2::xlim(0,0.05)
#     # print(Unscaled)
#     #scale variables
#     M<-median(testResults$dRSS,na.rm=TRUE)
#     V<-mad(testResults$dRSS,na.rm=TRUE)
#     #alternative scaling factor sig0-sq
#     altScale<-0.5*V/M
#     #filter out negative delta rss
#     testResults<-testResults %>% dplyr::filter(dRSS>0)
#     #effective degrees of freedom
#     ed1<-MASS::fitdistr(x=testResults$dRSS, densfun = "chi-squared", start = list(df=2))[["estimate"]]
#     ed2<-MASS::fitdistr(x=testResults$rss, densfun = "chi-squared", start = list(df=2))[["estimate"]]
#     #scale data
#     testScaled <-testResults %>%
#       dplyr::mutate(rssDiff = .$dRSS/altScale,
#                     rss1 =.$RSS/altScale,
#                     d1=ed1,
#                     d2=ed2)
#     #
#     #new F-test
#     testScaled<-testScaled %>% dplyr::mutate(Fvals=(dRSS/rss1)*(d2/d1))
#     Fvals<-testScaled$Fvals
#     d1<-testScaled$d1
#     d2<-testScaled$d2
#     
#     #scaled values
#     TestScaled<-ggplot(testScaled)+
#       geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
#       geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
#       theme_bw() + 
#       ggplot2::xlab("F-values")+
#       ggplot2::xlim(0,0.05)
#     # print(TestScaled)
#     #Define checked as filtered protein IDs
#     check<-testScaled$uniqueID
#     test<-testScaled 
#     test$d1<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
#     test$d2<-MASS::fitdistr(x=test$RSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
#     
#     testS<-ggplot(test)+
#       geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
#       geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
#       theme_bw() +
#       coord_cartesian(xlim=c(0,10))+
#       ggplot2::xlab("F-values")
#     # print(testS)
#     
#     
#     mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
#     mean1_1<-dplyr::bind_rows(mean1_1)
#     mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
#     mean3<-dplyr::bind_rows(mean3)
#     mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
#     
#     results<-dplyr::bind_rows(mean1,mean1_1,mean3) 
#     if(isTRUE(show_results)){
#       return(testResults)
#       if(isTRUE(scaled_dof)){
#         return(testScaled)
#       }
#     }else{
#       return(results)
#     }
#   }
#   
#   return(results)
# }
spstat<-function(DF,df,df1,Ftest=TRUE,show_results=TRUE,filters=TRUE,scaled_dof=FALSE,Peptide=FALSE){
  if(!any(names(df)=="C")&!any(names(df)=="temperature")){
    DF<-DF %>% dplyr::rename("C"="temperature")
    df<-df %>% dplyr::rename("C"="temperature")
    df1<-df1 %>% dplyr::rename("C"="temperature")
  }
  
  if(any(class(df)=="list")){
    df<-df %>% purrr::keep(function(x) is.data.frame(x))
    df1<-df1 %>% purrr::keep(function(x) is.data.frame(x))
    DF<-DF %>% purrr::keep(function(x) is.data.frame(x))
    
    df<-dplyr::bind_rows(df)
    df1<-dplyr::bind_rows(df1)
    DF<-dplyr::bind_rows(DF)
  }
  #plot spline results
  df1$treatment<-as.factor(df1$treatment)
  df$treatment<-as.factor(df$treatment)
  DF$treatment<-as.factor(DF$treatment)
  
  df1$I<-as.numeric(as.character(df1$I))
  df$I<-as.numeric(as.character(df$I))
  DF$I<-as.numeric(as.character(DF$I))
  
  df1$C<-as.numeric(as.character(df1$C))
  df$C<-as.numeric(as.character(df$C))
  DF$C<-as.numeric(as.character(DF$C))
  #mutate to get CV values
  DF<-DF %>% dplyr::group_split(C,uniqueID) 
  DF<- purrr::map(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df<-df %>% dplyr::group_split(C,uniqueID,treatment) 
  df<- purrr::map(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df1<-df1 %>% dplyr::group_split(C,uniqueID,treatment) 
  df1<- purrr::map(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  
  #convert to data frame
  
  df<-dplyr::bind_rows(df)
  df1<-dplyr::bind_rows(df1)
  DF<-dplyr::bind_rows(DF)
  
  #convert back to list
  if(any(names(DF)=="replicate")&any(names(DF)=="Fraction")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,Fraction,replicate)
    
  }else if(any(names(DF)=="Fraction")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,Fraction)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,Fraction)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,Fraction)
  }else if (any(names(DF)=="replicate")){
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID,replicate)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID,replicate)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID,replicate)
    
  }else{
    DF<-dplyr::bind_rows(DF) %>% dplyr::group_split(uniqueID)
    df<-dplyr::bind_rows(df) %>% dplyr::group_split(uniqueID)
    df1<-dplyr::bind_rows(df1) %>% dplyr::group_split(uniqueID)
  }
  
  #alternative spline fit method : Generalized Additive Models
  #fit penalized splines
  # if(any(names(df[[1]])=="time_point")){
  #   df<-purrr::map(df,function(x) dplyr::bind_rows(x) %>% dplyr::group_split(uniqueID,time_point))
  # }
  fit_gam<-function(x){
    y = x %>% dplyr::filter(!is.infinite(I)) %>% dplyr::mutate(M1 = list(tryCatch(mgcv::gam(x$I ~ s(x$C,by = treatment,k=5), data = x , method = "REML"),
                                                                                  error = function(cond) {
                                                                                    message("Here's the original error message:")
                                                                                    message(cond)
                                                                                    # Choose a return value in case of error
                                                                                    return(NA)})),
                                                               M2 = list(tryCatch(mgcv::gam(x$I ~ s(x$C,by = treatment,k=6), data = x , method = "REML"),
                                                                                  error = function(cond) {
                                                                                    message("Here's the original error message:")
                                                                                    message(cond)
                                                                                    # Choose a return value in case of error
                                                                                    return(NA)}))
    )
    return(y)
  }
  
  populate_fit<-function(x) {
    if(any(stringr::str_detect(names(x),"File.ID"))&!any(names(x)=="sample_id")){
      x<-x %>% dplyr::rename("sample_id"="File.ID")
    }
    x<-x %>% dplyr::group_by(sample_id) %>% dplyr::group_split()
    y<-purrr::map(x,function(x) x %>% dplyr::mutate(pr=list(predict(.$M1[[1]]))))
    y<-purrr::map(y,function(x) x %>% dplyr::mutate(Tm = stats::approx(y[[1]]$pr[[1]],y[[1]]$M1[[1]]$model$`x$C`, xout=min(y[[1]]$pr[[1]],na.rm=TRUE)+(0.5*(max(y[[1]]$pr[[1]], na.rm=TRUE)-min(y[[1]]$pr[[1]], na.rm=TRUE))))$y))
    y<-purrr::map(y,function(x)x %>% dplyr::mutate(uniqueID=.$uniqueID[1],
                                                   k_ = .$M1[[1]]$rank,
                                                   treatment=.$treatment[1],
                                                   rss=deviance(.$M1[[1]]),
                                                   CV_pct = ifelse(!is.null(.$CV_pct),.$CV_pct,NA),
                                                   AUC = pracma::trapz(.$M1[[1]]$fitted.values[(which(abs(.$M1[[1]]$fitted.values-0.5)==min(abs(.$M1[[1]]$fitted.values-0.5)))-1):(which(abs(.$M1[[1]]$fit-0.5)==min(abs(.$M1[[1]]$fitted.values-0.5)))+1)]),
                                                   rsq=summary(x$M1[[1]])$r.sq,
                                                   n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0),
                                                   sample_name=.$sample_name[1],
                                                   missing_pct=ifelse(any(names(x)=="missing_pct"),.$missing_pct[1],NA),
                                                   replicate=replicate) %>% 
                    dplyr::ungroup(.))
    
    return(y)
  }
  if (.Platform$OS.type=="windows"){
    m <- parallel::mclapply(df,fit_gam)
  }else{
    m <- parallel::mclapply(df,fit_gam,mc.cores = future::availableCores())
  }
  
  #remove fitted functions that failed
  m<-m %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(m)==0){
    warning("No fits were possible for the vehicle treatment")
  }
  #parallelize summary of results from gam fit 
  if (.Platform$OS.type=="windows"){
    m <- parallel::mclapply(m,populate_fit)
  }else{
    m <- parallel::mclapply(m,populate_fit,mc.cores = future::availableCores())
  }
  
  if(any(names(m[[1]])=="Annotated_Sequence")){
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,Annotated_Sequence,treatment,sample_id,replicate)
  }else{
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,treatment,sample_id,replicate)
  }
  
  m<-purrr::map(m,function(x)x %>%  
                  dplyr::mutate(fitted_values=list(data.frame(temperature=x$M1[[1]]$model$`x$C`,I=x$M1[[1]]$fitted.values))))
  
  
  if (.Platform$OS.type=="windows"){
    m1 <- parallel::mclapply(df1,fit_gam)
  }else{
    m1 <- parallel::mclapply(df1,fit_gam,mc.cores = future::availableCores())
  }
  m1<-m1 %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(m1)==0){
    warning("No fits were possible for the treated treatment")
  }
  
  #parallelize summary of results from gam fit 
  if (.Platform$OS.type=="windows"){
    m1 <-parallel::mclapply(m1,populate_fit)
  }else{
    m1 <-parallel::mclapply(m1,populate_fit,mc.cores = future::availableCores())
  }
  if(any(names(m1[[1]])=="Annotated_Sequence")){
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,Annotated_Sequence,treatment,sample_id,replicate)
  }else{
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,treatment,sample_id,replicate)
  }
  
  m1<-purrr::map(m1,function(x)x %>%  
                   dplyr::mutate(fitted_values=list(data.frame(temperature=x$M1[[1]]$model$`x$C`,I=x$M1[[1]]$fitted.values))))
  
  
  if (.Platform$OS.type=="windows"){
    mn <- parallel::mclapply(DF,fit_gam)
  }else{
    mn <- parallel::mclapply(DF,fit_gam,mc.cores = future::availableCores())
  }
  mn<-mn %>% purrr::keep(function(x)!any(class(x)=="try-error"))
  if(length(mn)==0){
    warning("No fits were possible for the null treatment")
  }
  #free up memory
  DF<-NA
  df<-NA
  df1<-NA
  #parallelize summary of results from gam fit 
  mn <-furrr::future_map(mn,function(x) populate_fit(x))
  mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
  if(length(mn)==0){
    warning("No fits were possible for the null treatment")
  }else{
    if(any(names(mn[[1]])=="Annotated_Sequence")){
      mn<-dplyr::bind_rows(mn) %>%
        dplyr::group_split(uniqueID,Annotated_Sequence,treatment,sample_id,replicate)
    }else{
      mn<-dplyr::bind_rows(mn) %>%
        dplyr::group_split(uniqueID,treatment,sample_id,replicate)
    }
    
    mn<-purrr::map(mn,function(x)x %>%  
                     dplyr::mutate(fitted_values=list(data.frame(temperature=x$M1[[1]]$model$`x$C`,I=x$M1[[1]]$fitted.values))))
  }
  
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
  
  m<-dplyr::bind_rows(m) %>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  m1<-dplyr::bind_rows(m1)%>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  mn<-dplyr::bind_rows(mn)%>%
    distinct(.)%>%
    dplyr::group_split(uniqueID)
  #calculate RSS
  m<-purrr::map2(m,m1,function(x,y) x %>% dplyr::mutate(RSSv=as.numeric(x$rss[x$treatment=="vehicle"][1]),
                                                        RSSt=as.numeric(y$rss[y$treatment=="treated"][1]),
                                                        RSS1=as.numeric(x$rss[x$treatment=="vehicle"][1])+as.numeric(y$rss[y$treatment=="treated"][1]),
                                                        treatment=x$treatment[1],
                                                        Tm=x$Tm[1]
  ))
  
  m1<-purrr::map2(m,m1,function(x,y) y %>% dplyr::mutate(RSSv=x$RSSv[1],
                                                         RSSt=x$RSSt[1],
                                                         RSS1=x$RSS1[1],
                                                         Tm=y$Tm[1]
  ))
  mn<-purrr::map(mn,function(x) x %>% dplyr::mutate(RSSv=NA,
                                                    RSSt=NA,
                                                    RSS0=x$rss,
                                                    Tm=x$Tm[1]
  ))
  if(any(names(m)=="Fraction")){
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,Fraction)
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,Fraction)
    mn<-dplyr::bind_rows(mn) %>%
      dplyr::group_split(uniqueID,Fraction)
  }else if(any(names(m)=="replicate")){
    m<-dplyr::bind_rows(m) %>%
      dplyr::group_split(uniqueID,replicate)
    m1<-dplyr::bind_rows(m1) %>%
      dplyr::group_split(uniqueID,replicate)
    mn<-dplyr::bind_rows(mn) %>%
      dplyr::group_split(uniqueID,replicate)
  }
  m<-purrr::map2(m,mn,function(x,y) x %>%
                   dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  m1<-purrr::map2(m1,mn,function(x,y)x %>%
                    dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  mn<-purrr::map2(m1,mn,function(x,y)y %>%
                    dplyr::mutate(rssDiff=y$RSS0[1]-x$RSS1[1]))
  #4031 proteins have rssDiff> 0 
  
  #convert to df and split by uniqueID 
  mean1<-dplyr::bind_rows(m)
  mean1_1<-dplyr::bind_rows(m1)
  mean3<-dplyr::bind_rows(mn)
  #equal column names
  mean1<-mean1 %>% dplyr::select(-RSSv,-RSSt,-RSS1)%>% distinct(.)
  if(any(names(mean1_1)=="RSS1")){
    mean1_1<-mean1_1  %>% dplyr::select(-RSSv,-RSSt,-RSS1)%>% distinct(.)
  }
  mean3<-mean3  %>% dplyr::select(-RSS0,-RSSv,-RSSt) %>% distinct(.)
  
  #filter out proteins with neg RSSdiff
  if(isTRUE(filters)){
    check<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
    
    mean1<-mean1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    mean3<-mean3 %>% dplyr::filter(uniqueID %in% check,rsq>0.8)
    
    
    mean1<-mean1 %>% dplyr::filter(uniqueID %in% check)
    mean1_1<- mean1_1 %>% dplyr::filter(uniqueID %in% check)
    mean3<-mean3 %>% dplyr::filter(uniqueID %in% check)
    
  }
  
  #split into lists by uniqueID
  mean1<-mean1%>% dplyr::group_split(uniqueID)
  mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
  mean3<-mean3 %>% dplyr::group_split(uniqueID)
  
  #Cliff
  results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>%
    dplyr::group_split(uniqueID,replicate)
  results<-purrr::map(results,function(x) x %>%
                        dplyr::mutate(dTm=x[which(x$treatment=="treated"),]$Tm[1]-x[which(x$treatment=="vehicle"),]$Tm[1]))
  
  #only keep data with vehicle and treated values
  results<-dplyr::bind_rows(results) %>%
    dplyr::group_by(uniqueID,replicate) %>%
    dplyr::group_split(.)
  #make sure we keep proteins with both treated and vehicle data
  results<-results%>%
    purrr::keep(function(x) length(unique(x$treatment))>1)
  if(!isTRUE(Peptide)){
    results_<-dplyr::bind_rows(results) %>% dplyr::group_split(uniqueID,sample_name,replicate)
    results_<-results_ %>% purrr::keep(function(x) nrow(x)>1)
    results_<-purrr::map(results_,function(x) x %>%dplyr::group_by(replicate) %>% dplyr::mutate(dTm=(.$Tm[which(.$treatment=="treated")]-.$Tm[which(.$treatment=="vehicle")])))
  }else{#if this is a peptide file,calculate dTm for each PSM
    results_<-dplyr::bind_rows(results) %>%
      dplyr::select(uniqueID,Annotated_Sequence,sample_name,sample_id,treatment,C,I,replicate) %>%
      distinct(.) %>% dplyr::group_split(uniqueID,Annotated_Sequence,sample_name,sample_id)
    
    results_<-dplyr::bind_rows(results_) %>%dplyr::group_split(uniqueID,replicate)
    results_<-purrr::map(results_,function(x) x %>% dplyr::mutate(dTm=(.$Tm[.$treatment=="treated"][1]-.$Tm[.$treatment=="vehicle"][1])))
    
  }
  
  results_2<-dplyr::bind_rows(results_) %>%
    dplyr::group_split(uniqueID,treatment)
  #nest the columns that dont involve Tm t-test calculation 
  nesting<-names(results_2[[1]])[!names(results_2[[1]]) %in% c("uniqueID","sample_id","treatment","Tm","replicate","sample_name","dTm","missing_pct")]
  
  results_2<-purrr::map(results_2,function(x) x %>% tidyr::nest(data=nesting))
  #only keep replicates with Tm values
  results_2<-dplyr::bind_rows(results_2) %>%
    dplyr::group_split(uniqueID)
  results_2<-purrr::map(results_2,function(x) x %>%
                          dplyr::filter(replicate %in% unique(x$replicate[duplicated(x$replicate)])))
  
  results_3<-dplyr::bind_rows(results_2) %>%  dplyr::group_split(uniqueID,replicate)
  
  results_<-results_3 %>% purrr::keep(function(x) any(!is.na(x$dTm)))
  if(length(results_)==0){
    warning(paste0(results_2[[1]]$sample_name[1], " has at least one Tm value missing"))
  }
  results_3<-purrr::map(results_3,function(x) x %>%
                          distinct(.) %>%
                          dplyr::mutate(
                            dTm=x$Tm[x$treatment %in% "treated"][1]-x$Tm[x$treatment %in% "vehicle"][1])%>% dplyr::ungroup(.)) 
  results_<-results_3 %>% purrr::keep(function(x) !is.na(x$dTm[1]))
  results_<-dplyr::bind_rows(results_) %>% dplyr::group_split(sample_name)
  results_<-results_ %>% purrr::keep(.,function(x) length(unique(x$treatment))>1)
  
  results_<-purrr::map(results_,function(x) 
    y<-x %>% 
      dplyr::mutate(hypothesis=ifelse(treatment=="null",as.factor("null"),as.factor("alternative")),
                    p_dTm= try(p.adjust(t.test(Tm ~ treatment, data = .,
                                               var.equal = FALSE,conf.level=0.975)$p.value[1],"BH"))))
  
  
  
  results1<-results_ %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
  
  results<-dplyr::bind_rows(results1) %>% tidyr::unnest(cols=data)
  #results<-results %>%  dplyr::mutate(p_dTm=calcP(uniqueID,Tm,dTm,30000))
  results2<-dplyr::bind_rows(results)$uniqueID
  if(!length(mean3)==length(mean1)){
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    
    CID<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-dplyr::intersect(CID,mean3$uniqueID)
    CID<-dplyr::intersect(CID,results2)
    
    #filter
    mean1 <-mean1 %>% dplyr::filter(uniqueID %in% CID)
    mean1_1<-mean1_1%>% dplyr::filter(uniqueID %in% CID)
    mean3<-mean3%>% dplyr::filter(uniqueID %in% CID)
    
    #split 
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1%>% dplyr::group_split(uniqueID)
    mean3<-mean3%>% dplyr::group_split(uniqueID)
    
    if(!nrow(mean3[[1]])==nrow(mean1[[1]])){
      mean1<-dplyr::bind_rows(mean1) %>% distinct(.)
      mean1_1<-dplyr::bind_rows(mean1_1)%>% distinct(.)
      mean3<-dplyr::bind_rows(mean3)%>% distinct(.)
      
      CID<-dplyr::intersect(mean1$uniqueID,mean1_1$uniqueID)
      CID<-dplyr::intersect(CID,mean3$uniqueID)
      #filter
      mean1 <-mean1 %>% dplyr::filter(uniqueID %in% CID)
      mean1_1<-mean1_1%>% dplyr::filter(uniqueID %in% CID)
      mean3<-mean3%>% dplyr::filter(uniqueID %in% CID)
      
      #split 
      mean1<-mean1 %>% dplyr::group_split(uniqueID)
      mean1_1<-mean1_1%>% dplyr::group_split(uniqueID)
      mean3<-mean3%>% dplyr::group_split(uniqueID)
    }
    results<-dplyr::bind_rows(mean1,mean1_1,mean3)
  }
  
  if(isTRUE(Ftest)){
    stats_summary_null<-function(x){
      y =data.frame(uniqueID=x$uniqueID[1],
                    sample_id=x$sample_id[1],
                    sample_name=x$sample_name[1],
                    replicate=x$replicate[1],
                    AUC=x$AUC[1],
                    rsq=x$rsq[1],
                    RSS0=x$rss[1],
                    se0=summary.gam(x$M1[[1]])$se[1],#standard error
                    edf0=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                    rsq0=summary.gam(x$M1[[1]])$r.sq[1],#r-squared
                    np0=summary.gam(x$M1[[1]])$np[1],#number of parameters
                    rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m=summary.gam(x$M1[[1]])$m[1]) %>% distinct(.)
    }
    #Calculate stats summary
    rss0<-purrr::map(mean3,function(x)stats_summary_null(x))
    #this is the alternative hypothesis (one curve per treatment)
    stats_summary_alt<-function(x,y){
      z =data.frame(uniqueID=x$uniqueID[1],
                    sample_id=x$sample_id[1],
                    sample_name=x$sample_name[1],
                    replicate=x$replicate[1],
                    AUC=x$AUC[1],
                    dTm=y$Tm[1]-x$Tm[1],
                    RSS=x$rss[1]+y$rss[1],
                    se=summary.gam(x$M1[[1]])$se[1],#standard error
                    pN1=summary.gam(x$M1[[1]])$edf[1],#effective degrees of freedom
                    rdf=summary.gam(x$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m=summary.gam(x$M1[[1]])$m[1],
                    se1=summary.gam(y$M1[[1]])$se[1],#standard error
                    pN2=summary.gam(y$M1[[1]])$edf[1],#effective degrees of freedom treated
                    rsq=mean(x$M1[[1]]$r.sq[1],y$M1[[1]]$r.sq[1],na.rm=TRUE),#r-squared
                    rdf1=summary.gam(y$M1[[1]])$residual.df[1],#residual degrees of freedom
                    m1=summary.gam(y$M1[[1]])$m[1],
                    npa=sum(summary.gam(x$M1[[1]])$np[1],summary.gam(y$M1[[1]])$np[1],na.rm=TRUE),
                    pNA=sum(summary.gam(x$M1[[1]])$edf[1],summary.gam(y$M1[[1]])$edf[1],na.rm=TRUE)) %>% distinct(.)
      return(z)
    }
    rss1<-purrr::map2(mean1,mean1_1,function(x,y)stats_summary_alt(x,y))
    
    # rss1<-purrr::map2(rss1,mean1,function(x,y) x %>% dplyr::mutate(fit_v=list(ifelse(class(try(mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",
    #                                                                                                              se.fit=TRUE,newdata.guaranteed = TRUE)))=='try-error',NA,
    #                                                                                  mgcv::predict.gam(y$M1[[1]],newdata=data.frame(C=y$newdata[[1]]),family="link",se.fit=TRUE,newdata.guaranteed = TRUE)))))
    int<-dplyr::intersect(names(mean3[[1]]),names(rss0[[1]]))
    int1<-dplyr::intersect(names(mean1[[1]]),names(rss1[[1]]))
    
    
    # mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
    # mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
    # 
    #bind predicted values to main data frame
    mean3<-purrr::map2(mean3,rss0,function(x,y)x %>% dplyr::right_join(y))
    mean1<-purrr::map2(mean1,rss1,function(x,y)x %>% dplyr::right_join(y))
    
    
    #params for null and alternative models
    f0<-lapply(mean3,function(x) data.frame(nfv=length(x$M1[[1]]$fitted.values)))#number of measurements
    if(length(mean1[[1]]$M1[[1]]$fitted.values)<length(mean1_1[[1]]$M1[[1]]$fitted.values)){
      f1<-purrr::map2(mean1,mean1_1,function(x,y) sum(length(x$M1[[1]]$fitted.values),
                                                      length(y$M1[[1]]$fitted.values)))#number of measurements alternative
    }else{
      f1<-purrr::map2(mean1,mean1_1,function(x,y) sum(length(x$M1[[1]]$fitted.values),
                                                      length(y$M1[[1]]$fitted.values)))#number of measurements alternative
    }
    f1<-purrr::map(f1,function(x) data.frame(nfvA=sum(x)))
    int<-dplyr::intersect(names(mean3[[1]]),names(f0[[1]]))
    int1<-dplyr::intersect(names(mean1[[1]]),names(f1[[1]]))
    
    # 
    # mean3<-purrr::map(mean3,function(x)x %>% dplyr::select(-all_of(int)))
    # mean1<-purrr::map(mean1,function(x) x %>% dplyr::select(-all_of(int1)))
    
    #bind data frames
    mean3<-purrr::map2(mean3,f0,function(x,y) merge(x,y,all=TRUE))
    mean1<-purrr::map2(mean1,f1,function(x,y) merge(x,y,all=TRUE))
    mean1_1<-purrr::map2(mean1_1,f1,function(x,y) merge(x,y,all=TRUE))
    #
    pN<-NA
    pA<-NA
    
    d1<-NA
    d2<-NA
    
    #calculate parameters 
    pN<-purrr::map(mean3,function(x) x[1,])
    pA<-purrr::map(mean1,function(x)x[1,])
    
    
    #degrees of freedom before
    d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pNA-y$np0))
    d2<-purrr::map(pA,function(x)data.frame(d2=x$npa-x$pNA))
    
    #delta RSS
    rssDiff<-purrr::map2(rss0,rss1,function(x,y) data.frame(dRSS=x$RSS0-y$RSS,#rss null minus rss alt
                                                            dTm=y$dTm[1],
                                                            rss1=y$RSS[1]))
    
    d2<-lapply(d2,function(x) x$d2[1])#edf2
    d1<-lapply(d1,function(x) x$d1[1])#edf1
    
    #RSS1 numerator
    Fvals<-purrr::map2(rssDiff,d1,function(x,y) data.frame(fNum=x$dRSS/y[1]))
    #Rss denominator
    Fd<-purrr::map2(rss1,d2,function(x,y) data.frame(fDen=x$RSS/y[1]))
    
    Fvals<-purrr::map2(Fvals,Fd,function(x,y) data.frame(fStatistic=x$fNum/y$fDen))
    Fvals<-purrr::map2(Fvals,d1,function(x,y) x %>% dplyr::mutate(df1=y[1]))
    Fvals<-purrr::map2(Fvals,d2,function(x,y) x %>% dplyr::mutate(df2=y[1]))
    
    Fvals<-purrr::map2(Fvals,rssDiff,function(x,y) cbind(x,y))
    Fvals<-purrr::map(Fvals,function(x) x %>% dplyr::mutate(Fvals=(x$dRSS/x$rss1)*(x$df2/x$df1)))
    #calculate p and p-adj vals
    
    Fvals<-purrr::map(Fvals,tryCatch({function(x) x %>% dplyr::mutate(pValue = 1 - pf(fStatistic, df1 = x$df1, df2 = x$df2),
                                                                      pAdj = p.adjust(pValue,"BH"))
    },error = function(cond){
      message(cond)
    }))
    Fvals<-purrr::map2(Fvals,mean1,function(x,y) x %>% dplyr::mutate(uniqueID=y$uniqueID[1]))
    
    
    mean1<-purrr::map(mean1,function(x) data.frame(x)%>% dplyr::select(-M1) )
    
    mean1<-purrr::map(mean1,function(x) x %>% distinct(.) )
    names1<-dplyr::intersect(names(mean1),names(Fvals))
    #convert results to list
    testResults<-purrr::map2(mean1,Fvals,function(x,y)x%>% dplyr::right_join(y,by=names1))
    testResults<-purrr::map(testResults,function(x) x[1,])
    testResults<-dplyr::bind_rows(testResults)
    mean1<-dplyr::bind_rows(mean1)
    Unscaled<-ggplot(testResults)+
      geom_density(aes(x=fStatistic),fill = "steelblue",alpha = 0.5) +
      geom_line(aes(x=fStatistic,y= df(fStatistic,df1=df1,df2=df2)),color="darkred",size = 1.5) +
      theme_bw() +
      ggplot2::xlab("F-values")+
      ggplot2::xlim(0,0.05)
    # print(Unscaled)
    #scale variables
    M<-median(testResults$dRSS,na.rm=TRUE)
    V<-mad(testResults$dRSS,na.rm=TRUE)
    #alternative scaling factor sig0-sq
    altScale<-0.5*V/M
    #filter out negative delta rss
    testResults<-testResults %>% dplyr::filter(dRSS>0)
    #effective degrees of freedom
    ed1<-tryCatch({MASS::fitdistr(x=testResults$dRSS, densfun = "chi-squared", start = list(df=2))[["estimate"]]},
                  error= function (cond){message(cond)})
    ed2<-tryCatch({MASS::fitdistr(x=testResults$rss, densfun = "chi-squared", start = list(df=2))[["estimate"]]},
                  error = function(cond){message(cond)})
    #scale data
    testScaled <-testResults %>%
      dplyr::mutate(rssDiff = .$dRSS/altScale,
                    rss1 =.$RSS/altScale,
                    d1=ed1,
                    d2=ed2)
    #
    #new F-test
    if(!class(ed1)=="NULL"&!class(ed2)=="NULL"){
      testResults<-testScaled %>% dplyr::mutate(Fvals=(dRSS/rss1)*(d2/d1))
      Fvals<-testResults$Fvals
      d1<-testResults$d1
      d2<-testResults$d2
      
      #scaled values
      TestScaled<-ggplot(testResults)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() + 
        ggplot2::xlab("F-values")+
        ggplot2::xlim(0,0.05)
      # print(TestScaled)
      #Define checked as filtered protein IDs
      check<-testResults$uniqueID
      test<-testResults
      test$d1<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test$d2<-MASS::fitdistr(x=test$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test<-test %>% dplyr::filter(test$pAdj<0.01)
      testS<-ggplot(test)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
      # print(testS)
      
      
      mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
      mean3<-dplyr::bind_rows(mean3)
      mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) 
      if(!isTRUE(Peptide)){
        results1<-dplyr::bind_rows(results) %>%
          dplyr::select(uniqueID,treatment,sample_id,p_dTm)
        nam<-dplyr::intersect(names(results),names(results1))
        results<-results %>% dplyr::right_join(results1,by=nam)
        
        nam<-dplyr::intersect(names(testResults),names(results1))
        testResults<-testResults %>% dplyr::right_join(results1,by=nam)
        
        names<-dplyr::intersect(names(test),names(results1))
        test<-test %>% dplyr::right_join(results1,by=names)
      }else{
        results1<-dplyr::bind_rows(results1) %>% dplyr::select(uniqueID,treatment,sample_id,p_dTm)
        names<-dplyr::intersect(names(results),names(results1))
        results<-results %>% dplyr::right_join(results1,by=names)
        
        names<-dplyr::intersect(names(testResults),names(results1))
        testResults<-testResults %>% dplyr::right_join(results1,by=names)
        
        names<-dplyr::intersect(names(test),names(results1))
        test<-test %>% dplyr::right_join(results1,by=names)
      }
    }
    
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) 
    # results<-results %>% dplyr::group_by(treatment) %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(performance_k5 = list(performance::model_performance(unlist(.$M1))),
    #                 performance_k6 = list(performance::model_performance(unlist(.$M2))))
    
  }
  if(isTRUE(show_results)&exists("testResults")){
    
    if(isTRUE(scaled_dof)){
      return(test)
    }
    
    return(testResults)
  }else{
    
    return(results)
  }
  
}

spf<-function(spresults,DFN,filters = TRUE){
  spresults<-spresults %>% 
    dplyr::mutate(replicate=as.numeric(replicate)) %>%
    dplyr::group_split(uniqueID)
  DFN<-DFN %>% dplyr::mutate(replicate=as.numeric(replicate)) 
  if(!isTRUE(filters)){
    sl<-purrr::map(seq_len(length(spresults)),function(x) as.numeric({paste(x)})) 
    sp<-purrr::map2(spresults,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))  
    sp<-dplyr::bind_rows(sp)  
    df1<-data.frame(uniqueID = unique(sp$uniqueID))  
    df2<-dplyr::bind_rows(DFN)  
    df2$C<-as.numeric(as.vector(df2$C)) 
    df2$I<-as.numeric(as.vector(df2$I))  
    
    if(any(names(df2)=="id")){
      df2<-df2 %>% dplyr::select(-id)
    }
    
    names<-intersect(names(df2),names(sp))
    df2<-sp %>% left_join(df2, by = names) 
    
    
    Df1<-spresults
  }else{
    #Apply filters 
    #keep the positive AUC differences
    
    spresults<-spresults %>% keep(function(x) mean(x$AUC[x$treatment=="treated"],na.rm=TRUE)>mean(x$AUC[!x$treatment=="vehicle"],na.rm=TRUE))
    
    spresults<-spresults %>% keep(function(x) max(x$lambda)<1)
    
    if (is.null(nrow(spresults))){
      return(warning("all proteins filtered out by AUC and lambda value"))
    }
    #get Tm and RSS differences
    sp<-purrr::map(spresults, function(x) x %>% dplyr::mutate(Tmd= x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "treated"),'Tm'][[1]] - x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "vehicle"),'Tm'][[1]],
                                                              RSSd = sum(x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']) - sum(x[!stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "null"),'rss']),
                                                              AUCd = x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "treated"),'AUC'])[[1]]- x[stringr::str_detect(tolower(data.frame(x)$treatment), pattern = "vehicle"),'AUC'][[1]])
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
      if(any(names(df2)=="id")){
        df2<-df2 %>% dplyr::select(-id)
      }
      names<-intersect(names(df2),names(sp))
      df2<-sp %>% left_join(df2, by = names) 
    }else{
      if(any(names(df2)=="id")){
        df2<-df2 %>% dplyr::select(-id)
      }
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

spCI<-function(i,df1,df2,Df1,df.temps,overlay=TRUE,alpha,residuals=FALSE,simulations=FALSE,CI=TRUE,Peptide=FALSE,CARRIER=TRUE,Protein=Protein,raw=FALSE){
  null<-data.frame()
  i<-i
  if(isTRUE(Peptide)){
    df2<-dplyr::bind_rows(df2)
    df2<-df2[,!stringr::str_detect(names(df2),"File.ID|Channel|RT|Confidence|Protein|p|Percolator|DeltaM|Tm|rsq|CC|k_|AUC")]
    df2<-df2[,!names(df2)=="rss"]
    df2<-df2%>% distinct(.) %>%
      dplyr::group_split(uniqueID,replicate,treatment,C)
    Df1<-dplyr::bind_rows(Df1) %>%
      distinct(.) %>%
      dplyr::group_split(uniqueID,treatment,C)
    
  }
  df2<-dplyr::bind_rows(df2) %>%
    distinct(.)
  Df1<-dplyr::bind_rows(Df1) %>%
    dplyr::group_split(uniqueID) 
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>% dplyr::mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-as.character(df2$uniqueID)
  
  i<-which(df1$uniqueID %in% Protein)
  #get original data
  ###########################################
  df1<-df1$uniqueID[i]
  DF1<-df2[which(df2$uniqueID %in% df1),]
  Df1<-data.frame(Df1[[i]])
  null<-Df1[which(Df1$uniqueID %in% df1 & Df1$treatment %in% "null"),]
  if(nrow(null)==0){
    null<-Df1 %>% dplyr::mutate(treatment="null")
  }
  ###########################################
  DF_f<-df2 %>% subset(uniqueID %in% df1 & treatment %in% "vehicle")
  vehicle<-Df1 %>% subset(uniqueID %in% df1 & treatment %in% "vehicle")
  ###########################################
  DF_f1<-df2%>% subset(uniqueID %in% df1 & treatment %in% "treated")
  treated<-Df1 %>% subset(uniqueID == df1 & treatment == "treated")
  
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
  df2<-df2 %>% dplyr::filter(uniqueID==as.character(df1))
  BSvarN<-df2 %>% dplyr::mutate(treatment=="null")
  BSvar1 <-df2 %>% dplyr::filter(treatment=="treated")
  BSvar <-df2 %>% dplyr::filter(treatment=="vehicle")
  if(nrow(BSvar)==0){
    return(warning(paste0("No vehicle data found for ",Protein)))
  }
  if(nrow(BSvar1)==0){
    return(warning(paste0("No treated data found for ",Protein)))
  }
  BSvarN$treatment<-as.factor(BSvarN$treatment)
  BSvar$treatment<-as.factor(BSvar$treatment)
  BSvar1$treatment<-as.factor(BSvar1$treatment)
  #####try GAM
  #fit penalized splines
  m <- mgcv::gam(I ~ s(C,k=5), data = BSvar , method = "ML")
  m1<-  mgcv::gam(I ~ s(C,k=5), data =BSvar1, method = "ML")
  mn<-  mgcv::gam(I ~ s(C,k=5), data = BSvarN, method = "ML")
  
  #####try GAM
  
  #get some parmeters
  Vb <- vcov(m)
  newd <- with(BSvar, data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 30)))%>% as.data.frame(.)
  pred <- predict(m, newd, se.fit = TRUE)%>% as.data.frame(.)#get confidence intervals
  se.fit <- pred$se.fit
  #get some parmeters
  Vb1<- vcov(m1) 
  newd1<- with(BSvar1,data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 30)))%>% as.data.frame(.)
  pred1<- predict(m1,newd1,se.fit = TRUE) %>% as.data.frame(.)
  se.fit1<- pred1$se.fit
  #generate seed for randomization
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
  pred$treatment<-"vehicle"
  pred$treatment<-as.factor(pred$treatment)
  pred$CI<-"vehicle"
  pred$CI<-as.factor(pred$CI)
  
  plot<-ggplot(pred,mapping= ggplot2::aes(x = C,y=fit,color=treatment))+
    geom_point(BSvar, mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(replicate)))+
    geom_ribbon(aes(ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle("")+ylim(-0.4,1.6)+xlim(37,68)+theme(legend.position="bottom")
  
  pred1<- transform(cbind(data.frame(pred1),newd1),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit))
  pred1$treatment<-"treated"
  pred1$treatment<-as.factor(pred1$treatment)
  pred1$CI<-"treated"
  pred1$CI<-as.factor(pred1$CI)
  pred1$AUC<-pracma::trapz(pred1$fit-pred$fit)
  pred1$AUC<-abs(round(pred1$AUC[1],3))
  
  pred1$RSS<- deviance(m1)+deviance(m)
  
  pred1$RSS<- round(pred1$RSS,3)
  #Residuals
  
  pred1$Tm<-round(treated$Tm[1]-vehicle$Tm[1],1)
  #missing values
  miss_v<-data.frame(NA)
  miss_t<-data.frame(NA)
  #max replicates
  
  miss_v<-DF1%>% dplyr::filter(treatment=="vehicle") %>% unique(.)
  miss_t<-DF1%>% dplyr::filter(treatment=="treated") %>% unique(.)
  
  Pred<-data.frame(m$fitted.values,m$residuals)
  names(Pred)<- c("fit","rn")
  Pred$treatment<-as.factor("vehicle")
  BSvar$treatment<-as.factor("vehicle")
  Pred1<-data.frame(m1$fitted.values,m1$residuals)
  names(Pred1)<-c("fit","rn")
  Pred1$treatment<-as.factor("treated")
  Preds<-rbind(Pred,Pred1)
  BSvar1$treatment<-as.factor("treated")
  #get fitted value data
  fitted.values<-data.frame(C=BSvar$M1[[1]]$model$`x$C`,fit=predict(BSvar$M1[[1]],se.fit=TRUE))
  names(fitted.values)<-c("C","fit","se.fit")
  fitted.values1<-data.frame(C=BSvar1$M1[[1]]$model$`x$C`,fit=predict(BSvar1$M1[[1]],se.fit=TRUE))
  names(fitted.values1)<-c("C","fit","se.fit")
  
  #append Tm values on predicted data
  pred1$Tm<-round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)
  if(isTRUE(residuals)){
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=treatment)) +ggplot2::geom_point()+ 
      ggplot2::ggtitle(paste(Df1$uniqueID[1]," ",stringr::str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    print(PLrs)
  }
  if(isTRUE(CI)){
    BSVar <-df2 %>% subset(uniqueID == df1 & treatment== "vehicle") # %>% dplyr::mutate(I=mean(I))
    BSVar1 <-df2 %>% subset(uniqueID == df1 & treatment== "treated")# %>% dplyr::mutate(I=mean(I))
    BSVarN<-df2 
    # BSVar<-BSVar[!is.na(BSVar$I),]
    # BSVar1<-BSVar1[!is.na(BSVar1$I),]
    
    BSVar<-BSVar %>% distinct(.)
    BSVar1<-BSVar1 %>% distinct(.)
    
    
    #fit penalized splines
    m <- mgcv::gam(I ~ s(C,k=5), data = BSVar , method = "ML")
    m1<-  mgcv::gam(I ~ s(C,k=5), data =BSVar1, method = "ML")
    mn<-  mgcv::gam(I ~ s(C,k=5), data = BSVarN, method = "ML")
    
    #####try GAM
    #get some parmeters
    Vb <- vcov(m)
    newd <- with(BSVar, data.frame(C = seq(from=min(BSVar$C,na.rm=TRUE), to=max(BSVar$C,na.rm=TRUE), length.out = 10)))%>% as.data.frame(.)
    BSVar <- BSVar %>% dplyr::mutate(fit=list(predict(m, newd, se.fit = TRUE)))
    
    #get some parmeters
    Vb1<- vcov(m1) 
    newd1<- with(BSVar1,data.frame(C = seq(min(C,na.rm=TRUE), max(C,na.rm=TRUE), length = 10)))%>% as.data.frame(.)
    BSVar1 <- BSVar1 %>% dplyr::mutate(fit=list(predict(m1, newd1, se.fit = TRUE)))
    # if (any(names(BSVar)=="sample.x")){
    #   BSVar<-BSVar %>% dplyr::rename("sample"="sample.x")
    #   
    # }
    # if (any(names(BSVar1)=="sample.x")){
    #   BSVar1<-BSVar1 %>% dplyr::rename("sample"="sample.x")
    #   
    # }
    # 
    #append missing value data
    if(!isTRUE(Peptide)){
      BSVar<-BSVar %>% dplyr::mutate(missing_v=round(BSVar$missing_pct[!is.na(BSVar$missing_pct)][1],0))
      BSVar<-BSVar %>% dplyr::mutate(missing_t=round(BSVar1$missing_pct[!is.na(BSVar1$missing_pct)][1],0))
    }
    p<-data.frame(BSVar$fit[[1]])
    p1<-data.frame(BSVar1$fit[[1]])
    
    fit_v<-p %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                               uprP=fit+(1.96*se.fit),
                               uprS = fit + (crit * se.fit),
                               lwrS = fit - (crit * se.fit),
                               C= seq(min(BSVar$C,na.rm=TRUE), max(BSVar$C,na.rm=TRUE), length = 10),
                               treatment=BSVar$treatment[1],
                               CI=BSVar$treatment[1])
    fit_t<-p1 %>% dplyr::mutate(lwrP=fit-(1.96*se.fit),
                                uprP=fit+(1.96*se.fit),
                                uprS = fit + (crit1 * se.fit),
                                lwrS = fit - (crit1 * se.fit),
                                C=seq(min(BSVar1$C,na.rm=TRUE), max(BSVar1$C,na.rm=TRUE), length = 10),
                                treatment=BSVar1$treatment[1],
                                CI=BSVar1$treatment[1])
    # 
    # id<-BSVar %>%dplyr::ungroup(.) %>%dplyr::select(sample,replicate)%>% distinct(.)
    # id1<-BSVar1 %>%dplyr::ungroup(.) %>%dplyr::select(sample,replicate)%>% distinct(.)
    # 
    # BSVar$replicate<-as.factor(BSVar$replicate)
    # BSVar1$replicate<-as.factor(BSVar1$replicate)
    # id$replicate<-as.factor(id$replicate)
    # id1$replicate<-as.factor(id1$replicate)
    # BSVar<-BSVar %>% dplyr::right_join(id,by=c("sample","replicate"))
    # BSVar1<-BSVar1%>% dplyr::right_join(id1,by=c("sample","replicate"))
    
    BSVar<-dplyr::bind_rows(BSVar) %>% distinct(.)
    BSVar1<-dplyr::bind_rows(BSVar1) %>% distinct(.)
    if(isTRUE(Peptide)){
      if(any(names(BSVar)=="replicate")&isTRUE(Peptide)){
        BSVar$Replicate<-as.factor(BSVar$Replicate)
        BSVar1$Replicate<-as.factor(BSVar1$Replicate)
      }else if(any(names(BSVar)=="Charge")){
        BSVar$PSMs_Charge<-as.factor(BSVar$Charge)
        BSVar1$PSMs_Charge<-as.factor(BSVar1$Charge)
      }
      if(any(names(BSVar)=="rank_l")){
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
      }else{
        BSVar$Stroke<-0
        BSVar1$Stroke<-0
      }
      if(!isTRUE(raw)&!all(BSVar$Stroke)==0){#if I isnt raw and there's no rank column for the data
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, fill=NA, size=4, colour="black", stroke=1)+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
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
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, fill=NA, size=4)+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.60,2)+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }else if(!isTRUE(raw)&BSVar$Stroke==0){#if I is raw and there's no rank column for the data
        BSVar$Num_PSMs<-as.factor(BSVar$Num_PSMs)#the number of PSMs found for this particular annotated_Sequence
        BSVar1$Num_PSMs<-as.factor(BSVar1$Num_PSMs)
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point()+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
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
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point()+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.60,2)+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }else{#if there is a rank column present, use the stroke to label low-ranking Peptides
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, size=4, colour="black", stroke=1)+
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_v),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", pred1$RSS[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
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
          ggplot2::geom_point(BSVar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, size=4, color=treatment)+
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=treatment), alpha = 0.2 ) +
          ggplot2::geom_ribbon(data.frame(fit_t),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(min(c(min(BSVar$I,na.rm=TRUE),min(BSVar1$I,na.rm=TRUE))),0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+xlim(37,68)+
          theme(legend.position="bottom")
        
        return(plot)
      }
      
    }else{#if this is a protein file
      if(!any(stringr::str_detect(names(BSVar),"replicate"))){
        BSVar$Replicate<-as.factor(BSVar$replicate)
        BSvar1$Replicate<-as.factor(BSVar1$replicate)
      }
      if(!any(stringr::str_detect(names(BSVar),"Replicate"))){#if there is no replicate column
        BSVar<-BSVar %>% dplyr::group_by(C,treatment) %>% distinct(.) %>%dplyr::group_split()
        BSVar<-lapply(BSVar,function(x) x %>% dplyr::mutate(replicate=row.names(x)))
        BSVar<-dplyr::bind_rows(BSVar) %>% dplyr::ungroup()
        
        BSvar1<-BSvar1  %>% dplyr::group_by(C,treatment) %>% distinct(.) %>%dplyr::group_split()
        BSvar1<-lapply(BSvar1,function(x) x %>% dplyr::mutate(replicate=row.names(x)))
        BSvar1<-dplyr::bind_rows(BSvar1) %>% dplyr::ungroup()
        
        BSVar$Replicate<-as.factor(BSVar$replicate)
        BSvar1$Replicate<-as.factor(BSvar1$replicate)
      }
      if(any(names(BSVar)=="rank_l")){#if this protein file has a rank column
        
        BSVar$Stroke<-ifelse(BSVar$rank_l==TRUE,1,0)
        BSVar1$Stroke<-ifelse(BSVar1$rank_l==TRUE,1,0)
        
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar[BSVar$Stroke==1,],
                     pch=21, size=4, color=treatment, stroke=1)+
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",signif(pred1$Tm[1],3),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        
        
        plot<-plot1+
          ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=factor(Replicate)))+
          geom_point(data=BSVar1[BSVar1$Stroke==1,],
                     pch=21, size=4, color=treatment, stroke=1)+
          ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred1),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = 2+round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   y = -0.10,
                   label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5,
                   yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                   xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0,
                   yend = 0.5,
                   colour = "red",linetype=2)+ ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.75,0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+xlim(37,68)+
          theme(legend.position="bottom")
        return(plot)
      }else{#if there is no rank column for this protein file
        
        plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x=C,y=I,color=treatment))+
          ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = treatment,shape=factor(Replicate)))+
          ggplot2::geom_point()+
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0)+
          ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")
        
        
        plot<-plot1+
          ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color=treatment,shape=Replicate))+
          ggplot2::geom_point()+
          ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
          ggplot2::geom_ribbon(data.frame(pred1),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrS, ymax = uprS ,fill=CI), alpha = 0.2,linetype=0) +
          ggplot2::labs(y = "Relative Solubility",
                        x = "Temperature (\u00B0C)")+
          annotate("text",
                   x = signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="red",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values1$C), xend = signif(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "red",linetype=2)+
          annotate("segment", x = BSVar1$Tm[1],
                   xend = BSVar1$Tm[1], y = 0, yend = 0.5,
                   colour = "red",linetype=2)+
          ggplot2::ggtitle(paste0(as.character(df1[1])," ",str_replace(df2$sample_name[1],"S",paste0("\u03A6"))))+
          ylim(-0.75,0.5+max(c(max(BSVar$I,na.rm=TRUE),max(BSVar1$I,na.rm=TRUE))))+
          xlim(37,68)+
          theme(legend.position="bottom")+
          ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]),size=3.5)+
          ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",signif(pred1$Tm[1],3),"\u00B0C"),size=3.5)+
          annotate("text",
                   x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = -0.15,
                   label=paste0(signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3)),
                   colour="blue",
                   size=3.5
          )+
          annotate("segment", x = min(fitted.values$C), xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   y = 0.5, yend = 0.5,
                   colour = "blue",linetype=2)+
          annotate("segment", x = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3),
                   xend = signif(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,3), y = 0, yend = 0.5,
                   colour = "blue",linetype=2)
        return(plot)
      }
    }
  }
}

spSim<-function(df2,species=9606,threshold=500,string_version=11.5){#df2 would be the data 
  
  #set C and I as numeric
  df2$temperature<-as.numeric(as.vector(df2$temperature))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-df2$Accession
  df1<-df2
  df1$uniqueID<-df1$Accession
  df2$uniqueID<-as.character(df2$uniqueID)
  
  ###########################################
  df1<-as.character(df1$uniqueID)
  
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
  
  BSVarN<-df2[df2$uniqueID %in% df1,]
  BSVar1 <-df2[df2$uniqueID %in% df1 & df2$treatment== "treated",]
  BSVar <-df2[df2$uniqueID %in% df1 & df2$treatment== "vehicle",]
  
  BSVarN<-BSVarN %>% dplyr::mutate(treatment="null") 
  BSVar<-BSVar %>% dplyr::mutate(treatment="vehicle")
  BSVar1<-BSVar1 %>% dplyr::mutate(treatment="treated")
  BSVarN$treatment<-as.factor(BSVarN$treatment)
  BSVar$treatment<-as.factor(BSVar$treatment)
  BSVar1$treatment<-as.factor(BSVar1$treatment)
  
  BSVar<-BSVar[!is.na(BSVar$I),]
  BSVar1<-BSVar1[!is.na(BSVar1$I),]
  BSVarN<-BSVarN[!is.na(BSVarN$I),]
  
  
  #load stringdb data
  if(species==9606){
    string_db <- STRINGdb$new( version="11.5", species=9606,score_threshold=threshold, input_directory="")
    Known_proteins<-string_db$mp(c("MAP2K1","MAP2K2"))
    
  }else if(species==7955){
    string_db <- STRINGdb$new( version="11.5", species=7955,score_threshold=treshold, input_directory="")
    
    Known_proteins<-string_db$mp(c("Q7ZUY3",#histone
                                   "Q6NV46",# stat3, 
                                   "Q68SP3",# stat5a, 
                                   "A4QNT9",# stat5b, 
                                   "Q90Y03",# aldh1a2,
                                   "Q0H2G3",# aldh1a3,
                                   "F1QZU7",# aldh1a2.2, 
                                   "A2BGR9",# aldh1a2.1, 
                                   "F1QLV5",#nqo1,
                                   "F1Q7F3"))
  }
  #CD81<-string_db$mp("CD81")
  #string_db$get_neighbors( c(MEK1, MEK2,CD81) ) [1:10]
  proteins<-string_db$get_proteins()
  graph<-string_db$get_graph()
  
  #get interactors and ensembl ids
  Interactors<-string_db$get_interactions( Known_proteins )
  
  if(length(Interactors)>5000){
    TP<-c(Known_proteins,Interactors,string_db$get_neighbors( Known_proteins ))[1:5000]
  }else{
    TP<-c(Known_proteins,Interactors,string_db$get_neighbors( Known_proteins ))
  }
  TP<-unlist(TP)
  #get protein data for interactors
  xx<-proteins[proteins$protein_external_id %in% TP,]
  symbols <- xx$preferred_name
  ext_id<-xx$protein_external_id
  #get uniprotID
  if(species=="9606"){
    ss<-AnnotationDbi::mapIds(org.Hs.eg.db,symbols,"UNIPROT",'SYMBOL')
    ss_names<-AnnotationDbi::mapIds(org.Hs.eg.db,symbols,"GENENAME",'SYMBOL')
  }else if (species ==7955){
    ss <- as.character(AnnotationDbi::mapIds(org.Dr.eg.db, symbols, "UNIPROT", 'SYMBOL'))
    ss_names<-AnnotationDbi::mapIds(org.Hs.eg.db,symbols,"GENENAME",'SYMBOL')
  }
  
  ss<-ss[!is.na(ss)]#remove missing values
  #plot ppi network
  example1_mapped <- string_db$map( data.frame(gene=ss), "gene", removeUnmappedRows = TRUE )
  #remove unmapped identifiers
  
  #hits
  if(nrow(example1_mapped)>1000){
    hits <- example1_mapped$STRING_id[1:1000]
  }else{
    hits <- example1_mapped$STRING_id
  }
  hits<-hits[!is.na(hits)]
  #enrich<-STRINGdb::STRINGdb$get_ppi_enrichment_full(hits,graph)
  
  # # see how many proteins do you have    
  # vcount(graph)
  # 
  # # find top 200 proteins with the highest degree
  # top.degree.verticies <- names(tail(sort(degree(graph)), 200))
  # 
  # # extract the relevant subgraph
  # top.subgraph <- induced_subgraph(graph, top.degree.verticies)
  # 
  # # count the number of proteins in it
  # vcount(top.subgraph)
  if(species==9606){
    ss<-unlist(c(ss,list("P36507","Q02750")))
  }
  
  BSVar_<-dplyr::bind_rows(BSVar) %>%
    dplyr::filter(uniqueID %in% ss) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(confidence=as.factor(ifelse(threshold>700,"high","low")))
  
  BSVar1_<-dplyr::bind_rows(BSVar1) %>%
    dplyr::filter(uniqueID %in% ss) %>%
    dplyr::mutate(confidence=as.factor(ifelse(threshold>700,"high","low"))) %>%
    dplyr::ungroup(.) 
  
  BSVarN_<-dplyr::bind_rows(BSVarN) %>%
    dplyr::filter(uniqueID %in% ss) %>%
    dplyr::mutate(confidence=as.factor(ifelse(threshold>700,"high","low"))) %>%
    dplyr::ungroup(.) 
  
  CID<-dplyr::intersect(unique(dplyr::bind_rows(BSVar_)$uniqueID),unique(dplyr::bind_rows(BSVar1_)$uniqueID))
  BSVar_<-BSVar_ %>% dplyr::filter(uniqueID %in% CID)%>%
    dplyr::group_split(uniqueID,sample_name)
  BSVar1_<-BSVar1_ %>% dplyr::filter(uniqueID %in% CID)%>%
    dplyr::group_split(uniqueID,sample_name)
  BSVarN_<-BSVarN_ %>% dplyr::filter(uniqueID %in% CID)%>%
    dplyr::group_split(uniqueID,sample_name)
  
  
  
  
  #add some noise to the original data
  set.seed(1)
  y_data <- purrr::map(BSVar_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$temperature), 0, 0.05)))
  y_data1 <- purrr::map(BSVar1_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$temperature), 0, 0.05)))
  y_dataN<-purrr::map(BSVarN_,function(x) x %>% dplyr::mutate(I=x$I + rnorm(length(x$temperature), 0, 0.05)))
  #show simulation results with TP
  fT<-TRUE
  show_results<-TRUE
  #show results with some noise
  spresults<-spstat(y_dataN,y_data,y_data1,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  TP<-spresults %>% dplyr::mutate(outcome="TP")
  #set false positives by overlaying vehicle curves
  spresults<-spstat(y_dataN,y_data,y_data,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  FP<-spresults %>% dplyr::mutate(outcome="FP")
  test<-rbind(TP,FP)
  
  if(any(names(test)=="dRSS")){
    test<-test %>% dplyr::rename("rssDiff"='dRSS')
  }
  # 
  # roc4 <- roc(test$outcome,
  #             test$Fvals, percent=TRUE,
  #             # arguments for ci
  #             ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
  #             # arguments for plot
  #             plot=TRUE,
  #             print.auc=TRUE, show.thres=TRUE)+title(test$sample_name[1])
  # 
  rocs <- pROC::roc(outcome ~Fvals + rssDiff + dTm+AUC,data = test)
  rocs$Fvals$auc[1]<-round(rocs$Fvals$auc[1],3)
  rocs$rssDiff$auc[1]<-round(rocs$rssDiff$auc[1],3)
  rocs$dTm$auc[1]<-round(rocs$dTm$auc[1],3)
  rocs$AUC$auc[1]<-round(rocs$AUC$auc[1],3)
  ROC<- ggroc(rocs)+ggtitle(paste0(str_replace(test$sample_name[1],"S",paste0("\u03A6"))," simulated data"))+
    scale_color_colorblind("Parameters",labels=c("F-stat",
                                                 expression(Delta*RSS),
                                                 expression(Delta*Tm),
                                                 expression(Delta*AUC)))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    ggplot2::annotate("text", x=0.25, y=0.25, label= paste("F-stat AUC = ", rocs$Fvals$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.15, label= paste("RSS AUC = ", rocs$rssDiff$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.05, label= paste("Tm AUC = ", rocs$dTm$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=-0.05, label= paste("AUC = ", rocs$AUC$auc[1]),size=3.5)
  
  
  
  
  #show results without
  spresults<-spstat(BSVarN_,BSVar_,BSVar1_,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  TP<-spresults %>% dplyr::mutate(outcome="TP")
  #set false positives by overlaying vehicle curves
  spresults<-spstat(BSVarN_,BSVar_,BSVar_,Ftest=fT,show_results=show_results,filters=FALSE,scaled_dof=FALSE,Peptide=FALSE)
  FP<-spresults %>% dplyr::mutate(outcome="FP")
  test<-rbind(TP,FP)
  if(any(names(test)=="dRSS")){
    test<-test %>% dplyr::rename("rssDiff"='dRSS')
  }
  
  rocs1<- roc(outcome ~Fvals + rssDiff + dTm+AUC,data = test)
  
  rocs1$Fvals$auc[1]<-round(rocs1$Fvals$auc[1],3)
  rocs1$rssDiff$auc[1]<-round(rocs1$rssDiff$auc[1],3)
  rocs1$dTm$auc[1]<-round(rocs1$dTm$auc[1],3)
  rocs1$AUC$auc[1]<-round(rocs1$AUC$auc[1],3)
  ROC1<- ggroc(rocs1)+ggtitle(paste0(str_replace(test$sample_name[1],"S",paste0("\u03A6"))," real data"))+
    scale_color_colorblind("Parameters",labels=c("F-stat",
                                                 expression(Delta*RSS),
                                                 expression(Delta*Tm),
                                                 expression(Delta*AUC)))+
    theme(legend.position="bottom", legend.box = "horizontal")+
    ggplot2::annotate("text", x=0.25, y=0.25, label= paste("F-stat AUC = ", rocs1$Fvals$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.15, label= paste("RSS AUC = ", rocs1$rssDiff$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=0.05, label= paste("Tm AUC = ", rocs1$dTm$auc[1]),size=3.5)+
    ggplot2::annotate("text", x=0.25, y=-0.05, label= paste("AUC = ", rocs1$AUC$auc[1]),size=3.5)
  
  ROC_list<-list(ROC,ROC1)
  
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

sigC<-function(df_,Protein,Peptide=FALSE,stats=FALSE){
  Protein<-as.character(Protein)
  df_$C<-as.numeric(as.character(df_$C))
  if(isTRUE(stats)){
    DFN<-df_
    df_1<-df_%>% dplyr::filter(treatment=="treated")
    df_<-df_%>% dplyr::filter(treatment=="vehicle")
  }else{
    DFN<-df_%>% dplyr::filter(uniqueID %in% Protein)
    df_1<-df_%>% dplyr::filter(uniqueID%in%Protein,treatment=="treated")
    df_<-df_%>% dplyr::filter(uniqueID%in%Protein,treatment=="vehicle")
  }
  
  if(isTRUE(Peptide)){
    df_<-df_ %>% dplyr::rename("value"="I3","temperature"="C")
    df_1<-df_1 %>% dplyr::rename("value"="I3","temperature"="C")
    DFN<-DFN %>% dplyr::rename("value"="I3","temperature"="C")
  }else{
    df_<-df_ %>% dplyr::rename("value"="I","temperature"="C")
    df_1<-df_1 %>% dplyr::rename("value"="I","temperature"="C")
    DFN<-DFN %>% dplyr::rename("value"="I","temperature"="C")
  }
  
  if(isTRUE(stats)){
    nlm2<-function(x){
      nlm2<-x %>% 
        dplyr::mutate(fit=list(cetsa_fit2(.,norm=FALSE)),
                      Tm=as.numeric(coef(fit[[1]][[1]])[3]))
      
      return(nlm2) 
    }
    df_<-df_ %>% dplyr::group_split(uniqueID)
    df_1<-df_1 %>% dplyr::group_split(uniqueID)
    DFN<-DFN %>% dplyr::group_split(uniqueID)
    if (.Platform$OS.type=="windows"){
      nlm2<-parallel::mclapply(df_,nlm2)
    }else{
      nlm2<-parallel::mclapply(df_,nlm2,mc.cores=future::availableCores())
    }
  }else{
    nlm2<-df_  %>%   
      dplyr::mutate(fit=list(cetsa_fit2(.,norm=FALSE)))
    nlm2$Tm<-NA
    nlm2$Tm<-as.numeric(coef(nlm2$fit[[1]][[1]])[3])
  }
  if(class(nlm2$fit[[1]][[1]])=='try-error'){
    warning("the function could not fit the data")
  }
  
  
  
  #ready confidence intervals for plot
  result<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(.$fit[[1]][[2]]$lwr.conf),
                                                        HI = list(.$fit[[1]][[2]]$upr.conf),
                                                        CP=list(.$fit[[1]][[2]]$temperature),
                                                        IP=list(.$fit[[1]][[2]]$value),
                                                        dof=summary(.$fit[[1]][[1]])$df[2])# lower CI
  result<-result %>% dplyr::rename("I"="value","C"="temperature")
  
  result <-result %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]][[1]]))
  
  
  nlm1<-list()
  CT<-list()
  dfc<-list()
  
  #sigmoidal fit for treated
  nlm2<-df_1  %>%   
    dplyr::mutate(fit=list(try(cetsa_fit2(.,norm=FALSE))))
  nlm2$Tm<-NA
  nlm2$Tm<-as.numeric(coef(nlm2$fit[[1]][[1]])[3])
  
  if(class(nlm2$fit[[1]][[1]])=='try-error'){
    warning("the function could not fit the data")
  }
  
  #ready confidence intervals for plot
  result1<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(.$fit[[1]][[2]]$lwr.conf),
                                                         HI = list(.$fit[[1]][[2]]$upr.conf),
                                                         CP=list(.$fit[[1]][[2]]$temperature),
                                                         IP=list(.$fit[[1]][[2]]$value),
                                                         dof=summary(.$fit[[1]][[1]])$df[2])# lower CI
  result1<-result1 %>% dplyr::rename("I"="value","C"="temperature")
  
  result1 <-result1 %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]][[1]]))
  
  
  
  return(rbind(Pred,Pred1))
}

sigfit<-function(SigF,Peptide=FALSE){
  
  #
  if(isTRUE(Peptide)){
    
    Pred<-SigF %>%
      subset(treatment=="vehicle") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,
                    LOW,HI,sample_name,sample_id,missing_pct)
    
    Pred1<-SigF%>%
      subset(treatment=="treated") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,
                    LOW,HI,sample_name,sample_id,missing_pct)
    
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    
    #unnest data
    Pred<-Pred %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred1<-Pred1 %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred$AUC<-round(pracma::trapz(Pred$IP),2)
    Pred1$AUC<-round(pracma::trapz(Pred1$IP),2)
    Pred1$dAUC<-abs(Pred1$AUC[1]-Pred$AUC[1])
    
    Pred1$dAUC<-pracma::trapz(Pred1$IP[(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))-1):(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))+1)])-pracma::trapz(Pred$IP[(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))-1):(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))+1)])
    Pred1$dAUC<-abs(round(Pred1$dAUC[1],3))
    
    #max replicates
    roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
      if(length(x) != 1) stop("'x' must be of length 1")
      10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
    }
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    
    #append missing value data
    if(isTRUE(Peptide)){
      
      Pred<-Pred %>%
        distinct(.) %>%
        dplyr::group_by(uniqueID,treatment,sample_id,Annotated_Sequence) %>%
        dplyr::group_split(.)
      Pred<-purrr::map(Pred,function(x) x %>% dplyr::mutate(missing_v=100*(roundUpNice(length(unique(x$I)))-length(unique(x$I)))/roundUpNice(length(unique(x$I)))))
      Pred<-dplyr::bind_rows(Pred)
      Pred<-Pred %>% dplyr::group_split(uniqueID,treatment)
      Pred<-purrr::map(Pred,function(x) x %>%
                         dplyr::mutate(missing_v=getmode(x$missing_v)))
      
      Pred1<-Pred1 %>%
        distinct(.) %>%
        dplyr::group_by(uniqueID,treatment,sample_id,Annotated_Sequence) %>%
        dplyr::group_split(.)
      Pred1<-purrr::map(Pred1,function(x) x %>%
                          dplyr::mutate(missing_t=100*(roundUpNice(length(unique(x$I)))-length(unique(x$I)))/roundUpNice(length(unique(x$I)))))
      Pred1<-dplyr::bind_rows(Pred1)
      Pred1<-Pred1 %>%
        dplyr::group_split(uniqueID,treatment)
      Pred1<-purrr::map(Pred1,function(x) x %>%
                          dplyr::mutate(missing_t=getmode(x$missing_t)))
      
      
      Pred<-dplyr::bind_rows(Pred)
      Pred1<-dplyr::bind_rows(Pred1)
    }else{
      Pred$missing_v<-round(Pred$missing_v[1],0)
      Pred1$missing_t<-round(Pred1$missing_t[1],0)
    }
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x=C,y=IP,color=treatment))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred1,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      ggplot2::annotate("text", x=45, y=-0.65, label= paste("missing: ",round(Pred$missing_v[1],0),"%"),colour="#00BFC4")+ 
      ggplot2::annotate("text", x=45, y=-0.75, label= paste("missing: ",round(Pred1$missing_t[1],0),"%"),colour="#F8766D")+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = 37, xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +
      ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = 37, xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2)+ylim(-0.8,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
    print(P1)
  }else{
    Pred<-SigF %>%
      subset(treatment=="vehicle") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,Cell_Component,Bio_Process,Coverage,MW_kDa,
                    LOW,HI,sample_name,sample_id,missing_pct,rank)
    Pred1<-SigF%>%
      subset(treatment=="treated") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,Cell_Component,Bio_Process,Coverage,MW_kDa,
                    LOW,HI,sample_name,sample_id,missing_pct,rank)
    
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #unnest data
    Pred<-Pred %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred1<-Pred1 %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    
    Pred$AUC<-round(pracma::trapz(Pred$IP),2)
    Pred1$AUC<-round(pracma::trapz(Pred1$IP),2)
    Pred1$dAUC<-abs(Pred1$AUC[1]-Pred$AUC[1])
    
    
    Pred1$dAUC<-pracma::trapz(Pred1$IP[(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))-1):(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))+1)])-pracma::trapz(Pred$IP[(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))-1):(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))+1)])
    Pred1$dAUC<-abs(round(Pred1$dAUC[1],3))
    
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x=C,y=IP,color=treatment))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred1,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      ggplot2::annotate("text", x=45, y=-0.65, label= paste("missing: ",round(Pred$missing_v[1],0),"%"),colour="#00BFC4")+ 
      ggplot2::annotate("text", x=45, y=-0.75, label= paste("missing: ",round(Pred1$missing_t[1],0),"%"),colour="#F8766D")+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = 37, xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +
      ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = 37, xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2)+ylim(-0.8,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
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
    test<-test%>% dplyr::select(treatment,CV_pct,treatment,C,sample_name) %>% unique(.)
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
    test$treatment<-as.factor(test$treatment)
    
    ggplot(test,aes(y=CV_pct,x=treatment,fill=treatment))+
      facet_grid(c(~temperature,~carrier),scales="free_x")+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample_id")+
      theme(axis.title.x = element_text(face="bold",size="14",colour="white"),
            axis.title.y = element_text(face="bold",size="14",colour="black"),
            axis.text.x = element_text(angle = 90,face="bold",size="14",colour="black"),
            axis.text.y = element_text(face="bold",size="14",colour="black"),
            legend.text = element_text(face="bold",size="14",colour="black"),
            legend.title = element_text(face="bold",size="14",colour="black"),
            strip.text.x = element_text(size = 14, colour = "black"))+
      ggplot2::ylim(0,200)
  }else{  
    
    test$treatment<-test$treatment
    ggplot(test,aes(y=CV_pct,x=treatment,fill=treatment))+
      facet_grid(~temperature)+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA,aes(alpha=0.2))+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("RSD%")+
      ggplot2::xlab("sample_id")+
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
    test<-test%>% dplyr::select(treatment,CV_pct,treatment,C,sample_name) %>% unique(.)
    test<-test %>% dplyr::rename("temperature"="C")
    test<-test %>% dplyr::left_join(df.temps,by="temperature")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum.File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$treatment<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
    test<-test %>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temp_ref","S_N"="Average_Reporter_S/N","PEP"="Percolator_PEP",
                                 "MissedCleavages"="#_MissedCleavages","DeltaM"="DeltaM_[ppm]","IonInjTime"="Ion_Inject_Time_[ms]",
                                 "I_Interference"="Isolation_Interference_[%]")
    
    
  }else{
    test<-df_raw
    test<-test %>% dplyr::left_join(df.temps,by="temp_ref")
    test <-test %>% dplyr::mutate(CC=ifelse(stringr::str_detect(Spectrum.File,"DMSO")==TRUE,0,1))#concentration values are defined in uM
    
    test$treatment<-ifelse(test$CC==0,"vehicle","treated")
    
    
    test$sample_name<-paste0(ifelse(str_detect(test$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(str_detect(test$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(test$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                             ifelse(str_detect(test$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(str_detect(test$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
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
MSStats_converter<-function(df_raw,solvent,ref){
  editPSMs2<-data.frame()
  ref<-as.character(ref)
  editPSMs2<-df_raw %>% dplyr::rename("ProteinName"="Accession",
                                      "PeptideSequence"="Annotated_Sequence",
                                      "Run"="Spectrum.File",
                                      #"PSM"="PSMs_Peptide_ID",
                                      "Channel"="temp_ref",
                                      "Intensity"="value")
  #Condition, Bioreplicate and TechRepMixture need to be filled
  
  editPSMs2$Condition<-ifelse(editPSMs2$Channel==ref,"Norm",0)
  if(nchar(editPSMs2$Run[1])>8){
    editPSMs2$BioReplicate<-ifelse(stringr::str_detect(editPSMs2$Run,solvent)=="TRUE","vehicle","treated")
    editPSMs2<-editPSMs2 %>% 
      dplyr::mutate(Mixture=paste0(ifelse(stringr::str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))
  }else{
    editPSMs2$BioReplicate<-ifelse(stringr::str_detect(editPSMs2$Run,"01"),1,2)
    editPSMs2<-editPSMs2 %>% 
      dplyr::mutate(Mixture=paste0(ifelse(stringr::str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                   ifelse(stringr::str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))
    
  }
  editPSMs2$TechRepMixture<-1
  editPSMs2<-editPSMs2 %>% dplyr::select(ProteinName,PeptideSequence,Charge,PSM,Mixture,TechRepMixture,Run,Channel,Condition,BioReplicate,Intensity)
  
  Annotation<-editPSMs2 %>% dplyr::select(Run,TechRepMixture,Channel,Condition,Mixture,BioReplicate)
  Annotation$Fraction<-1
  Annotation<-Annotation %>% dplyr::group_split(Run)
}



#TPP TR Reader
TPPbenchmark<-function(f,overlaps=NA,volcano=TRUE,filters="TPP"){
  f<-f
  f<-list.files(f)
  #read data
  df_TPP<-lapply(f,function(x) read_excel(x,.name_repair = "unique"))
  #extract experiment names
  f<-str_extract_all(f,c("C_F_E","C_F_S","C_nF_E","C_nF_S","nC_F_E","nC_F_S","nC_nF_E","nC_nF_S"))
  f<-lapply(f,function(x) data.frame(sample_name=as.factor(x)))
  
  #join data
  df_TPP<-purrr::map2(df_TPP,f,function(x,y)cbind(x,y))
  
  
  #select columns of interest
  if(filters=="HQ"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(!is.na(meltPoint_MEKi_2) & !is.na(meltPoint_MEKi_2) & !is.na(meltPoint_DMSO_2) &!is.na(meltPoint_DMSO_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_MEKi_2 >0.8 & R_sq_MEKi_1 >0.8 & R_sq_DMSO_2>0.8 &R_sq_DMSO_1>0.8 &
                      plateau_MEKi_2<0.3 & plateau_MEKi_1<0.3 & plateau_DMSO_2<0.3 & plateau_DMSO_1<0.3)
  }else if(filters=="Sigmoidal"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    model_converged_MEKi_2,model_converged_MEKi_1,model_converged_DMSO_2,model_converged_DMSO_1,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(model_converged_DMSO_1=="Yes",model_converged_DMSO_2=="Yes",model_converged_MEKi_1=="Yes",model_converged_MEKi_2=="Yes")
  }else{
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_MEKi_2,R_sq_MEKi_1,R_sq_DMSO_2,R_sq_DMSO_1,
                    plateau_MEKi_2,plateau_MEKi_1,plateau_DMSO_2,plateau_DMSO_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_MEKi_2,meltPoint_MEKi_1,meltPoint_DMSO_2,meltPoint_DMSO_1,
                    sample_name,
                    pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2)%>%
      dplyr::filter(!is.na(meltPoint_MEKi_2) & !is.na(meltPoint_MEKi_2) & !is.na(meltPoint_DMSO_2) &!is.na(meltPoint_DMSO_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_MEKi_2 >0.8 & R_sq_MEKi_1 >0.8 & R_sq_DMSO_2>0.8 &R_sq_DMSO_1>0.8 &
                      plateau_MEKi_2<0.3 & plateau_MEKi_1<0.3 & plateau_DMSO_2<0.3 & plateau_DMSO_1<0.3 & fulfills_all_4_requirements=="Yes")
  }
  df_TPP$Protein_ID<-as.factor(df_TPP$Protein_ID)
  df_TPP$sample_name<-str_replace(df_TPP$sample_name,"S","\u03A6")
  #get names
  
  df_TPP<-df_TPP %>% dplyr::group_split(sample_name) 
  df_TPP1<-purrr::map(df_TPP,function(x)x %>% dplyr::select(Protein_ID,sample_name,pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2) %>% 
                        pivot_longer(c(pVal_adj_MEKi_1_vs_DMSO_1,pVal_adj_MEKi_2_vs_DMSO_2),
                                     names_to = c("hi","dTm"),
                                     names_pattern = c("(.+)pVal_(.+)"),
                        ) %>% dplyr::rename("p_dTm"="value"))
  df_TPP1<-purrr::map(df_TPP1,function(x)x %>% dplyr::select(Protein_ID,sample_name,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2,p_dTm) %>% 
                        pivot_longer(c(diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_MEKi_2_vs_DMSO_2),
                                     names_to = c("hi","set"),
                                     names_pattern = "(.+)diff_(.+)"
                        ) %>% dplyr::rename("dTm"="value") %>% dplyr::select(-hi,-set))
  df_TPP1<-purrr::map(df_TPP1,function(x) x[!is.na(x$dTm),])
  if(isTRUE(volcano)){
    
    df_TPP2<-dplyr::bind_rows(df_TPP1) 
    df_TPP2$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_TPP2$diffexpressed[df_TPP2$dTm > 2 & df_TPP2$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_TPP2$diffexpressed[df_TPP2$dTm < -2 & df_TPP2$p_dTm < 0.05] <- "Destabilized"
    df_TPP2$delabel <- NA
    df_TPP2$delabel[df_TPP2$Protein_ID %in%c("P36507","Q02750")] <- as.character(df_TPP2$Protein_ID[df_TPP2$Protein_ID %in%c("P36507","Q02750")])
    #get overlaps in data
    #df_TPP2$Stroke <- NA
    #df_TPP2$Stroke<-ifelse(df_TPP2$diffexpressed=="Stabilized" | df_TPP2$diffexpressed=="Destabilized",ifelse(df_TPP2$Protein_ID %in% overlaps,TRUE,FALSE),FALSE)
    
    
    
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% distinct(.) %>%  dplyr::group_split(sample_name,Protein_ID)
    df_TPP2<-purrr::map(df_TPP2,function(x) x[1,])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name)
    flevels<-data.frame(colors=c("green","black","red"))
    df_TPP3<-purrr::map(df_TPP2,function(x)
      ggplot(data=x,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+ 
        labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
        xlim(-15,15)+
        theme(legend.position="bottom", legend.box = "horizontal")+
        geom_label(aes(16,y=40),label=paste0("S = ",nrow(x[x$diffexpressed=="Stabilized",])),show.legend=FALSE)+
        geom_label(aes(-16,y=40),label=paste0("DS = ",nrow(x[x$diffexpressed=="Destabilized",])),show.legend=FALSE)+
        scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
        #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
        geom_label_repel(aes(dTm, -log10(p_dTm),label=delabel),color="black",  nudge_x = 1,
                         force = 3,
                         box.padding = 3,
                         segment.alpha = .5)+
        ggtitle(x$sample_name[1])+ylim(-0.1,55))#+
    #geom_point(data=df_TPP2[df_TPP2$Stroke==1,],
    #pch=21, fill=NA, size=4, colour="black", stroke=1)
    
    
    
    return(df_TPP3)
  }
  
  #for number of curves
  
  df_1<-dplyr::bind_rows(df_TPP1)
  
  df_1$sample_name<-str_replace(df_1$sample_name,"S","\u03A6")
  df_1<-dplyr::bind_rows(df_1) 
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-df_1 %>% dplyr::select(p_dTm,uniqueID,sample_name,dTm) %>% distinct(.)
  df_1<-df_1[!duplicated(df_1),]
  df_1<-df_1[!is.na(df_1$dTm),]
  df_1$p_dTm<-base::round(df_1$p_dTm,3) 
  df_1<-dplyr::bind_rows(df_1) %>% distinct(.)%>% dplyr::group_split(uniqueID,sample_name)
  df_1<-purrr::map(df_1,function(x) x[1,])
  df_1<-dplyr::bind_rows(df_1)
  
  colors<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_1)$sample_name)))
  colors$sample_name<-as.character(colors$sample_name)
  colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
  
  check1<-df_1 %>% count(sample_name) %>%
    mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
    ggplot() +
    ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
    ggforce::theme_no_axes()+
    scale_fill_manual(values = colors$hex,aesthetics="fill")+
    xlim(-1.1,1.45)+
    ggplot2::geom_label(mapping=aes(x=colors$x,y=colors$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
    ggtitle("Number of fitted curves")
  if(!isTRUE(filters)){
    pdf("Number_of_fitted_curves_TPP_Protein_only_converging.pdf",encoding="CP1253.enc",compress=TRUE,width=6.12,height=4.02)
    check1
    dev.off()
  }else{
    pdf("Number_of_fitted_curves_TPP_iMAATSA_4_req_Filters.pdf",encoding="CP1253.enc",compress=TRUE,width=6.12,height=4.02)
    check1
    dev.off()
  }
  # 
  # 'Size'=(
  #   intersection_size(
  #     mode=mode,
  #     mapping=aes(fill=exclusive_intersection),
  #     size=0,
  #     text=list(check_overlap=TRUE)
  #   ) + scale_fill_venn_mix(
  #     data=abc_data,
  #     guide=FALSE,
  #     colors=c('A'='red', 'B'='blue', 'C'='green3')
  # 
  df_1<-dplyr::bind_rows(df_TPP)
  df_1$uniqueID<-df_1$Protein_ID
  df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
    pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
  
  
  
  
  #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
  # df_ <- df_ %>%
  #   mutate_if(sapply(df_, is.factor), as.numeric)
  IDs<-df_1$uniqueID
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
  df_1<-cbind(IDs,df_1)
  colors<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
  colors$sample_name<-as.character(colors$sample_name)
  colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
  rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                   values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
  df_1<-na.omit(df_1)
  check<-list()
  check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID")],
               #min_degree=6,
               set_sizes=FALSE,
               n_intersections=10,
               min_degree=1,
               encode_sets=FALSE,
               stripes='white',
               # intersections=list(
               #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
               #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
               #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
               #   
               # ),
               # mode = 'inclusive_union',
               #mode=mode,
               # intersections=list(
               #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
               #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
               #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
               #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
               #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
               #   
               # ),
               #sort_intersections_by="cardinality",
               base_annotations=list(
                 '# of fitted curves'=
                   intersection_size(
                     # mode='inclusive_intersection',
                     counts=TRUE,
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   #scale_color_brewer(palette="Dark2")+
                   theme(legend.position = 'none',
                         # axis.text=element_text(size=20, face='bold'),
                         # axis.title=element_text(size=20, face='bold')
                   )+
                   ggplot2::ylab("# of fitted curves")
               ),               # base_annotations=list(
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
               height_ratio=0.8
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) ",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP" | filters=="HQ","filtered","unfiltered")))
  
  #check upset plots for intersections
  check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
  #filter exclusive intersection colors
  colors<-colors %>% dplyr::filter(sample_name %in% check_intersections$IDs)
  
  #add query
  check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID")],
               #min_degree=6,
               set_sizes=FALSE,
               n_intersections=10,
               min_degree=1,
               encode_sets=FALSE,
               stripes='white',
               # intersections=list(
               #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
               #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
               #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
               #   
               # ),
               # mode = 'inclusive_union',
               #mode=mode,
               # intersections=list(
               #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
               #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
               #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
               #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
               #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
               #   
               # ),
               #sort_intersections_by="cardinality",
               base_annotations=list(
                 '# of fitted curves'=
                   intersection_size(
                     # mode='inclusive_intersection',
                     counts=TRUE
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   ggplot2::ylab("# of fitted curves")
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
               queries<-list(
                 upset_query(
                   intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
                   color='black',
                   fill='black',
                   only_components='# of fitted curves'
                 ),
                 upset_query(
                   intersect=colors$sample_name[1],
                   color=colors$hex[1],
                   fill=colors$hex[1],
                   only_components='# of fitted curves'
                 ),
                 upset_query(
                   intersect=colors$sample_name[2],
                   color=colors$hex[2],
                   fill=colors$hex[2],
                   only_components='# of fitted curves'
                 )
                 # ,
                 # upset_query(
                 #   intersect=colors$sample_name[3],
                 #   color=colors$hex[3],
                 #   fill=colors$hex[3],
                 #   only_components='# of fitted curves'
                 # )
                 # ,
                 # upset_query(
                 #   intersect=colors$sample_name[4],
                 #   color=colors$hex[4],
                 #   fill=colors$hex[4],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect=colors$sample_name[5],
                 #   color=colors$hex[5],
                 #   fill=colors$hex[5],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect=colors$sample_name[6],
                 #   color=colors$hex[6],
                 #   fill=colors$hex[6],
                 #   only_components='# of fitted curves'
                 # ),
                 # upset_query(
                 #   intersect= as.character(colors$sample_name[7]),
                 #   color=as.character(colors$hex[7]),
                 #   fill=as.character(colors$hex[7]),
                 #   only_components='# of fitted curves'
                 # )
                 # ,
                 # upset_query(
                 #   intersect=as.character(colors$sample_name[8]),
                 #   color=as.character(colors$hex[8]),
                 #   fill=as.character(colors$hex[8]),
                 #   only_components='# of fitted curves'
                 # )
               )
               
               
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) (",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP","additional filters)",ifelse(filters=="HQ","high-quality filtered)","unfiltered)"))))
  if(!isTRUE(filters)){
    pdf("Upset_fitted_curves_TPP_Protein_only_converging.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
    check
    dev.off()
  }else{
    pdf("Upset_fitted_curves_TPP_iMAATSA_4_req_additional_Filters.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
    check
    dev.off()
  }
  return(list(check,check1))
  
  
}
TPPbenchmark_generic<-function(f,overlaps=NA,volcano=TRUE,filters="TPP",Peptide=FALSE,filter=TRUE){
  f<-f
  f<-list.files(f)
  #read data
  df_TPP<-lapply(f,function(x) read_excel(x,.name_repair = "unique"))
  #extract experiment names
  f<-str_extract_all(f,c("C_F_E","C_F_S","C_nF_E","C_nF_S","nC_F_E","nC_F_S","nC_nF_E","nC_nF_S"))
  
  f<-lapply(f,function(x) data.frame(sample_name=as.factor(x)))
  f<-f %>% purrr::keep(function(x) nrow(x)>0)
  #join data
  df_TPP<-purrr::map2(df_TPP,f,function(x,y)cbind(x,y))
  
  
  #select columns of interest
  if(filters=="HQ"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Vehicle_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(!is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) &!is.na(meltPoint_Vehicle_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_Treatment_2 >0.8 & R_sq_Treatment_1 >0.8 & R_sq_Treatment_2>0.8 &R_sq_Vehicle_1>0.8 &
                      plateau_Treatment_2<0.3 & plateau_Treatment_1<0.3 & plateau_Treatment_2<0.3 & plateau_Vehicle_1<0.3)
    
  }else if(filters=="Sigmoidal"){
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    model_converged_Treatment_2,model_converged_Treatment_1,model_converged_Treatment_2,model_converged_Vehicle_1,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Treatment_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(model_converged_Vehicle_1=="Yes",model_converged_Treatment_2=="Yes",model_converged_Treatment_1=="Yes",model_converged_Treatment_2=="Yes")
    
  }else{
    df_TPP<-dplyr::bind_rows(df_TPP) %>% 
      dplyr::select(Protein_ID,fulfills_all_4_requirements,
                    minSlopes_less_than_0.06,
                    R_sq_Treatment_2,R_sq_Treatment_1,R_sq_Treatment_2,R_sq_Vehicle_1,
                    plateau_Treatment_2,plateau_Treatment_1,plateau_Treatment_2,plateau_Vehicle_1,
                    diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,
                    meltP_diffs_have_same_sign,
                    meltPoint_Treatment_2,meltPoint_Treatment_1,meltPoint_Treatment_2,meltPoint_Vehicle_1,
                    sample_name,
                    pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2)%>%
      dplyr::filter(!is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) & !is.na(meltPoint_Treatment_2) &!is.na(meltPoint_Vehicle_1) &
                      minSlopes_less_than_0.06=="Yes"
                    & R_sq_Treatment_2 >0.8 & R_sq_Treatment_1 >0.8 & R_sq_Treatment_2>0.8 &R_sq_Vehicle_1>0.8 &
                      plateau_Treatment_2<0.3 & plateau_Treatment_1<0.3 & plateau_Treatment_2<0.3 & plateau_Vehicle_1<0.3 & fulfills_all_4_requirements=="Yes")
  }
  df_TPP$Protein_ID<-as.factor(df_TPP$Protein_ID)
  df_TPP$sample_name<-str_replace(df_TPP$sample_name,"S","\u03A6")
  #get names
  
  df_TPP<-df_TPP %>% dplyr::group_split(sample_name) 
  df_TPP1<-purrr::map(df_TPP,function(x)x %>% dplyr::select(Protein_ID,sample_name,pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2) %>% 
                        pivot_longer(c(pVal_adj_Treatment_1_vs_Vehicle_1,pVal_adj_Treatment_2_vs_Vehicle_2),
                                     names_to = c("hi","dTm"),
                                     names_pattern = c("(.+)pVal_(.+)"),
                        ) %>% dplyr::rename("p_dTm"="value"))
  df_TPP1<-purrr::map(df_TPP1,function(x)x %>% dplyr::select(Protein_ID,sample_name,diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2,p_dTm) %>% 
                        pivot_longer(c(diff_meltP_Treatment_1_vs_Vehicle_1,diff_meltP_Treatment_2_vs_Vehicle_2),
                                     names_to = c("hi","set"),
                                     names_pattern = "(.+)diff_(.+)"
                        ) %>% dplyr::rename("dTm"="value") %>% dplyr::select(-hi,-set))
  df_TPP1<-purrr::map(df_TPP1,function(x) x[!is.na(x$dTm),])
  if(isTRUE(volcano)){
    
    df_TPP2<-dplyr::bind_rows(df_TPP1) 
    df_TPP2$diffexpressed <- "No"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_TPP2$diffexpressed[df_TPP2$dTm > 2 & df_TPP2$p_dTm < 0.05] <- "Stabilized"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_TPP2$diffexpressed[df_TPP2$dTm < -2 & df_TPP2$p_dTm < 0.05] <- "Destabilized"
    df_TPP2$delabel <- NA
    df_TPP2$delabel[df_TPP2$Protein_ID %in%c("P36507","Q02750")] <- as.character(df_TPP2$Protein_ID[df_TPP2$Protein_ID %in%c("P36507","Q02750")])
    #get overlaps in data
    #df_TPP2$Stroke <- NA
    #df_TPP2$Stroke<-ifelse(df_TPP2$diffexpressed=="Stabilized" | df_TPP2$diffexpressed=="Destabilized",ifelse(df_TPP2$Protein_ID %in% overlaps,TRUE,FALSE),FALSE)
    
    
    
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% distinct(.) %>%  dplyr::group_split(sample_name,Protein_ID)
    df_TPP2<-purrr::map(df_TPP2,function(x) x[1,])
    df_TPP2<-dplyr::bind_rows(df_TPP2) %>% dplyr::group_split(sample_name)
    flevels<-data.frame(treatment=c("No","Destabilized","Stabilized"),
                        colors=c("black","green","purple"))
    
    df_TPP3<-purrr::map(df_TPP2,function(x)
      ggplot(data=x,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+ 
        labs(y=expression(-log["10"]*(P-value)),x=expression(Delta*T["m"]))+
        xlim(-15,15)+
        theme(legend.position="bottom", legend.box = "horizontal")+
        geom_label(aes(16,y=40),label=paste0("S = ",nrow(x[x$diffexpressed=="Stabilized",])),show.legend=FALSE)+
        geom_label(aes(-16,y=40),label=paste0("DS = ",nrow(x[x$diffexpressed=="Destabilized",])),show.legend=FALSE)+
        scale_color_manual("Stabilization",values=flevels$colors[which(flevels$treatment %in% x$diffexpressed)],labels = flevels$treatment[which(flevels$treatment %in% x$diffexpressed)])+
        #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
        geom_label_repel(aes(dTm, -log10(p_dTm),label=delabel),color="black",  nudge_x = 1,
                         force = 3,
                         box.padding = 3,
                         segment.alpha = .5)+
        ggtitle(x$sample_name[1])+ylim(-0.1,55))#+
    #geom_point(data=df_TPP2[df_TPP2$Stroke==1,],
    #pch=21, fill=NA, size=4, colour="black", stroke=1)
    
    
    
    return(df_TPP3)
  }
  
  #for number of curves
  
  df_1<-dplyr::bind_rows(df_TPP1)
  
  
  df_1<-dplyr::bind_rows(df_1) 
  df_1$sample_name<-as.factor(df_1$sample_name)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-df_1 %>% dplyr::select(p_dTm,uniqueID,sample_name,dTm) %>% distinct(.)
  df_1<-df_1[!duplicated(df_1),]
  df_1<-df_1[!is.na(df_1$dTm),]
  df_1$p_dTm<-base::round(df_1$p_dTm,3) 
  df_1<-dplyr::bind_rows(df_1) %>% distinct(.)%>% dplyr::group_split(uniqueID,sample_name)
  df_1<-purrr::map(df_1,function(x) x[1,])
  df_1<-dplyr::bind_rows(df_1)
  
  colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_1)$sample_name)))
  colors1$sample_name<-as.character(colors1$sample_name)[1:length(unique(colors1$sample_name))]
  colors1$sample_name<-as.factor(colors1$sample_name)
  colors1$sample_name<-levels(colors1$sample_name)
  colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
  colors1$x<-c(1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
  colors1$y<-c(1,0.42,-0.38,-0.88,-0.88,-0.38,0.42,1)[1:length(unique(colors1$sample_name))]
  df_TPP1<-df_1
  
  if(!nrow(colors1)==8){
    check1<-df_TPP1 %>% count(sample_name) %>%
      #mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
      ggplot() +
      ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name), stat = "pie") +
      ggforce::theme_no_axes()+
      scale_fill_manual(values = colors1$hex,aesthetics="fill")+
      xlim(-1.1,1.45)+
      ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
      ggtitle("Number of fitted curves")
  }else{
    check1<-df_TPP1 %>% count(sample_name) %>%
      mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
      ggplot() +
      ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name,explode=focus), stat = "pie") +
      ggforce::theme_no_axes()+
      scale_fill_manual(values = colors1$hex,aesthetics="fill")+
      xlim(-1.1,1.45)+
      ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
      ggtitle("Number of fitted curves")
    return(check1)
  }
  # 
  # 'Size'=(
  #   intersection_size(
  #     mode=mode,
  #     mapping=aes(fill=exclusive_intersection),
  #     size=0,
  #     text=list(check_overlap=TRUE)
  #   ) + scale_fill_venn_mix(
  #     data=abc_data,
  #     guide=FALSE,
  #     colors=c('A'='red', 'B'='blue', 'C'='green3')
  # 
  df_1<-dplyr::bind_rows(df_TPP)
  df_1$uniqueID<-df_1$Protein_ID
  
  df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
    pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
  
  
  
  
  
  #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
  # df_ <- df_ %>%
  #   mutate_if(sapply(df_, is.factor), as.numeric)
  IDs<-df_1$uniqueID
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
  df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
  df_1<-cbind(IDs,df_1)
  
  rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                   values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
  df_1<-na.omit(df_1)
  check<-list()
  check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID")],
               #min_degree=6,
               set_sizes=FALSE,
               n_intersections=10,
               min_degree=1,
               encode_sets=FALSE,
               stripes='white',
               # intersections=list(
               #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
               #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
               #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
               #   
               # ),
               # mode = 'inclusive_union',
               #mode=mode,
               # intersections=list(
               #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
               #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
               #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
               #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
               #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
               #   
               # ),
               #sort_intersections_by="cardinality",
               base_annotations=list(
                 '# of fitted curves'=
                   intersection_size(
                     # mode='inclusive_intersection',
                     counts=TRUE,
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   #scale_color_brewer(palette="Dark2")+
                   theme(legend.position = 'none',
                         # axis.text=element_text(size=20, face='bold'),
                         # axis.title=element_text(size=20, face='bold')
                   )+
                   ggplot2::ylab("# of fitted curves")
               ),               # base_annotations=list(
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
               height_ratio=0.8
               
  )+ggtitle(paste0
            ("Top 10 number of fitted curves (TPP) ",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP" | filters=="HQ","filtered","unfiltered")))
  
  #check upset plots for intersections
  check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
  #filter exclusive intersection colors
  colors<-colors1 %>% dplyr::filter(sample_name %in% check_intersections$IDs)
  
  level_data=rev(levels(check$data$intersection))
  #colors$sample_name<-as.character(levels(colors$sample_id))
  queries<-list(
    upset_query(
      intersect=colors$sample_name[1],
      color=colors$hex[1],
      fill=colors$hex[1],
      only_components='# of fitted curves'
    ),
    upset_query(
      intersect=colors$sample_name[2],
      color=colors$hex[2],
      fill=colors$hex[2],
      only_components='# of fitted curves'
    )
    ,
    upset_query(
      intersect=as.character(colors$sample_name[3]),
      color=colors$hex[3],
      fill=colors$hex[3],
      only_components='# of fitted curves'
    )
    ,
    upset_query(
      intersect=colors$sample_name[4],
      color=colors$hex[4],
      fill=colors$hex[4],
      only_components='# of fitted curves'
    )
    ,
    upset_query(
      intersect=colors$sample_name[5],
      color=colors$hex[5],
      fill=colors$hex[5],
      only_components='# of fitted curves'
    )
    ,
    upset_query(
      intersect=colors$sample_name[6],
      color=colors$hex[6],
      fill=colors$hex[6],
      only_components='# of fitted curves'
    ),
    upset_query(
      intersect= as.character(colors$sample_name[7]),
      color=as.character(colors$hex[7]),
      fill=as.character(colors$hex[7]),
      only_components='# of fitted curves'
    )
    ,
    upset_query(
      intersect=as.character(colors$sample_name[8]),
      color=as.character(colors$hex[8]),
      fill=as.character(colors$hex[8]),
      only_components='# of fitted curves'
    ))[1:length(unique(colors$sample_name))]
  if(stringr::str_count(check_intersections$IDs[1],"-")<7){
    queries<-queries
  }else{
    queries<-c(queries,list(
      upset_query(
        intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E','C_F_E'),
        color='black',
        fill='black',
        only_components='# of fitted curves'
      )))
  }
  #add query
  check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
               #min_degree=6,
               set_sizes=FALSE,
               n_intersections=10,
               min_degree=1,
               encode_sets=FALSE,
               stripes='white',
               # intersections=list(
               #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
               #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
               #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
               #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
               #   
               # ),
               # mode = 'inclusive_union',
               #mode=mode,
               # intersections=list(
               #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
               #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
               #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
               #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
               #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
               #   
               # ),
               #sort_intersections_by="cardinality",
               base_annotations=list(
                 '# of fitted curves'=
                   intersection_size(
                     # mode='inclusive_intersection',
                     counts=TRUE
                     # mapping=aes(fill='inclusive_intersection')
                   )+
                   ggplot2::ylab("# of fitted curves")
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
               queries=queries
               # ,
               # upset_query(
               #   intersect=colors$sample_name[4],
               #   color=colors$hex[4],
               #   fill=colors$hex[4],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect=colors$sample_name[5],
               #   color=colors$hex[5],
               #   fill=colors$hex[5],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect=colors$sample_name[6],
               #   color=colors$hex[6],
               #   fill=colors$hex[6],
               #   only_components='# of fitted curves'
               # ),
               # upset_query(
               #   intersect= as.character(colors$sample_name[7]),
               #   color=as.character(colors$hex[7]),
               #   fill=as.character(colors$hex[7]),
               #   only_components='# of fitted curves'
               # )
               # ,
               # upset_query(
               #   intersect=as.character(colors$sample_name[8]),
               #   color=as.character(colors$hex[8]),
               #   fill=as.character(colors$hex[8]),
               #   only_components='# of fitted curves'
               # )
               
               
               
               
  )+ggtitle(paste0("Top 10 number of fitted curves (TPP) (",ifelse(isTRUE(Peptide),"peptide","protein"),"-level ",ifelse(filters=="TPP","additional filters)",ifelse(filters=="HQ","high-quality filtered)",ifelse(isFALSE(filter),"unfiltered)","filtered)")))))
  
  return(list(check,check1))
  
  
}
#CV<-his_sp(res_sp[[3]],df.temps,MD=FALSE)
UpSet_curves<-function(f,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE){
  f<-f %>% purrr::keep(function(x) !class(x)=='try-error')
  if(isTRUE(Peptide)){
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=as.factor(sample_name),treatment=as.factor(treatment)) %>%
      dplyr::select(-C,-I,-temp_ref,-CV_pct,-missing_pct) %>% dplyr::filter(!is.na(sample_name)) %>%
      dplyr::group_split(uniqueID,treatment,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>% group_by(sample_id,treatment) %>% dplyr::summarise(uniqueID=uniqueID,
                                                                                         treatment=treatment,
                                                                                         sample_name=sample_name,
                                                                                         p_dTm=p_dTm,
                                                                                         sample_id=sample_id,
                                                                                         Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                                                                         rss=mean(rss,na.rm=TRUE),
                                                                                         rsq=mean(rsq,na.rm=TRUE),
                                                                                         AUC=mean(AUC,na.rm=TRUE)) %>% 
                    ungroup(.) %>% distinct(.))
  }else if(!isTRUE(Peptide)){
    f<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=sample_name,treatment=as.factor(treatment)) %>% 
      dplyr::select(-rank,-C,-I,-temp_ref,-CV_pct,-missing_pct) %>%
      dplyr::group_split(uniqueID,treatment,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>%
                    group_by(sample_id,treatment) %>%
                    dplyr::summarise(uniqueID=uniqueID,
                                     treatment=treatment,
                                     sample_name=sample_name,
                                     sample_id=sample_id,
                                     p_dTm=p_dTm,
                                     Coverage=Coverage,
                                     MW_kDa=MW_kDa,
                                     Tm=mean(Tm,na.rm=TRUE),#for peptide groups, caculate averages for parameters
                                     rss=mean(rss,na.rm=TRUE),
                                     rsq=mean(rsq,na.rm=TRUE),
                                     AUC=mean(AUC,na.rm=TRUE)) %>%
                    ungroup(.) %>% distinct(.)) 
  }else if (isTRUE(Trilinear)){#if this is a trilinear result
    #f<-f %>% dplyr::group_split(uniqueID,treatment)
    # f<-lapply(f,function(x) dplyr::bind_rows(x))
    # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(sample_name=sample_name,treatment=as.factor(treatment)) %>%
      dplyr::select(-rsq,-data) %>%
      dplyr::group_split(uniqueID,treatment,sample_name)
    
    
    f<-purrr::map(f,function(x) x %>% dplyr::summarise(uniqueID=uniqueID,
                                                       treatment=treatment,
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
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$treatment=="treated"][1]>x$Tm[x$treatment=="vehicle"][1],"Stabilized","Destabilized")),
                                                    dTm=x$Tm[x$treatment=="treated"][1]-x$Tm[x$treatment=="vehicle"][1]))
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    
    df_1<-dplyr::bind_rows(f)%>% dplyr::mutate(sample_name=as.character(sample_name)) %>% dplyr::group_split(uniqueID,sample_name)
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,model_converged,stabilized,rsq_greater_than_0.8) %>% 
      pivot_wider(names_from=sample_name,values_from=c(model_converged)) %>% distinct(.)
    
  }else if(isTRUE(Splines) & !isTRUE(Peptide)){
    
    f<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(Coverage))%>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name,sample_id)
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm>0 & p_dTm<0.05,"Stabilized","Destabilized")))
    
    f<-purrr::map(f,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    df_TPP<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(sample_name)) %>% 
      dplyr::select(uniqueID,sample_name,Tm,rss,AUC,p_dTm,stabilized,sample_id) %>% 
      distinct(.)
    
    df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,sample_id,Tm) %>% distinct(.)
    df_1<-dplyr::bind_rows(df_1)
    df_2<-df_1
    colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors1$sample_name<-as.factor(colors1$sample_name)
    colors1$sample_name<-levels(colors1$sample_name)
    colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
    colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
    colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
    df_1<-df_1[!duplicated(df_1),]
    
    if(length(unique(df_1$sample_name))==1){
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = 0) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }else{
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = ifelse(sample_name == "C_F_E", 0.2, 0)) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }
    
    # pdf("Number_of_curves_upset_splines_Protein.pdf",encoding="CP1253.enc",compress=FALSE,width=6.12,height=4.02)
    # check1
    # dev.off()
    
    df_1<-dplyr::bind_rows(df_2)
    
    df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name) %>% 
      pivot_wider(names_from=sample_name,values_from=sample_name) %>% distinct(.)
    
    
    
    
    #df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
    
    IDs<-as.factor(df_1$uniqueID)
    if(!ncol(df_1)==2){
      df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., !is.na(.), "TRUE"))
      df_1 <- mutate_all(df_1[,2:length(df_1)], ~replace(., is.na(.), "FALSE"))
    }
    df_1<-cbind(IDs,df_1)
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    df_1<-na.omit(df_1)
    
    #keep the first row out of redundant peptide group data
    
    #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
    # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
    #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
    #                           RSQ=mean(Rsq,na.rm=TRUE)))
    # f<-dplyr::bind_rows(f) 
    # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
    #                        rss_n=ifelse(treatment=="null",RSS,NA))
    # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
    
    #select columns of interest
    
    
    #,rsq_greater_than_0.8,stabilized
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    
    if(isTRUE(filter)){
      check<-list()
      
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=FALSE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   height_ratio=0.8
                   
      )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(filter=="TRUE","filtered)","unfiltered)")))
      
      level_data=rev(levels(check$data$intersection))
      #colors$sample_name<-as.character(levels(colors$sample_id))
      colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors$sample_name<-as.factor(colors$sample_name)
      colors$sample_name<-levels(colors$sample_name)
      colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')
      colors$sample_name<-levels(colors$sample_name)
      colors<-dplyr::bind_rows(colors) %>% dplyr::filter(sample_name %in% level_data)
      
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      #check upset plots for intersections
      check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
      
      if(str_count(check_intersections$IDs[1],"-")<7){
        queries<-queries
      }else{
        c(queries,list(
          upset_query(
            intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
            color='black',
            fill='black',
            only_components='# of fitted curves'
          )))
      }
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   queries<-queries
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines","trilinear")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      
      df_<-df_[apply(df_!=0, 1, all),]
      return(list(check,check1,df_))
    }else{#if the data isnt filtered
      
      f<-dplyr::bind_rows(f) %>% 
        distinct(.) %>% dplyr::group_split(uniqueID,sample_name,sample_id)
      
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm>0 & p_dTm<0.01,"Stabilized","Destabilized")))
      
      f<-purrr::map(f,function(x) x[1,])
      
      f<-dplyr::bind_rows(f) %>% dplyr::select(-sample_id) %>% 
        f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      df_TPP<-dplyr::bind_rows(f) %>% dplyr::filter(!is.na(sample_name)) %>% 
        dplyr::select(uniqueID,sample_name,Tm,rss,AUC,p_dTm,stabilized) %>% 
        distinct(.)
      
      df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,Tm) %>% distinct(.)
      df_1<-dplyr::bind_rows(df_1)
      df_2<-df_1
      colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors1$sample_name<-as.factor(colors1$sample_name)
      colors1$sample_name<-levels(colors1$sample_name)
      colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
      colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
      colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      
      df_1<-df_1[!duplicated(df_1),]
      if(length(unique(df_1$sample_name))==1){
        check1<-df_1 %>% count(sample_name) %>%
          mutate(focus = 0) %>%
          ggplot() +
          ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
          ggforce::theme_no_axes()+
          scale_fill_manual(values = colors1$hex,aesthetics="fill")+
          xlim(-1.1,1.45)+
          ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
          ggtitle("Number of fitted curves")+
          theme(legend.position="bottom", legend.box = "horizontal")
        return(check1)
      }else{
        check1<-df_1 %>% count(sample_name) %>%
          mutate(focus = ifelse(sample_name == "C_F_E", 0.2, 0)) %>%
          ggplot() +
          ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
          ggforce::theme_no_axes()+
          scale_fill_manual(values = colors1$hex,aesthetics="fill")+
          xlim(-1.1,1.45)+
          ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
          ggtitle("Number of fitted curves")+
          theme(legend.position="bottom", legend.box = "horizontal")
      }
      
      df_1<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,Tm) %>% 
        pivot_wider(names_from=sample_name,values_from=c(Tm)) %>% distinct(.)
      
      df_1<-df_1 %>% dplyr::mutate(uniqueID=as.character(uniqueID))
      df_1 <- df_1 %>%
        mutate_if(sapply(df_1, is.factor), as.numeric)
      df_1$uniqueID<-as.factor(df_1$uniqueID)
      df_1<-mutate_all(df_1, ~replace(., !is.na(.) & is.numeric(.), "TRUE"))
      df_1 <- mutate_all(df_1, ~replace(., is.na(.), "FALSE"))
      df_1 <- mutate_all(df_1, ~replace(., is.character(.), as.logical(.)))
      
      #keep the first row out of redundant peptide group data
      
      #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
      # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
      #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
      #                           RSQ=mean(Rsq,na.rm=TRUE)))
      # f<-dplyr::bind_rows(f) 
      # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
      #                        rss_n=ifelse(treatment=="null",RSS,NA))
      # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
      
      #select columns of interest
      
      
      #,rsq_greater_than_0.8,stabilized
      rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                       values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
      
      check<-list()
      
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% "uniqueID"],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=FALSE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   height_ratio=0.8
                   
      )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(filter=="TRUE","filtered)","unfiltered)")))
      
      level_data=rev(levels(check$data$intersection))
      #colors$sample_name<-as.character(levels(colors$sample_id))
      colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
      colors$sample_name<-as.factor(colors$sample_name)
      colors$sample_name<-levels(colors$sample_name)
      colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors$sample_name))]
      
      colors<-dplyr::bind_rows(colors) %>% dplyr::filter(sample_name %in% level_data)
      
      queries<-list(
        upset_query(
          intersect=colors$sample_name[1],
          color=colors$hex[1],
          fill=colors$hex[1],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect=colors$sample_name[2],
          color=colors$hex[2],
          fill=colors$hex[2],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[3]),
          color=colors$hex[3],
          fill=colors$hex[3],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[4],
          color=colors$hex[4],
          fill=colors$hex[4],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[5],
          color=colors$hex[5],
          fill=colors$hex[5],
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=colors$sample_name[6],
          color=colors$hex[6],
          fill=colors$hex[6],
          only_components='# of fitted curves'
        ),
        upset_query(
          intersect= as.character(colors$sample_name[7]),
          color=as.character(colors$hex[7]),
          fill=as.character(colors$hex[7]),
          only_components='# of fitted curves'
        )
        ,
        upset_query(
          intersect=as.character(colors$sample_name[8]),
          color=as.character(colors$hex[8]),
          fill=as.character(colors$hex[8]),
          only_components='# of fitted curves'
        ))[1:length(colors$sample_name)]
      #check upset plots for intersections
      check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
      
      if(str_count(check_intersections$IDs[1],"-")<7){
        queries<-queries
      }else{
        c(queries,list(
          upset_query(
            intersect=c('C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
            color='black',
            fill='black',
            only_components='# of fitted curves'
          )))
      }
      check<-upset(df_1,colnames(df_1)[!colnames(df_1) %in% c("sample_name","stabilized","uniqueID","IDs")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   queries<-queries
                   
                   
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines","trilinear")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      
      df_<-df_[apply(df_!=0, 1, all),]
    }
    return(list(check,check1,df_))
  }else{#if this is a peptide spline file
    
    f1<-dplyr::bind_rows(f) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_id,sample_name)
    
    
    f<-purrr::map(f1,function(x) x[1,])
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    df_TPP<-dplyr::bind_rows(f) %>%dplyr::filter(!is.na(sample_name)) %>% 
      dplyr::select(uniqueID,sample_name,rss,AUC,p_dTm,Tm) %>% 
      distinct(.)
    
    df_1<-df_TPP %>% dplyr::select(p_dTm,uniqueID,sample_name,Tm) %>% distinct(.)
    
    df_1<-dplyr::bind_rows(df_1)
    df_2<-df_1
    colors1<-data.frame(sample_name=as.character(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors1$sample_name<-as.factor(colors1$sample_name)
    colors1$sample_name<-levels(colors1$sample_name)
    colors1$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')[1:length(unique(colors1$sample_name))]
    colors1$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors1$sample_name))]
    colors1$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors1$sample_name))]
    
    
    
    df_1<-df_2[!duplicated(df_2),]
    if(length(unique(df_1$sample_name))==1){
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = 0) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
      return(check1)
    }else{
      check1<-df_1 %>% count(sample_name) %>%
        mutate(focus = ifelse(sample_name == "C_F_Φ", 0.2, 0)) %>%
        ggplot() +
        ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n, fill = sample_name, explode = focus), stat = "pie") +
        ggforce::theme_no_axes()+
        scale_fill_manual(values = colors1$hex,aesthetics="fill")+
        xlim(-1.1,1.45)+
        ggplot2::geom_label(mapping=aes(x=colors1$x,y=colors1$y,label=n,color=sample_name),inherit.aes=TRUE,vjust="top",show.legend = FALSE)+
        ggtitle("Number of fitted curves")+
        theme(legend.position="bottom", legend.box = "horizontal")
    }
    
    
    # pdf("Number of fitted_curves_Peptide_filt.pdf",encoding="CP1253.enc",compress=TRUE,width=5.31,height=4.02)
    # check1
    # dev.off()
    
    df_<-dplyr::bind_rows(df_1) %>% dplyr::select(uniqueID,sample_name,Tm) %>% 
      pivot_wider(names_from=sample_name,values_from=c(Tm)) %>% distinct(.)
    
    df_<-df_ %>% dplyr::mutate(uniqueID=as.character(uniqueID))
    df_ <- df_ %>%
      mutate_if(sapply(df_, is.factor), as.numeric)
    df_$uniqueID<-as.factor(df_$uniqueID)
    df_<-mutate_all(df_, ~replace(., !is.na(.) & is.numeric(.), "TRUE"))
    df_ <- mutate_all(df_, ~replace(., is.na(.), "FALSE"))
    df_ <- mutate_all(df_, ~replace(., is.character(.), as.logical(.)))
    
    #keep the first row out of redundant peptide group data
    
    #f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=ifelse(x$Tm[x$treatment=="treated"]>x$Tm[x$treatment=="vehicle"],1,0)))
    # f<-lapply(f,function(x) x %>% dplyr::group_by(uniqueID,treatment) %>% 
    #             dplyr::mutate(RSS=sum(rss,na.rm=TRUE),
    #                           RSQ=mean(Rsq,na.rm=TRUE)))
    # f<-dplyr::bind_rows(f) 
    # f<-f %>% dplyr::mutate(rss_a=ifelse(treatment=="vehicle"|treatment=="treated",RSS,NA),
    #                        rss_n=ifelse(treatment=="null",RSS,NA))
    # f<-f %>% dplyr::ungroup(.) %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(rss_a=sum(unique(rss_a),na.rm=TRUE))
    
    #select columns of interest
    
    
    #,rsq_greater_than_0.8,stabilized
    rating_scale = scale_fill_manual(name="Stabilization (Tm-based)",
                                     values=c("Stabilized" ='#fee6ce', "Destabilized" ='#fdae6b', "NA"  = '#e6550d'))
    
    check<-list()
    
    check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                 #min_degree=6,
                 set_sizes=FALSE,
                 n_intersections=10,
                 min_degree=1,
                 encode_sets=FALSE,
                 stripes='white',
                 # intersections=list(
                 #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                 #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                 #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                 #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                 #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ")
                 #   
                 # ),
                 # mode = 'inclusive_union',
                 #mode=mode,
                 # intersections=list(
                 #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                 #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                 #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                 #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                 #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                 #   
                 # ),
                 #sort_intersections_by="cardinality",
                 base_annotations=list(
                   '# of fitted curves'=
                     intersection_size(
                       # mode='inclusive_intersection',
                       counts=TRUE
                       # mapping=aes(fill='inclusive_intersection')
                     )+
                     ggplot2::ylab("# of fitted curves")
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
                 height_ratio=0.8
                 
    )+ggtitle(paste0("Top 10 number of fitted curves splines ",ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",ifelse(isTRUE(filter),"filtered)","unfiltered)")))
    
    
    level_data=rev(levels(check$data$intersection))
    check_intersections<-data.frame(IDs=check[[1]]$data$intersection,inclusive=check[[1]]$data$inclusive_intersection_size) %>% distinct(.) %>% dplyr::arrange(inclusive) %>% head(10)
    
    #colors$sample_name<-as.character(levels(colors$sample_id))
    colors<-data.frame(sample_name=as.factor(unique(dplyr::bind_rows(df_TPP)$sample_name)))
    colors$sample_name<-levels(colors$sample_name)
    colors$hex<-c('#d07884','#ffb12c','#7adf68','#40bc39','#12a7c8','#404898','#ac5180','#ec5481')
    colors$x<-c(1.1,1.3,1.3,0.7,-0.8,-1,-1.1,-0.75)[1:length(unique(colors$sample_name))]
    colors$y<-c(1,0.42,-0.48,-0.88,-0.78,-0.48,0.42,1)[1:length(unique(colors$sample_name))]
    
    colors<-colors %>% dplyr::filter(sample_name %in% check_intersections$IDs)
    queries=list(
      upset_query(
        intersect=colors$sample_name[1],
        color=colors$hex[1],
        fill=colors$hex[1],
        only_components='# of fitted curves'
      ),
      upset_query(
        intersect=colors$sample_name[2],
        color=colors$hex[2],
        fill=colors$hex[2],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=as.character(colors$sample_name[3]),
        color=colors$hex[3],
        fill=colors$hex[3],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[4],
        color=colors$hex[4],
        fill=colors$hex[4],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[5],
        color=colors$hex[5],
        fill=colors$hex[5],
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=colors$sample_name[6],
        color=colors$hex[6],
        fill=colors$hex[6],
        only_components='# of fitted curves'
      ),
      upset_query(
        intersect= as.character(colors$sample_name[7]),
        color=as.character(colors$hex[7]),
        fill=as.character(colors$hex[7]),
        only_components='# of fitted curves'
      )
      ,
      upset_query(
        intersect=as.character(colors$sample_name[8]),
        color=as.character(colors$hex[8]),
        fill=as.character(colors$hex[8]),
        only_components='# of fitted curves'
      ))[1:length(colors$sample_name)]
    #check upset plots for intersections
    
    if(str_count(check_intersections$IDs[1],"-")<7){
      queries<-queries
    }else{
      queries=c(queries,list(
        upset_query(
          intersect=c('C_F_E','C_F_Φ', 'nC_F_Φ','nC_F_E',"C_nF_Φ",'nC_nF_Φ','C_nF_E','nC_nF_E'),
          color='black',
          fill='black',
          only_components='# of fitted curves'
        )))
    }
    if(!isTRUE(filter)){#if the data is not filtered
      check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   queries=queries
                   
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines ","trilinear ")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      df_<-df_[apply(df_!=0, 1, all),]
      
      return(list(check,check1,df_))
      
    }else{#if the peptide data is unfiltered
      check<-upset(df_,colnames(df_)[!colnames(df_) %in% c("sample_name","stabilized","uniqueID")],
                   #min_degree=6,
                   set_sizes=FALSE,
                   n_intersections=10,
                   min_degree=1,
                   encode_sets=TRUE,
                   stripes='white',
                   # intersections=list(
                   #   c("C_F_E","C_F_Φ","C_nF_E","nC_F_E","nC_nF_E","nC_nF_Φ","C_nF_Φ","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","C_nF_E","C_nF_Φ"),
                   #   c("nC_nF_E","nC_nF_Φ","nC_nF_E","nC_F_Φ"),
                   #   c("C_F_E","C_F_Φ","nC_F_E","nC_F_Φ"),
                   #   c("C_nF_E","C_nF_Φ","nC_nF_E","nC_nF_Φ"),
                   #   
                   # ),
                   # mode = 'inclusive_union',
                   #mode=mode,
                   # intersections=list(
                   #   c('C_F_E', paste0('C_F_','\u03A6'),'C_nF_E',paste0('C_nF_','\u03A6'),'nC_F_E',paste0('nC_F_','\u03A6'),'nC_nF_E',paste0('nC_F_','\u03A6')),
                   #   c(paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6'),paste0('nC_F_','\u03A6'),paste0('nC_F_','\u03A6')),
                   #   c('C_F_E','C_nF_E','nC_F_E','nC_nF_E'),
                   #   c('C_F_E','C_nF_E',paste0('C_F_','\u03A6'),paste0('C_nF_','\u03A6')),
                   #   c('C_F_E','nC_F_E',paste0('C_F_','\u03A6'),paste0('nC_F_','\u03A6'))
                   #   
                   # ),
                   #sort_intersections_by="cardinality",
                   base_annotations=list(
                     '# of fitted curves'=
                       intersection_size(
                         # mode='inclusive_intersection',
                         counts=TRUE
                         # mapping=aes(fill='inclusive_intersection')
                       )+
                       ggplot2::ylab("# of fitted curves")
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
                   queries=queries
      )+ggtitle(paste0("Top 10 number of fitted curves ", ifelse(isTRUE(Splines),"splines ","trilinear ")
                       ,ifelse(isTRUE(Peptide),"(peptide","(protein"),"-level ",
                       ifelse(isTRUE(filter),"filtered)","unfiltered)")))
      df_<-df_[apply(df_!=0, 1, all),]
      
      return(list(check,check1,df_))
    }
  } 
  
}

volcano_data<-function(f,Trilinear=FALSE,Splines=FALSE,Sigmoidal=TRUE,Peptide=FALSE,benchmark=TRUE,Fhist=TRUE,labels=FALSE,type=NA){
  f<-f %>% purrr::keep(function(x) !class(x)=='try-error')
  f<-data.frame(dplyr::bind_rows(f))
  if(any(names(f)=="sample_id")){
    if(any(names(f)=="sample_id")){
      f<-f %>% dplyr::select(-sample_id)
    }
    f<-f %>% dplyr::rename("sample_id"="sample_id")
  }
  if(any(names(f)=="I3")){
    f<-f %>% dplyr::rename("I"="I3")
  }
  if(isTRUE(benchmark)){
    if(isTRUE(Peptide)){
      
      f<-f %>% 
        dplyr::select(uniqueID,sample_id,sample_name,dTm,p_dTm,AUC,rsq,fStatistic,pValue,pAdj,treatment,Fvals,p_dTm,dRSS,d1,d2,rss1,RSS)
      f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      f<-dplyr::bind_rows(f)
      f$diffexpressed <- "Not Shifted"
      f$diffexpressed[f$dTm > 0 & f$p_dTm<0.05] <- "Stabilized Tm"
      f$diffexpressed[f$dTm < 0& f$p_dTm<0.05] <- "Destabilized Tm"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm > 0 & f$p_dTm<0.05] <- "Stabilized Shift"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm < 0 & f$p_dTm<0.05] <- "Destabilized Shift"
      f$diffexpressed<-as.factor(f$diffexpressed)
      f$targets<-NA
      f$delabel <- NA
      others<-as.character(f$uniqueID[which(f$pAdj < 0.05)])
      f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
      f$delabel[f$uniqueID %in%c("P36507","Q02750",others)] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750",others)])
      flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                          labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
      flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
      
      
      f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f)
      if(isTRUE(benchmark)){
        if(!isTRUE(Fhist)){
          check<-ggplot(data=f,mapping=aes(x=dRSS,y=-log10(pAdj),color=diffexpressed))+geom_point()+
            geom_hline(yintercept=-log10(0.05), col="red")+
            scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
            labs(x=expression(RSS["0"]-RSS["1"]),y=expression(-log["10"]*p["adj"]-value))+
            #geom_label(aes(dTm, pAdj, label = delabel), data = f,color="black")+
            # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
            #                 nudge_x = 1,
            #                 force = 2,
            #                 box.padding = 2,
            #                 segment.alpha = .5)+
            ggtitle(f$sample_name[1])+
            theme(legend.position="bottom", legend.box = "horizontal")+
            xlim(0,max(f$Fvals,na.rm=TRUE))#+
          #ylim(0,4) 
          if(isTRUE(labels)){
            if(type=="targets" | is.na(type)){
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)+
                xlim(-max(sqrt(f$dRSS),na.rm=TRUE),max(sqrt(f$dRSS),na.rm=TRUE))
            }else{
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }
          }
          
        }else{
          
          check<-ggplot(f)+
            geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) +
            theme_bw() +
            coord_cartesian(xlim=c(0,10))+
            ggplot2::xlab("F-values")+
            xlim(0,max(f$Fvals))
        }
        
      }
      
    }else if(!isTRUE(Peptide)){
      
      f<-dplyr::bind_rows(f) %>% 
        dplyr::select(uniqueID,sample_id,sample_name,dTm,AUC,rsq,fStatistic,pValue,pAdj,treatment,Fvals,p_dTm,dRSS,d1,d2,rss1,RSS)
      f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
      
      f<-dplyr::bind_rows(f)
      f$diffexpressed <- "Not Shifted"
      f$diffexpressed[f$dTm > 0] <- "Stabilized Tm"
      f$diffexpressed[f$dTm < 0] <- "Destabilized Tm"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm > 0] <- "Stabilized Shift"
      f$diffexpressed[f$Fvals < f$fStatistic & f$pAdj < 0.05 & f$dTm < 0] <- "Destabilized Shift"
      f$diffexpressed<-as.factor(f$diffexpressed)
      f$targets<-NA
      f$delabel <- NA
      others<-as.character(f$uniqueID[which(f$pAdj < 0.05)])
      f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
      f$delabel[f$uniqueID %in%c("P36507","Q02750",others)] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750",others)])
      flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                          labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
      flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
      
      f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f)
      if(isTRUE(benchmark)){
        if(!isTRUE(Fhist)){
          check<-ggplot(data=f,mapping=aes(x=dRSS,y=-log10(pAdj),color=diffexpressed))+geom_point()+
            geom_hline(yintercept=-log10(0.05), col="red")+
            scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
            labs(x=expression(RSS["0"]-RSS["1"]),y=expression(-log["10"]*p["adj"]-value))+
            #geom_label(aes(dTm, pAdj, label = delabel), data = f,color="black")+
            # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
            #                 nudge_x = 1,
            #                 force = 2,
            #                 box.padding = 2,
            #                 segment.alpha = .5)+
            ggtitle(f$sample_name[1])+
            theme(legend.position="bottom", legend.box = "horizontal")+
            coord_cartesian(xlim=c(0,max(f$dRSS,na.rm=TRUE)),ylim=c(0,max(-log10(f$pAdj),na.rm=TRUE)))
          if(isTRUE(labels)){
            if(type=="targets" | is.na(type)){
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }else{
              check<-check+
                geom_label_repel(f,mapping=aes(dRSS, -log10(pAdj),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                                 nudge_x = 1,
                                 force = 3,
                                 box.padding = 3,
                                 segment.alpha = .5)
            }
          }
          
        }else{
          
          check<-ggplot(f)+
            geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) +
            theme_bw() +
            coord_cartesian(xlim=c(0,10))+
            ggplot2::xlab("F-values")+
            xlim(0,max(f$Fvals))
          return(check)
        }
        
      }
      return(check)
      
    }
  }else if(isFALSE(benchmark)){
    if(isTRUE(Peptide)){
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-f %>% 
        dplyr::group_split(uniqueID,treatment,sample_id)
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(Tm=try(stats::approx(.$I,.$C, xout=(min(.$I,na.rm=TRUE)+(0.5*(max(.$I, na.rm=TRUE)-min(.$I, na.rm=TRUE)))),n=100)$y)))
      f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
      f<-purrr::map(f,function(x) x[1,])
      f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID)
      
      
    }else if(!isTRUE(Peptide)){
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-data.frame(f)
      f$sample_name<-as.factor(f$sample_name)
      f$sample_id<-as.factor(f$sample_id)
      f<-f %>%
        dplyr::group_split(uniqueID,treatment,sample_id)
      f<-f %>% purrr::keep(function(x) nrow(x)>0)
      
      f<-dplyr::bind_rows(f)
      
      
    }else if (isTRUE(Trilinear)){#if this is a trilinear result
      #f<-f %>% dplyr::group_split(uniqueID,treatment)
      # f<-lapply(f,function(x) dplyr::bind_rows(x))
      # f<-lapply(f,function(x) x %>% dplyr::mutate(sample_name=x$data[[1]]$sample_name[1]))
      f<-f %>% dplyr::filter(!is.na(sample_name))
      hi<-f %>% dplyr::group_by(sample_name) %>% dplyr::summarise(uniqueID=uniqueID,sample_id=sample_id,sample_name=sample_name,ptest=calcP(uniqueID,Tm,dTm,30))
      f<-f %>% dplyr::right_join(hi,by=dplyr::intersect(names(hi),names(f)))
      f<-dplyr::bind_rows(f) %>%
        dplyr::group_by(uniqueID,treatment,sample_id) %>% # group individual curve data to calculate individual Tm
        dplyr::mutate(sample_name=as.factor(sample_name),
                      treatment=as.factor(treatment),
                      Tm=with(.,stats::approx(I,C, xout=min(I,na.rm=TRUE)+(0.5*(max(I, na.rm=TRUE)-min(I, na.rm=TRUE))))$y)) %>% 
        dplyr::select(-rsq,-CI,-data) %>% dplyr::ungroup(.)
      f<-f %>% 
        dplyr::group_split(uniqueID)
      
      f<-purrr::map(f,function(x) x %>% dplyr::mutate(
        # v_Tm=mean(x$Tm[x$treatment == "vehicle"],na.rm=TRUE),
        #             t_Tm=mean(x$Tm[x$treatment == "treated"],na.rm=TRUE),
        dTm=mean(x$Tm[x$treatment == "treated"],na.rm=TRUE)-mean(x$Tm[x$treatment == "vehicle"],na.rm=TRUE),
        #FC=log2(v_Tm/t_Tm),
        variance_equal_vt = var.test(x$Tm ~ x$treatment)$p.value,
        p_dTm= try(ifelse(all(variance_equal_vt < 0.05),
                          t.test(Tm ~ treatment, data = x,
                                 var.equal = ifelse(all(variance_equal_vt<0.05),FALSE,TRUE))$p.value[1],NA)),
        p_dTm = p.adjust(p_dTm,"BH")))
      
      f<-f %>% purrr::keep(function(x) !class(x$p_dTm)=='try-error')
      f<-dplyr::bind_rows(f) 
      # f<-purrr::map(f,function(x) x %>% dplyr::mutate(p_dTm=calcP(uniqueID,Tm,dTm,300)))
    }
  }else{
    f$sig<-sign(f$dTm)
    fpos<-f %>% dplyr::filter(sig>0)
    fneg<-f %>% dplyr::filter(sig<0)
    f1<-fpos %>%
      dplyr::mutate(pV = (1-stats::pf(sqrt(dRSS),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f1<-f1 %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    
    f2<-fneg %>%
      dplyr::mutate(pV = (1-stats::pf(sqrt(dRSS),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f2<-f2 %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    f<-dplyr::bind_rows(f1,f2)
    f<-f %>%
      dplyr::mutate(pV = (1-stats::pf(log2(Fvals+1),df1=as.numeric(f$d1[1]),df2=as.numeric(f$d2[1]))))
    f<-f %>% dplyr::mutate(pVAdj = p.adjust(.$pV,method="BH"))
    
    # ggplot(f)+
    #   geom_density(aes(x=log2(Fvals+1),fill = "steelblue",alpha = 0.5))+ 
    #   geom_line(aes(x=log2(Fvals+1),y= df(log2(Fvals+1),df1=4,df2=8)),color="darkred",size = 1.5)
    
    # 
    # M<-median(f$dRSS,na.rm=TRUE)
    # V<-mad(f$dRSS,na.rm=TRUE)
    # #alternative scaling factor sig0-sq
    # altScale<-0.5*V/M
    # #filter out negative delta rss
    # f<-f %>% dplyr::filter(dRSS>0)
    # #effective degrees of freedom
    # ed1<-MASS::fitdistr(x=f$dRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # ed2<-MASS::fitdistr(x=f$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # #scale data
    # f <-f %>% 
    #   dplyr::mutate(rssDiff = .$dRSS/altScale,
    #                 rss1 =.$rss1/altScale,
    #                 d1=ed1,
    #                 d2=ed2)
    # #
    # #new F-test
    # f<-f %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
    # 
    # d1<-f$d1
    # d2<-f$d2
    # 
    # #RSS1 numerator
    # f$fNum<-f$dRSS/d1
    # #Rss denominator
    # f$fDen<-f$RSS/d2
    # #     
    # f$fStatistic=f$fNum/f$fDen
    
    f$diffexpressed <- "Not Shifted"
    f$diffexpressed[f$dTm > 2 & f$p_dTm<0.01] <- "Stabilized Tm"
    f$diffexpressed[f$dTm < -2 & f$p_dTm<0.01] <- "Destabilized Tm"
    f$diffexpressed[f$Fvals > f$fStatistic & f$pV<0.05 & f$dTm > 2 & f$p_dTm<0.01] <- "Stabilized Shift"
    f$diffexpressed[f$Fvals > f$fStatistic &f$pV<0.05 & f$dTm < -2 & f$p_dTm<0.01] <- "Destabilized Shift"
    f$diffexpressed<-as.factor(f$diffexpressed)
    f$targets<-NA
    f$delabel <- NA
    others<-as.character(f$uniqueID[which(f$pAdj < 0.01)])
    f$targets[f$uniqueID %in%c("P36507","Q02750")] <- as.character(f$uniqueID[f$uniqueID %in%c("P36507","Q02750")])
    f$delabel[f$uniqueID %in%others] <- as.character(f$uniqueID[f$uniqueID %in%others])
    flevels<-data.frame(colors=c("#762a83","#af8dc3","#b35806","#7fbf7b","#1b7837"),
                        labels=c("Destabilized Shift","Destabilized Tm","Not Shifted","Stabilized Tm","Stabilized Shift"))
    flevels<-flevels %>% dplyr::filter(labels %in% unique(f$diffexpressed))
    f<-data.frame(f) %>% dplyr::group_split(uniqueID,sample_id)
    f<-purrr::map(f,function(x) x[1,])
    f<-dplyr::bind_rows(f)
    if(!isTRUE(Fhist)){
      check<-ggplot(data=f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),color=diffexpressed))+geom_point()+
        #geom_hline(yintercept=-log10(0.01), col="red")+ 
        labs(y=expression(log["2"]*F["vals"]+1),x=expression(sign("k")*sqrt(RSS["0"]-RSS["1"])))+
        scale_color_manual("Stabilization",values=as.character(flevels$colors),labels = flevels$labels)+
        ggtitle(f$sample_name[1])+
        theme(legend.position="bottom", legend.box = "horizontal")+
        xlim(-max(sqrt(f$dRSS),na.rm=TRUE),max(sqrt(f$dRSS),na.rm=TRUE))
      if(isTRUE(labels)){
        if(type=="targets" | is.na(type)){
          check<-check+
            geom_label_repel(f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                             nudge_x = 1,
                             force = 3,
                             box.padding = 3,
                             segment.alpha = .5)
        }else{
          check<-check+
            geom_label_repel(f,mapping=aes(y=log2(Fvals+1),x=sig*sqrt(dRSS),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                             nudge_x = 1,
                             force = 3,
                             box.padding = 3,
                             segment.alpha = .5)
        }
      }
      
      
    }else{
      check<-ggplot(f)+
        geom_density(aes(x=Fvals),fill = "black",alpha = 0.5) + 
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
    }
    return(check)
  }
  
  if(isTRUE(Trilinear)){
    
    f<-dplyr::bind_rows(f) %>% dplyr::group_split(uniqueID,sample_name)
    f<-f %>% purrr::keep(function(x) any(class(x$M1[[1]])=="lm"))
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(x$Tm[x$treatment=="treated"][1]>x$Tm[x$treatment=="vehicle"][1],"Stabilized","Destabilized"))))
    
    
    
    f<-dplyr::bind_rows(f) %>% dplyr::mutate(stabilized=as.factor(stabilized))
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    f<-purrr::map(f,function(x) x[1,])
    
    f<-f %>% dplyr::mutate(model_converged=as.factor(ifelse(class(M1[[1]])=="lm",1,0)),
                           rsq_greater_than_0.8=as.factor(ifelse(rsq>0.8,1,0)))
    df_<-dplyr::bind_rows(f) %>% 
      select(uniqueID,FC,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,rsq_greater_than_0.8,model_converged,stabilized) %>% 
      distinct(.)
    df_$diffexpressed <- "Not shifted"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 2 & df_$p_dTm < 0.05] <- "Stabilized Shift"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -2 & df_$p_dTm < 0.05] <- "Destabilized Shift"
    df_$delabel <- NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750","P60033")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750","P60033")])
    
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
      labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
      #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_)+
      geom_text_repel(df_,mapping=aes(Fvals, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")
    
  }else if(isTRUE(Splines)){
    
    f<-dplyr::bind_rows(f) %>% 
      distinct(.) %>% dplyr::group_split(uniqueID,sample_name)
    
    f<-purrr::map(f,function(x) x %>% dplyr::mutate(stabilized=as.factor(ifelse(dTm>0 & ptest<0.01,"Stabilized Shift",ifelse(dTm<0 & ptest<0.01,"Destabilized Shift","No")))))
    
    f<-dplyr::bind_rows(f) 
    f$sample_name<-str_replace(f$sample_name,"S","\u03A6")
    
    
    df_<-f %>% 
      dplyr::select(uniqueID,p_dTm,sample_name,dTm,Tm,rss,AUC,rsq,stabilized,ptest,pAdj) %>% 
      distinct(.)
    df_$diffexpressed <- "Not Shifted"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    df_$diffexpressed[df_$dTm > 2 & df_$ptest < 0.01] <- "Stabilized Shift"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    df_$diffexpressed[df_$dTm < -2 & df_$ptest < 0.01] <- "Destabilized Shift"
    
    df_$delabel <- NA
    Names<-as.character(df_$uniqueID[which(df_$pAdj<0.05 & !df_$diffexpressed=="Not Shifted")])
    
    df_$targets<-NA
    df_$delabel[df_$uniqueID %in%c("P36507","Q02750",Names)] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750",Names)])
    df_$targets[df_$uniqueID %in%c("P36507","Q02750")] <- as.character(df_$uniqueID[df_$uniqueID %in%c("P36507","Q02750")])
    df_$diffexpressed<-as.factor(df_$diffexpressed)
    flevels<-data.frame(colors=c("blue","black","green"))
    df_<-df_ %>% dplyr::select(-Tm,-rss,-AUC,-rsq,-stabilized)
    #df_<-df_[!duplicated(df_),]
    check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(ptest),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
      geom_hline(yintercept=-log10(0.01), col="red")+
      scale_color_manual("Stabilization",values=flevels$colors[1:length(levels(df_$diffexpressed))],labels = levels(df_$diffexpressed))+
      labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
      
      # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
      #                 nudge_x = 1,
      #                 force = 2,
      #                 box.padding = 2,
      #                 segment.alpha = .5)+
      ggtitle(df_$sample_name[1])+
      theme(legend.position="bottom", legend.box = "horizontal")+
      geom_label(aes(x=10,y=max(-log10(ptest),na.rm=TRUE)),label=paste0("S = ",nrow(df_[df_$diffexpressed=="Stabilized Shift",])),show.legend=FALSE)+
      geom_label(aes(x=-10,y=max(-log10(ptest),na.rm=TRUE)),label=paste0("DS = ",nrow(df_[df_$diffexpressed=="Destabilized Shift",])),show.legend=FALSE)+
      xlim(-15,15)
    
    if(isTRUE(labels)){
      if(type=="targets" | is.na(type)){
        check<-check+
          geom_label_repel(df_,mapping=aes(dTm, -log10(ptest),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           nudge_x = 1,
                           force = 3,
                           box.padding = 3,
                           segment.alpha = .5)
      }else{
        check<-check+
          geom_label_repel(df_,mapping=aes(dTm, -log10(ptest),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           nudge_x = 1,
                           force = 3,
                           box.padding = 3,
                           segment.alpha = .5)
      }
    }
    
    
  }
  return(check)
}
#     df_$diffexpressed<-as.factor(df_$diffexpressed)
#     flevels<-data.frame(colors=c("blue","black","red"))
#     check<-ggplot(data=df_,mapping=aes(x=dTm,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
#       geom_hline(yintercept=-log10(0.05), col="red")+ scale_color_manual("Stabilization",values=flevels$colors[1:length(unique(df_$diffexpressed))],labels = unique(df_$diffexpressed))+
#       labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
#       #geom_text(aes(dTm, -log10(p_dTm), label = delabel), data = df_,color="black")+
#       # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
#       #                 nudge_x = 1,
#       #                 force = 2,
#       #                 box.padding = 2,
#       #                 segment.alpha = .5)+
#       ggtitle(df_$sample_name[1])+
#       theme(legend.position="bottom", legend.box = "horizontal")+
#       geom_text(aes(x=10,y=max(-log10(df_$p_dTm),na.rm=TRUE)),label=paste0("S = ",nrow(df_[df_$diffexpressed=="Stabilized Shift",])),show.legend=FALSE)+
#       geom_text(aes(x=-10,y=max(-log10(df_$p_dTm),na.rm=TRUE)),label=paste0("DS = ",nrow(df_[df_$diffexpressed=="Destabilized Shift",])),show.legend=FALSE)+
#       xlim(-15,15)+
#       geom_label_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
#                        nudge_x = 1,
#                        force = 3,
#                        box.padding = 3,
#                        segment.alpha = .5) 
#     if(isTRUE(NPARC)){
#       df_$diffexpressed<-as.factor(df_$diffexpressed)
#       flevels<-data.frame(colors=c("blue","black","red"))
#       check<-ggplot(data=df_,mapping=aes(x=fStatistic,y=-log10(p_dTm),color=diffexpressed))+geom_point()+ geom_vline(xintercept=c(-2, 2), col="red") +
#         geom_hline(yintercept=-log10(0.01), col="red")+ scale_color_manual("Stabilization",values=flevels$colors[1:length(unique(df_$diffexpressed))],labels = unique(df_$diffexpressed))+
#         labs(y=expression(-log["10"]*(p["adj"]-value)),x=expression(Delta*T["m"]))+
#         geom_label(aes(fStatistic, -log10(p_dTm), label = delabel), data = df_,color="black")+
#         # geom_text_repel(df_,mapping=aes(dTm, -log10(p_dTm),label=targets),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
#         #                 nudge_x = 1,
#         #                 force = 2,
#         #                 box.padding = 2,
#         #                 segment.alpha = .5)+
#         ggtitle(df_$sample_name[1])+
#         theme(legend.position="bottom", legend.box = "horizontal")+
#         geom_text(aes(x=50,y=max(-log10(df_$p_dTm),na.rm=TRUE)),label=paste0("S = ",nrow(df_[df_$diffexpressed=="Stabilized Shift",])),show.legend=FALSE)+
#         geom_text(aes(x=50,y=max(-log10(df_$p_dTm)-1,na.rm=TRUE)),label=paste0("DS = ",nrow(df_[df_$diffexpressed=="Destabilized Shift",])),show.legend=FALSE)+
#         geom_label_repel(df_,mapping=aes(fStatistic, -log10(p_dTm),label=delabel),color="black",max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
#                          nudge_x = 1,
#                          force = 3,
#                          box.padding = 3,
#                          segment.alpha = .5) 
#       
#     }
#   }
#   
#   return(check)
#   
# }

calcP <-function(uniqueID, Tm, mpDiffs, binWidth){
  
  if (binWidth > length(Tm)){
    stop(simpleError(paste("Error in p-value computation: Assigned bin width (",binWidth,") is larger than maximum number of proteins that passed quality control.", sep="")))
  }
  
  
  i <- !is.na(Tm)
  Tm <- Tm[i]
  mpDiffs   <- mpDiffs[i]
  idsValid  <- uniqueID[i]
  
  mpNum       <- length(Tm)
  binWidthRel <- binWidth/mpNum
  binProp <- sort(unique(c(seq(1, 0, by=-binWidthRel), 0)))
  bounds  <- quantile(Tm, binProp, na.rm=TRUE)
  bins    <- .bincode(Tm, bounds, include.lowest=TRUE, right=TRUE)
  
  ## If bin with lowest values is smaller than the others, include it into the
  ## succeeding bin:
  if (sum(bins==1) < binWidth) bins[bins==1] <- 2
  
  ## Compute p-values for each bin
  pVals <- rep(NA_real_, length(bins))
  for (b in unique(bins)){
    iBin <- which(bins==b)
    pVals[i] <- pValsBin(mpDiffs[i])
  }
  
  ## Perform Benjamini-Hochberg correction (over all bins)
  pVals <- p.adjust(pVals, "BH")
  
  #plot(abs(minSlope)~ mpDiff, data=subset(mpDiffPvals, pVals<=0.05), col="blue")
  #points(abs(minSlope)~ mpDiff, data=subset(mpDiffPvals, pVals>0.05), col="red")
  
  ## Output vector
  pOut <- rep(NA_real_, mpNum)
  pOut <- pVals
  return(pOut)
} 
mpPval <- function(x, r_1, r0, r1){
  ## Compute p-value for dmp
  if(!is.na(x)){
    
    z<-ifelse(x>r0,(x-r0)/(r1-r0),
              (r0-x)/(r0-r_1))
    p <-0.5*VGAM::erfc(z/sqrt(2))
    
  } else {
    p <-NA
  }
  return(p)
}
#get p-values
pValsBin <- function(mpDiffs){
  ## Determine meltpoints beyond the 16th and 85th quantiles per bin
  dmp <- stats::quantile(mpDiffs, probs = c(0.1587, 0.5, 0.8413), na.rm=TRUE)
  r_1 <- dmp[1]
  r0 <- dmp[2]
  r1 <- dmp[3]
  pV <- apply(as.matrix(mpDiffs), MARGIN=1, FUN=mpPval, r_1=r_1, r0=r0, 
              r1=r1)
  return(pV)
}
replicate_labels<-function(x){
  if(any(names(x)=="Annotated_Sequence")&any(names(x)=="Accession")){
    x<-x%>% dplyr::group_by(Accession,Annotated_Sequence,treatment,sample_id) %>% 
      dplyr::mutate(replicate=row.names(.))
  }else if (any(names(x)=="uniqueID")&any(names(x)=="Annotated_Sequence")){
    x<-x%>% dplyr::group_by(uniqueID,Annotated_Sequence,treatment,sample_id) %>% 
      dplyr::mutate(replicate=row.names(.))
  }else if (any(names(x)=="uniqueID")){
    x<-x%>% dplyr::group_by(uniqueID,treatment,sample_id) %>% 
      dplyr::mutate(replicate=row.names(.))
  }else if (any(names(x) == "Accession")){
    x<-x%>% dplyr::group_by(Accession,treatment,sample_id) %>% 
      dplyr::mutate(replicate=row.names(.))
  }
}
rename_TPP<-function(x,temps=df.temps){#rename script data to run TPP
  string_db <- STRINGdb$new( version="10", species=9606,score_threshold=400, input_directory="")
  
  TPP_Cliff<-x %>% dplyr::rename("Condition"="treatment")
  if(any(names(x)=="I3")){
    TPP_Cliff<-TPP_Cliff %>% dplyr::rename("I"="I3")
    
    # ss <- data.frame(gene_name=proteins$preferred_name,uniqueID=as.character(AnnotationDbi::mapIds(org.Hs.eg.db, proteins$preferred_name, 'UNIPROT', 'SYMBOL')))
    # TPP_Cliff<-TPP_Cliff %>% dplyr::right_join(ss,by="uniqueID")
    TPP_Cliff$gene_name<-TPP_Cliff$uniqueID
    TPP_Cliff$gene_name<-as.character(TPP_Cliff$gene_name)
    
    TPP_Cliff<-dplyr::bind_rows(TPP_Cliff) %>% 
      dplyr::select(sample_id,Condition,Annotated_Sequence,gene_name,temp_ref,I) %>%
      distinct(.)
    data<-TPP_Cliff %>% 
      dplyr::select(sample_id,Condition,gene_name,Annotated_Sequence) %>% 
      distinct(.) %>% 
      dplyr::group_split(Condition,gene_name)
    if (.Platform$OS.type=="windows"){
      data<-parallel::mclapply(data,replicate_labels)
    }else{
      data<-parallel::mclapply(data,replicate_labels,mc.cores=future::availableCores())
    }
    
    
    data<-dplyr::bind_rows(data)
    TPP_Cliff<-TPP_Cliff %>%
      dplyr::right_join(data,by=names(data)[!names(data) %in% "replicate"])
    
    TPP_Cliff<-pivot_wider(
      TPP_Cliff,
      id_cols = NULL,
      names_from = temp_ref,
      names_prefix = "rel_fc_",
      names_sep = "_",
      names_repair = "minimal",
      values_from = I
    )
  }else{
    
    TPP_Cliff$gene_name<-TPP_Cliff$uniqueID
    TPP_Cliff$gene_name<-as.character(TPP_Cliff$gene_name)
    
    
    TPP_Cliff<-dplyr::bind_rows(TPP_Cliff) %>% 
      dplyr::select(sample_id,Condition,gene_name,temp_ref,I) %>%
      distinct(.)
    
    
    data<-TPP_Cliff %>% 
      dplyr::select(sample_id,Condition,gene_name) %>% 
      distinct(.) %>% 
      dplyr::group_split(Condition,gene_name) 
    if (.Platform$OS.type=="windows"){
      data<-parallel::mclapply(data,replicate_labels)
    }else{
      data<-parallel::mclapply(data,replicate_labels,mc.cores=future::availableCores())
    }
    
    data<-dplyr::bind_rows(data)
    TPP_Cliff<-TPP_Cliff %>%
      dplyr::right_join(data,by=names(data)[!names(data) %in% "replicate"])
    
    TPP_Cliff<-pivot_wider(
      TPP_Cliff,
      id_cols = NULL,
      names_from = temp_ref,
      names_prefix = "rel_fc_",
      names_sep = "_",
      names_repair = "minimal",
      values_from = I
      
    )
  }
  TPP_Cliff<-TPP_Cliff %>% distinct(.)
  check<-names(TPP_Cliff)
  check1<-check[str_detect(check,"[:digit:][:upper:]")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)
  #replace C or N with L and H
  check1<-str_replace(check1,"C","H")
  check1<-str_replace(check1,"N","L")
  #replace names
  check[data2]<-check1
  names(TPP_Cliff)<-check
  TPP_Cliff$Condition<-ifelse(TPP_Cliff$Condition=="vehicle","Vehicle","Treatment")
  TPP_Cliff$Experiment<-ifelse(TPP_Cliff$Condition=="Treatment",
                               paste0("Treatment_",TPP_Cliff$replicate),paste0(TPP_Cliff$Condition,"_",TPP_Cliff$replicate))
  
  TPP_Cliff$ComparisonVT1<-NA
  TPP_Cliff$ComparisonVT2<-NA
  
  TPP_Cliff$ComparisonVT1<-ifelse(TPP_Cliff$replicate==1,"x","")
  TPP_Cliff$ComparisonVT2<-ifelse(TPP_Cliff$replicate==2,"x","")
  
  check1<-check[str_detect(check,"rel_fc_[[:digit:]]+|rel_fc_[[:digit:]]+[:upper:]")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)
  
  config<-TPP_Cliff
  names(config)<-str_replace(names(config),"rel_fc_","")
  check<-c(config %>% dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2),config[data2])
  temp_ref<-str_replace(df.temps$temp_ref,"C","H")
  temp_ref<-str_replace(temp_ref,"N","L")
  
  temps<-df.temps %>% dplyr::mutate(temp_ref=temp_ref) %>% distinct(.)
  temps<-pivot_wider(temps,names_from=temp_ref,values_from=temperature)
  temps<-purrr::map_dfr(seq_len(nrow(TPP_Cliff)), ~temps)
  temp_ref<-str_replace(names(temps),"C","H")
  temp_ref<-str_replace(temp_ref,"N","L")
  names(temps)<-temp_ref
  TPP_Cliff<-cbind(TPP_Cliff,temps)
  #keep two replicates
  TPP_Cliff<-TPP_Cliff %>% dplyr::filter(ComparisonVT1=="x" | ComparisonVT2=="x")
  # TPP_Cliff$qssm<-as.integer(sample_id(0:125,nrow(TPP_Cliff),replace=TRUE))
  # TPP_Cliff$qupm<-as.integer(sample_id(4:40,nrow(TPP_Cliff),replace=TRUE))
  TPP_Cliff$qssm<-as.integer(5)
  TPP_Cliff$qupm<-as.integer(10)
  return(TPP_Cliff)
}

runTPP<-function(x,df.temps){
  
  TPP<- rename_TPP(x,df.temps)
  
  TPPconfig<-TPP %>% distinct(.)
  
  temp_ref<-str_replace(df.temps$temp_ref,"C","H")
  temp_ref<-str_replace(temp_ref,"N","L")
  hi<-TPP %>% 
    dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,temp_ref)
  hi<-hi[,1:14] %>% distinct(.)
  hi<-hi[!is.na(hi$Experiment),]
  row.names(hi)<-seq(nrow(hi))
  
  gene_names<-data.frame(gene_name=dplyr::bind_rows(TPP)$gene_name,
                         sample_id=dplyr::bind_rows(TPP)$sample_id,
                         Condition=dplyr::bind_rows(TPP)$Condition)
  
  gene_names$n<-row.names(gene_names)
  
  TPP<-TPP %>% dplyr::right_join(gene_names,by=c("gene_name","sample_id","Condition")) %>% distinct(.)
  
  TPP<-TPP %>% 
    dplyr::mutate(n=make.names(TPP$n,unique=TRUE)) 
  
  
  TPPdata<-TPP %>% 
    dplyr::select(gene_name,Experiment,qssm,qupm,n,tidyr::starts_with("rel_fc"),-tidyr::ends_with("rel_fc_NA")) %>%
    distinct(.) %>% 
    dplyr::filter(!is.na(Experiment)) %>% 
    dplyr::group_split(Experiment) %>% 
    setNames(unique(dplyr::bind_rows(.)$Experiment)) %>% 
    purrr::map(.,function(y) y %>%
                 magrittr::set_rownames(stringr::str_extract(y$n,"X")))
  
  
  resultPath<-file.path(getwd())
  #sort data
  TPPconfig<-TPP %>% dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,tidyr::starts_with("126"),tidyr::starts_with("127L"),tidyr::starts_with("127H"),tidyr::starts_with("128L"),tidyr::starts_with("128H"),tidyr::starts_with("129L"),tidyr::starts_with("129H"),tidyr::starts_with("130"),tidyr::starts_with("130H"),tidyr::starts_with("131L")) %>% distinct(.)
  TPPdata<-TPPdata[order(unique(TPPconfig$Experiment))]
  TPPdata<-lapply(TPPdata,function(x) as.data.frame(x) %>%
                    dplyr::ungroup(.) %>% dplyr::select(-n,-Experiment))
  
  #import TPPtr
  if(any(!isTRUE(TPPconfig$Experiment==names(TPPdata)))){
    TPPdata<-TPPdata[order(unique(TPPconfig$Experiment))]
  }
  if(any(!isTRUE(TPPconfig$Experiment==names(TPPdata)))){
    TPPdata<-TPPdata[order(unique(TPPconfig$Experiment))]
  }
  trData <- tpptrImport(configTable = TPPconfig, data = TPPdata)
  TRresults <- analyzeTPPTR(configTable = TPPconfig, 
                            methods = "meltcurvefit",
                            data = TPPdata, 
                            nCores = future::availableCores(),
                            resultPath = resultPath, 
                            plotCurves = FALSE,
                            normalize = FALSE)
  return(TRresults)
}

df.t <- function(n,temperatures,protein_path,sample_mapping_name=NA){
  if(!is.logical(sample_mapping_name)){
    TMT<-read_xlsx(sample_mapping_name) %>% 
      dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample_id"="Sample","sample_name"="MS_sample_number","time_point"="Time_point")  
    TMT$time_point<-as.factor(TMT$time_point)
    
  }else{
    
    TMT<-data.frame(NA)
    if(n==10){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = temperatures, stringsAsFactors = FALSE)
    }else if (n==11){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), temperature = temperatures, stringsAsFactors = FALSE)
    }else if (n == 16){
      TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C','132N','132C','133N','133C','134N'), temperature = temperatures,stringsAsFactors = FALSE)
    }
  }
  return(TMT)
}

medianPolish <- function(intensities, num_channels){
  wide <- matrix(intensities, byrow = TRUE, ncol = num_channels)
  tmp_fit <- stats::medpolish(wide, na.rm = TRUE, trace.iter = FALSE)
  tmp_fit$overall + tmp_fit$col
}

#Renaming the protein, peptide, psms



#create a function to parallelize gam fits

###############################
memory.limit(175921900000)#set for 16 GB RAM
plan(multicore,workers=future::availableCores())
options(future.globals.maxSize = 8000 * 1024^2)
##################################
#benchmark data with TPP
#hi<-purrr::map(df_norm1[2],function(x) try(runTPP(x,df.temps)))

df.temps<-df.t(11,temperatures=c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67,68))

# #Zebra
df.temps<-df.t(10,temperatures=c(34,37.3,40.6,43.9,47.2,50.5,53.8,57.1,60.4,64))
df.temps$temperature<-df.temps$temperature[order(df.temps$temperature,decreasing=TRUE)]
# #Covid
df.temps<-df.t(16,temperatures=NA,sample_mapping_name="sample_mapping.xlsx")

df_raw <- read_cetsa("E:/Zebrafish","E:/Zebrafish","_Proteins",Peptide="PSMs",Frac=TRUE,CFS=FALSE,solvent="Control",CARRIER=FALSE,rank=TRUE,sub=10,temperatures=df.temps,baseline="min",NORM="QUANTILE",keep_shared_proteins==FALSE)     
df_raw <- read_cetsa("/work/ivanovlab/figueroa-navedo.a/Scripts/Files/Zebra/Napabucasin/Trembl","/work/ivanovlab/figueroa-navedo.a/Scripts/Files/Zebra/Napabucasin/Trembl","_Proteins",Peptide="PSMs",Frac=TRUE,CFS=FALSE,solvent="Control",CARRIER=FALSE,rank=TRUE,sub=10,temperatures=df.temps,baseline="min",NORM="QUANTILE",keep_shared_proteins==FALSE)     

df_raw <- read_cetsa("~/Files/Scripts/Files/Covid","~/Files/Scripts/Files/Covid","_Proteins",Peptide=FALSE,CFS=FALSE,Frac=TRUE,solvent="AM",CARRIER=FALSE,rank=TRUE,sub=NA,temperatures=df.temps,baseline="min",NORM="QUANTILE",keep_shared_proteins==FALSE)                                                              
df_raw <- read_cetsa("~/CS7290/Protein_analysis","~/CS7290/Protein_analysis","_Proteins",Peptide="PG",CFS=TRUE,Frac=TRUE,solvent="Control",CARRIER=FALSE,rank=FALSE,sub=1000,temperatures=df.temps,baseline="min",NORM="QUANTILE",keep_shared_proteins==FALSE)                                                              

df_raw <- read_cetsa("/work/ivanovlab/figueroa-navedo.a/Scripts/Files/2.4/CFS_vs_CFE/Fractions_I/Shared",
                     "/work/ivanovlab/figueroa-navedo.a/Scripts/Files/2.4/CFS_vs_CFE/Fractions_I/Shared",
                     Prot_Pattern = "_Proteins",Peptide="PSMs",Frac=TRUE,solvent="DMSO",CARRIER=TRUE,
                     CFS=TRUE,rank=TRUE,sub=NA,temperatures=df.temps,baseline="min",NORM="QUANTILE",
                     keep_shared_proteins=FALSE)                                      


#df_raw <- read_cetsa("~/Files/Scripts/Files/CONSENSUS/Unshared","~/Files/Scripts/Files/CONSENSUS/Unshared","_Proteins",Peptide=FALSE,Frac=FALSE,solvent="DMSO",CARRIER=TRUE)                                                              
#df_raw <- read_cetsa("~/CONSENSUS11","~/CONSENSUS11","_Proteins",Peptide=FALSE,Frac=FALSE,solvent="DMSO")                                                              
#saveRDS(df_raw,"df_raw.RDS")
filter_Peptides<-function(df_,S_N,PEP,XCor,Is_Int,Missed_C,Mods,Charg,DeltaMppm,Occupancy,filter_rank=FALSE,keep_shared_proteins=FALSE,CFS=TRUE,Frac=FALSE){
  
  if(any(stringr::str_detect(names(df_),"Accession"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Accession")]
    df_<-df_ %>% dplyr::rename("Accession"=pep)
  }
  if(any(stringr::str_detect(names(df_),"Ion_Inject_Timems_"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Ion_Inject_Timems_")]
    df_<-df_ %>% dplyr::rename("IonInjTime"=pep)
  }
  if(any(stringr::str_detect(names(df_),"Isolation_Interference_"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Isolation_Interference_")]
    df_<-df_ %>% dplyr::rename("I_Interference"=pep)
  }
  if(!any(names(df_)=="sample_name")){
    df_$sample_name<-df_$Spectrum.File[1]
  }
  if(any(stringr::str_detect(names(df_),"PEP"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"PEP")]
    df_<-df_ %>% dplyr::rename("P_PEP"=pep)
  }
  if(any(stringr::str_detect(names(df_),"Charge"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Charge")]
    df_<-df_ %>% dplyr::rename("Charge"=pep)
  }
  if(any(stringr::str_detect(names(df_),"DeltaM"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"DeltaM")]
    df_<-df_ %>% dplyr::rename("DeltaM"=pep)
  }
  if(any(stringr::str_detect(names(df_),"XCorr"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"XCorr")]
    df_<-df_ %>% dplyr::rename("XCorr"=pep)
  }
  
  if(any(stringr::str_detect(names(df_),"Average_Reporter"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Average_Reporter")]
    df_<-df_ %>% dplyr::rename("AR_S_N"=pep)
  }
  if(any(stringr::str_detect(names(df_),"value"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"value")]
    df_<-df_ %>% dplyr::mutate("I"=pep)
  }
  
  if(any(names(df_)=="Fraction")&!any(names(df_)=="replicate")){
    df_$replicate<-as.numeric(df_$Fraction)
    df_$Fraction<-as.numeric(df_$Fraction)
  }
  if(any(stringr::str_detect(df_$Accession,";"))&!isTRUE(keep_shared_proteins)){
    df_<-df_[!stringr::str_detect(df_$Accession,";"),]
  }
  if(any(stringr::str_detect(names(df_),"PEP"))=="FALSE"){
    df_<-df_ %>% dplyr::mutate(Percolator_PEP=0)
  }
  #rename problematic headers
  check<-names(head(df_))
  ch<-stringr::str_replace_all(check,paste0("[","[:punct:]","]"),paste0("_"))
  ch<-stringr::str_replace_all(ch,paste0("[:punct:]","[:punct:]","[:punct:]"),paste0(""))
  ch<-stringr::str_replace_all(ch,paste0("[:punct:]","[:punct:]"),paste0(""))
  #set new names
  names(df_)<-ch
  #filter
  
  df_<-df_ %>% dplyr::filter("AR_S_N">S_N,P_PEP<PEP,Charge<Charg,MissedCleavages<Missed_C,abs(DeltaM)<DeltaMppm)
  
  
  if(any(names(df_)=="Channel_Occupancy_")){
    df_<-df_ %>% dplyr::filter(Channel_Occupancy_>Occupancy)
  }
  #rank by the highest intensity channel
  if(isTRUE(Frac)&filter_rank==TRUE){
    rank<-df_%>% dplyr::filter(temp_ref=="126") %>% 
      dplyr::mutate(rank=dplyr::ntile(.$S_N,3)) %>%
      dplyr::filter(!is.na(rank),!is.na(uniqueID)) 
    #convert to data.table
    rank<-data.table::data.table(rank)
    df_<-data.table::data.table(df_)
    
    
    #Use I column to distinguish duplicated rank values/PSM
    col_n<-dplyr::intersect(names(rank),names(df_))
    #set join columns
    data.table::setkeyv(rank,cols=col_n)
    data.table::setkeyv(df_,cols=col_n)
    #right_join
    df_ <- data.frame(merge(df_,rank, all.y=TRUE))
    
    df_<-df_ %>% dplyr::mutate(rank_l=ifelse(rank==3,TRUE,FALSE)) 
    
    df_<-df_ %>% dplyr::group_split(sample_id,uniqueID)
    df_<-df_ %>% purrr::keep(function(x) any(x$XCor_l==TRUE & x$rank_l==TRUE,na.rm=TRUE))
  }else if(!isTRUE(Frac)&filter_rank==TRUE){
    rank<-df_%>% dplyr::filter(temp_ref=="126") %>% 
      dplyr::mutate(rank=dplyr::ntile(.$S_N,3)) %>%
      dplyr::filter(!is.na(rank),!is.na(uniqueID)) 
    #convert to data.table
    rank<-data.table::data.table(rank)
    df_<-data.table::data.table(df_)
    
    
    #Use I column to distinguish duplicated rank values/PSM
    col_n<-dplyr::intersect(names(rank),names(df_))
    #set join columns
    data.table::setkeyv(rank,cols=col_n)
    data.table::setkeyv(df_,cols=col_n)
    #right_join
    df_ <- data.frame(merge(df_,rank, all.y=TRUE))
    df_<-df_ %>% dplyr::mutate(rank_l=ifelse(rank==3,TRUE,FALSE)) 
    
    df_<-df_ %>% dplyr::group_split(sample_id,uniqueID)
    df_<-df_ %>% purrr::keep(function(x) any(x$XCor_l==TRUE & x$rank_l==TRUE,na.rm=TRUE))
    
  }
  if(any(names(df_)=="value")&any(names(df_)=="Accession")){
    df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","PEP"="P_PEP")
  }
  if(length(XCor)==2){
    df_<-df_ %>% dplyr::filter(XCorr>XCor[1],XCorr<XCor[2])
  }else{
    df_<-df_ %>% dplyr::mutate(XCor_l=ifelse(Charge==2 & XCorr > 1.8,TRUE,ifelse(Charge>2 & XCorr > XCor,TRUE,FALSE)))
  }
  
  
  df_<-dplyr::bind_rows(df_) %>% distinct(.)
  
  
  
  return(df_)
  
}

####Zebra
df_raw1<-furrr::future_map(df_raw,function(x) try(filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,80,filter_rank=FALSE,keep_shared_proteins=FALSE,CFS=FALSE,Frac=TRUE)))

#######
df_raw1<-furrr::future_map(df_raw,function(x) try(filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,80,filter_rank=FALSE,keep_shared_proteins=FALSE,CFS=FALSE,Frac=TRUE)))
df_raw1<-furrr::future_map(df_raw1,function(x) try(filter_Peptides(x,20,0.01,2.3,30,2,1,7,5,80,filter_rank=FALSE,keep_shared_proteins=TRUE,CFS=TRUE,Frac=FALSE)))
df_raw1<-df_raw1 %>% purrr::keep(function(x) !nrow(x)==0)

#average peptides to the protein level
Mean_Ab<-function(x){
  
  x<-x %>% 
    dplyr::group_split(uniqueID,treatment,temp_ref,sample_id) 
  
  x<-purrr::map(x,function(x) x %>% distinct(.) %>% 
                  dplyr::mutate(I=mean(.$I,na.rm=TRUE))) 
  
  x<-dplyr::bind_rows(x) %>% dplyr::ungroup(.) %>% distinct(.)
  return(x)
  
  
}
#sum Peptides to the protein level
Sum_Ab<-function(x){
  
  x<-x %>% 
    dplyr::group_split(uniqueID,treatment,temp_ref,sample_id) 
  
  x<-purrr::map(x,function(x) x %>% distinct(.) %>% 
                  dplyr::mutate(I=sum(.$I,na.rm=TRUE))) 
  
  
  x<-dplyr::bind_rows(x) %>% dplyr::ungroup(.) %>% distinct(.)
  return(x)
}
#Add up the fractions into peptide groups
if (.Platform$OS.type=="windows"){
  df_raw_sum<-parallel::mclapply(df_raw1,Sum_Ab)
}else{
  df_raw_sum<-parallel::mclapply(df_raw1,Sum_Ab,mc.cores=future::availableCores())
}

#average up the fractions into peptide groups
if (.Platform$OS.type=="windows"){
  df_raw_mean<-parallel::mclapply(df_raw1,Mean_Ab)
}else{
  df_raw_mean<-parallel::mclapply(df_raw1,Mean_Ab,mc.cores=future::availableCores())
}

df_raw1<-df_raw1 %>% purrr::keep(function(x) !length(x)==0)
#new for PSMs
#df_raw<-df_raw %>% dplyr::rename("sample_name"="Spectrum.File")
#annotate protein data with missing values
MID<-df_raw[is.na(df_raw$value),]
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = unique(df_raw$temp_ref),temperature = c(40, 42.1, 43.8, 46.5, 50, 54, 57.3, 60.1, 62, 64), stringsAsFactors = FALSE)


# df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(67, 64, 60.4, 57.1, 53.8, 50.5, 47.2, 43.9, 40.6, 37.3), stringsAsFactors = FALSE)
#df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4'), sample_name = c('DMSO_1','DMSO_2', '655_1','655_2'), stringsAsFactors = FALSE)

df.s <- function(data_path,n,rep_,bio_,vehicle_name,treated_name,Frac=FALSE,PSM=FALSE){#n is df_raw, rep is tech rep
  
  if(any(names(n)=="sample_name") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,sample_name) %>% unique(.)
    return(n)
  }
  if(any(names(n)=="Spectrum.File") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,Spectrum.File) %>% unique(.)
    return(n)
  }
  if (!isTRUE(PSM)){
    if(isTRUE(Frac)){
      n<-data.frame(sample_name=unique(n$sample_name),sample_id=unique(n$sample_id))
      return(n)
    }
    f<-data_path
    find<-c('[:upper:][[:digit:]]+')
    check<-list()
    check<-purrr::map(seq(f),function(x){
      data.frame(sample_name= stringr::str_extract_all(names(read_xlsx(f[x]))[str_detect(names(read_xlsx(f[x])),"Abundance")],find)[[1]])
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

df.samples<-df.s(f,dplyr::bind_rows(df_raw),3,2,"DMSO","TREATED",Frac=TRUE,PSM=TRUE)
#Peptides
df_raw<-dplyr::bind_rows(df_raw) %>% dplyr::group_by(sample_name) %>% group_split(.)
#Zebra
#df_raw1<-list(dplyr::bind_rows(df_raw1))
df_clean_mean <- furrr::future_map(df_raw_mean,function(x) try(clean_cetsa(x, temperatures = df.temps,samples = df.samples,Peptide=TRUE,solvent="Control",CFS=FALSE,CARRIER=FALSE,baseline="min")))#assgns temperature and replicate values

df_clean_sum <- furrr::future_map(df_raw_sum,function(x) try(clean_cetsa(x, temperatures = df.temps,samples = df.samples,Peptide=TRUE,solvent="Control",CFS=FALSE,CARRIER=FALSE,baseline="min")))#assgns temperature and replicate values

df_raw <- furrr::future_map(df_raw,function(x) try(clean_cetsa(x, temperatures = df.temps,samples = df.samples,Peptide=TRUE,solvent="Control",CFS=FALSE,CARRIER=FALSE,baseline="min")))#assgns temperature and replicate values

#Covid data
#df_clean<-purrr::map(seq_along(df_clean),function(x) rbind(df_clean[[1]],x))
#df_clean<-dplyr::bind_rows(df_clean) %>% dplyr::group_split(sample_name)
#zebra data
#df_clean<-list(dplyr::bind_rows(df_clean))

df_norm <- furrr::future_map(df_clean,function(x) try(normalize_cetsa(x, temperatures=df.temps,Peptide=TRUE,filters=FALSE,CARRIER=TRUE,baseline="min"))) #normalizes according to Franken et. al. without R-squared filter
#normalize data
df_norm <- furrr::future_map(df_clean,function(x) try(normalize_cetsa(x, temperatures=df.temps,Peptide=FALSE,filters=FALSE,CARRIER=FALSE,baseline="min"))) #normalizes according to Franken et. al. without R-squared filter
Int_plot<-function(x,Peptide=FALSE,raw=FALSE){
  if(any(x$temp_ref=="131C")){
    x<-x %>% dplyr::filter(!temp_ref=="131C")
  }
  x<-data.frame(x)
  x$treatment<-as.factor(x$treatment)
  if(!any(stringr::str_detect(names(x),"uniqueID"))){
    x<-x %>% dplyr::rename("uniqueID"="Accession")
    x$uniqueID<-as.factor(x$uniqueID)
  }
  if(any(stringr::str_detect(x$sample_name,"S"))&any(stringr::str_detect(x$Spectrum.File,"FAIMS"))){
    x$sample_name<-stringr::str_replace(x$sample_name,"S","\u03A6")
  }
  if(any(names(x)=="I3")){
    x<-x %>% dplyr::rename("I"="I3")
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x<-x %>% dplyr::rename("I"="value")
  }
  x<-x %>% dplyr::filter(!is.na(x$I))
  if(!any(names(x)=="C")){
    x<-x %>% dplyr::rename("C"="temperature")
  }
  x$C<-as.factor(x$C)
  x$I<-as.numeric(x$I)
  if(!isTRUE(Peptide)){#if this is a protein file
    if(!isTRUE(raw)){
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Normalized intensity protein")+
        ggtitle(x$sample_name[1])+
        ylim(-0.1,5)+
        theme(legend.position="bottom")+ labs(colour = "Temperature (\u00B0C)")
    }else{
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,10000000)
    }
  }else{#if this is a peptide file
    if(!isTRUE(raw)){
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Normalized intensity")+
        ggtitle(x$sample_name[1])+
        ylim(-0.1,2)+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")
    }else if(max(x$I,na.rm=TRUE)>900000){#if this is a PSM file
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,10000000)
    }else{#if this is a PG file
      list<-ggplot2::ggplot(x,mapping=aes(x=C,y=I))+
        geom_jitter(position=position_jitter(2),alpha=0.5)+
        geom_boxplot(mapping=aes(color=C))+xlab('Temperature (\u00B0C)')+
        facet_wrap(~treatment)+
        ylab("Raw intensity")+
        ggtitle(x$sample_name[1])+
        theme(legend.position="bottom")+labs(colour = "Temperature (\u00B0C)")+
        ylim(0,2)
    }
  }
  
}
plot_I<-purrr::map(df_clean_mean,function(x) try(Int_plot(x,Peptide=TRUE,raw=TRUE)))
plot_I1<-purrr::map(df_clean_sum,function(x) try(Int_plot(x,Peptide=TRUE,raw=TRUE)))
check<-ggplot2::ggplot_build(plot_I[[1]])
y<-get_legend(check$plot)
data<-unlist(lapply(plot_I,function(x) x$labels$title))
plot_I<-plot_I[order(data)]
P<-ggpubr::ggarrange(plotlist=c(plot_I,plot_I1),ncol=1,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

P1<-ggpubr::ggarrange(plotlist=c(plot_I1),ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
pdf("Zebra_after_mean_sum_PG.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P
dev.off()
##Generate upset plots for missing value data###
df_<-df_clean %>% dplyr::rename("sample_id"="sample_id")%>% dplyr::select(-missing_pct,-value,-missing,-rank)

df_<-df_raw%>% dplyr::right_join(df_,by=c("Accession","sample_id"))


##
listUP<-upMV(df_norm,"C_F_S",plot_multiple=TRUE,PSMs=FALSE)
pdf("UpsetMV_as_Fractions_shared.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
listUP
dev.off()

df_norm_sum<-df_clean_mean

#df_norm<-purrr::map(df_norm1,function(x) x[str_detect(x$uniqueID,c("P36507","Q02750")),])

# df_norm<-dplyr::bind_rows(df_norm) %>% dplyr::group_split(uniqueID)
# df_norm<-df_norm %>% purrr::keep(function(x) nrow(x)>1)
df_norm<-purrr::map(df_norm_sum,function(x)x %>% dplyr::filter(uniqueID %in% c("P36507","Q02750")))
df_norm1<-purrr::map(df_norm1,function(x)x %>% dplyr::filter(uniqueID %in% c("P36507","Q02750")))

#Zebra
df_norm<-purrr::map(df_norm,function(x) x %>% dplyr::mutate(uniqueID=as.character(uniqueID)))
df_norm<-purrr::map(df_norm,function(x)x %>% dplyr::filter(uniqueID %in% c("Q7ZUY3",#histone
                                                                           "Q6NV46",# stat3, 
                                                                           "Q68SP3",# stat5a, 
                                                                           "A4QNT9",# stat5b, 
                                                                           "Q90Y03",# aldh1a2,
                                                                           "Q0H2G3",# aldh1a3,
                                                                           "F1QZU7",# aldh1a2.2, 
                                                                           "A2BGR9",# aldh1a2.1, 
                                                                           "F1QLV5",#nqo1,
                                                                           "F1Q7F3")#pora
                                                           
                                                           
))                                                                    #A0A1S6L757;A0A1S6L752;Q6NV46;O93599;C0SPC7;A0A0R4IGF2;Q9DDJ8
df_norm<-purrr::map(df_norm1,function(x)x %>% dplyr::filter(Accession %in% c("A0A0R4IGF2","A0A1S6L752","A0A1S6L757","Q6NV46","O93599","C0SPC7","Q7ZTS5","Q9DDJ8","Q8JFU8","C5J410","B0S789","Q8JFS5","Q9DE49"#stat3
                                                                             ,"A4QNT9","F1QWX2","Q68SP2","Q4V9B2"#stat5
                                                                             ,"Q6P0H2","A0A0R4IM15","F1QLV5","E9QEA9"#NQO1
                                                                             ,"Q90XS8","Q8QGQ1","B3DKM0","Q90Y03"#aldh1a2
                                                                             ,"Q0H2G3",#aldh1a3
                                                                             "Q499B1","F1Q7F3","Q1XB72"#PORA
                                                                             
)))
PlotTrilinear<-function(df_norm,target,df.temps,Ft,filt,Peptide=FALSE,show_results=FALSE){
  df_norm$CC<-ifelse(df_norm$treatment=="vehicle",0,1)
  if(isTRUE(Peptide) & any(names(df_norm)=="XCorr")){
    #remove columns not needed for curve fitting 
    df_norm<-df_norm %>% dplyr::select(-XCorr,-temp_ref,-Modifications,-MissedCleavages,-DeltaM,-tidyr::contains("PEP"),-Charge)
    if(any(names(df_norm)=="I_Interference")){
      df_norm<-df_norm %>% dplyr::select(-I_Interference,-IonInjTime,-S_N,-Spectrum.File)
    }
    
  }
  if(any(names(df_norm)=="I3")){
    df_norm<-df_norm %>% dplyr::mutate(I=I3)%>% dplyr::select(-I3,-I5,-I10)
  }
  ##SCRIPT STARTS HERE
  DF<-df_norm %>% dplyr::group_split(uniqueID) #split null treatment only by protein ID
  d_<-df_norm %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,treatment) #split vehicle treatment
  d_1<-df_norm %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,treatment) #split treated treatment
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
  #split treatment into equal-sized lists
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

plot<-purrr::map(df_norm,function(x) try(PlotTrilinear(x,"P36507",df.temps,Ft=FALSE,filt=FALSE,Peptide=TRUE,show_results=FALSE)))
#plot<-purrr::map(df_norm,function(x) try(PlotTrilinear(x,"P48506",df.temps,Ft=FALSE,filt=FALSE,Peptide=FALSE,show_results=FALSE)))

check<-ggplot2::ggplot_build(plot[[1]])
y<-get_legend(check$plot)
data<-order(unlist(lapply(plot,function(x) x$labels$title)))
plot<-plot[data]
P1<-ggarrange(plotlist=plot,ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plot2<-purrr::map(df_norm,function(x) try(PlotTrilinear(x,"Q02750",df.temps,Ft=FALSE,filt=FALSE,Peptide=FALSE,show_results=FALSE)))
check<-ggplot2::ggplot_build(plot2[[1]])
y<-get_legend(check$plot)
data<-order(unlist(lapply(plot,function(x) x$labels$title)))
plot2<-plot2[data]
P2<-ggarrange(plotlist=plot2,ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

#saveIDs filtered

pdf("Target_curves_trilinear_Protein_as_Fractions_shared.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P1
P2
dev.off()

#plot Number of curves
Check<-UpSet_curves(plotS2,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE)
pdf("Protein_TechrepFrac_Number_of_curves_upset_splines_.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
Check[[1]]
dev.off()

pdf("PEPFILT_Number_of_curves_splines_Peptide_filt.pdf",encoding="CP1253.enc",compress=FALSE,width=6.02,height=6.02)
Check[[2]]
dev.off()


###############################
plot_Splines<-function(x,Protein="Q02750",df.temps,Filters=FALSE,fT=TRUE,show_results=TRUE,Peptide=TRUE,CI=TRUE,simulations=FALSE,CARRIER=FALSE,Frac=TRUE,raw=FALSE){
  Filters=Filters
  fT=fT
  
  if(any(names(x)=="Charge")){
    x$Charge<-as.factor(x$Charge)
  }
  if(any(names(x)=="value")&!any(names(x)=="C")){
    x<-x %>% dplyr::rename("C"="temperature")
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x<-x %>% dplyr::rename("I"="value")
    
  }
  if(any(names(x)=="Fraction")){
    x<-x %>% dplyr::select(-Fraction) %>% distinct(.)
    x<-x %>% dplyr::mutate(Fraction=1)
    
  }
  if(any(names(x)=="Spectrum.File")){
    x<-x %>% dplyr::select(-Spectrum.File) %>% distinct(.)
    
  }
  if(any(names(x)=="Accession")&!any(names(x)=="uniqueID")){
    x<-x %>% dplyr::rename("uniqueID"="Accession")
    
  }
  
  if(length(unique(x$treatment))==1){
    return(warning("Please check that you have vehicle and treated samples in your data.  Only one of the conditions read."))
  }
  
  if(any(names(x)=="Annotated_Sequence")&any(names(x)=="Spectrum_File")){
    x<-dplyr::bind_rows(x) %>% dplyr::select(-Spectrum_File) %>%
      distinct(.) %>% 
      dplyr::group_by(uniqueID,Annotated_Sequence,temp_ref,treatment,sample_name) %>% 
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(replicate=row.names(.)))
  }else if(any(names(x)=="Annotated_Sequence")&any(names(x)=="Spectrum.File")){
    x<-dplyr::bind_rows(x) %>% dplyr::select(-Spectrum.File) %>%
      distinct(.) %>% 
      dplyr::group_by(uniqueID,Annotated_Sequence,temp_ref,treatment,sample_name) %>% 
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(replicate=row.names(.)))
  }else if(any(names(x)=="Spectrum_File")){
    x<-dplyr::bind_rows(x) %>% dplyr::select(-Spectrum_File) %>%
      distinct(.) %>% 
      dplyr::group_by(uniqueID,temp_ref,treatment,sample_name) %>% 
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(replicate=row.names(.)))
  }else if(any(names(x)=="Spectrum.File")){
    x<-dplyr::bind_rows(x) %>% dplyr::select(-Spectrum.File) %>%
      distinct(.) %>% 
      dplyr::group_by(uniqueID,temp_ref,treatment,sample_name) %>% 
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(replicate=row.names(.)))
  }
  if(isTRUE(Peptide)){
    x<-dplyr::bind_rows(x)
  }
  
  if(isTRUE(CARRIER)){
    df.temps<-df.temps %>% dplyr::filter(temperature<68)
  }
  
  if(!isTRUE(Peptide)){#if this is a protein file
    DFN<-x 
    df_<-x %>% dplyr::filter(treatment=="vehicle")
    df_1<-x %>% dplyr::filter(treatment=="treated")
    
    if(nrow(df_1)==0){
      df_1<-df_
    }
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    if(isTRUE(show_results)){
      spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
      return(spresults)
    }else{
      spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
      
    }
    if(any(class(spresults)=="list")){
      spresults<-dplyr::bind_rows(spresults)
      spresults<-spresults[!is.na(spresults$C),]
      res_sp<-spf(spresults,DFN,filters=Filters)
    }else{
      res_sp<-spf(spresults,DFN,filters=FALSE)
    }
    
    if(isTRUE(simulations)){
      ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
      return(ROCs)
    }
    
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
    
    return(Pred1)
  }else{#if this is a peptide file
    DFN<-dplyr::bind_rows(x)
    df_<-dplyr::bind_rows(x) %>% dplyr::filter(treatment=="vehicle")
    df_1<-dplyr::bind_rows(x) %>% dplyr::filter(treatment=="treated")
    
    if(nrow(df_1)==0){#if there are no treated files, make one up to avoid errors in curve fitting
      df_1<-df_
    }
    #get spline results
    spresults<-list()
    spresults_PI<-list()
    
    spresults<-tryCatch(spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters),
                        error = function(cond){
                          message("Error found in fitting and summarise function")
                          message(cond)
                        })
    
    if(isTRUE(show_results)){
      return(spresults)
    }
    if(any(class(spresults)=="list")){
      spresults<-dplyr::bind_rows(spresults)
      spresults<-spresults[!is.na(spresults$C),]
      res_sp<-spf(spresults,DFN,filters=Filters)
    }else{
      
      DFN$temperature<-DFN$C
      res_sp<-spf(spresults,DFN,filters=FALSE)
    }
    
    if(isTRUE(simulations)){
      ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
      return(ROCs)
    }
    
    #saveIDs filtered
    i<-which(res_sp[[1]]$uniqueID %in% Protein)
    #generate 95%CI for splines
    Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
    
    return(Pred1)
  }
  
  DFN<-x %>% dplyr::filter(uniqueID %in% as.character(Protein))
  df_<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),treatment=="vehicle")
  df_1<-x %>% dplyr::filter(uniqueID %in% as.character(Protein),treatment=="treated")
  #get spline results
  spresults<-list()
  spresults_PI<-list()
  
  spresults<-spstat(DFN,df_,df_1,Ftest=fT,show_results=show_results,filters=Filters)
  if(isTRUE(show_results)){
    return(spresults)
  }
  
  if(any(class(spresults)=="list")){
    spresults<-dplyr::bind_rows(spresults)
    spresults<-spresults[!is.na(spresults$C),]
    res_sp<-spf(spresults,DFN,filters=Filters)
  }else{
    DFN$temperature<-DFN$C
    res_sp<-spf(spresults,DFN,filters=Filters)
  }
  if(isTRUE(simulations)){
    ROCs<-spSim(res_sp[[1]],res_sp[[2]],res_sp[[3]])
    return(ROCs)
  }
  #saveIDs filtered
  i<-which(res_sp[[1]]$uniqueID %in% Protein)
  #generate 95%CI for splines
  Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],df.temps,overlay=TRUE,alpha=0.05,residuals=FALSE,simulations=FALSE,Peptide=Peptide,CI=CI,CARRIER=CARRIER,Protein=Protein,raw=raw)
  return(Pred1)
  
  
}

#For zebra
# "Q7ZUY3",#histone
# "Q6NV46",# stat3, 
# "Q68SP3",# stat5a, 
# "A4QNT9",# stat5b, 
# "Q90Y03",# aldh1a2,
# "Q0H2G3",# aldh1a3,
# "F1QZU7",# aldh1a2.2, 
# "A2BGR9",# aldh1a2.1, 
# "F1QLV5",#nqo1,
# "F1Q7F3"#pora

#Zebra
#"Q499B1","F1Q7F3","Q1XB72","Q0H2G3")))#GCLC
df_norm1<-list(dplyr::bind_rows(df_norm1))
df_norm<-list(dplyr::bind_rows(df_norm))
plotS2 <- purrr::map(df_norm,function(x) try(plot_Splines(x,"F1QLV5",df.temps,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=TRUE,simulations=FALSE,CARRIER=FALSE,Frac=TRUE,raw=TRUE)))
P1<-plotS2
#saveRDS(plotS2,"Napabucasin_Protein_unique_data.RDS")
#A4QNT9 F1Q7F3 F1QLV5 Q0H2G3 Q6NV46 Q90Y03
check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
#P1<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2 <- purrr::map(df_norm,function(x) try(plot_Splines(x,"P36507",df.temps,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE,simulations=FALSE,CARRIER=TRUE,Frac=FALSE,raw=FALSE)))
P3<-plotS2

check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
#ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2 <- purrr::map(df_clean,function(x) try(plot_Splines(x,"Q02750",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=TRUE,simulations=FALSE,CARRIER=TRUE,Frac=TRUE,raw=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P2<-plotS2
#P3<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2 <- purrr::map(df_clean_ch,function(x) try(plot_Splines(x,"Q02750",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=TRUE,simulations=FALSE,CARRIER=TRUE,Frac=TRUE,raw=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P4<-plotS2
P5<-ggarrange(plotlist=c(P1,P3),ncol=2,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
#Zebra
df_norm1<-list(dplyr::bind_rows(df_norm1))
plotS <- furrr::future_map(df_norm,function(x) try(plot_Splines(x,"Q0H2G3",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE,simulations=FALSE,CARRIER=FALSE,Frac=TRUE)))
check<-ggplot2::ggplot_build(plotS[[1]])
y<-get_legend(check$plot)
data<-unlist(lapply(plotS,function(x) x$labels$title))
plotS<-plotS[order(data)]
P3<-ggarrange(plotlist=plotS,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
#For covid treatment
# check<-dplyr::bind_rows(df_norm) %>% dplyr::group_split(time_point)
# plotS2 <- purrr::map(check,function(x) try(plot_Splines(x,"P0DTC2",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE)))

plotS2 <- purrr::map(df_norm,function(x) try(plot_Splines(x,"Q0H2G3",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=FALSE,simulations=FALSE,CARRIER=FALSE)))
saveRDS(plotS2,"CFS_repPeptide_filt_data.RDS")
check<-ggplot2::ggplot_build(plotS2[[2]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P2<-ggarrange(plotlist=plotS2,ncol=4,nrow=4,font.label = list(size = 14, color = "black", face = "bold"),labels = c("I","J","K","L","M","N","O","P"),legend.grob = y)
pdf("CFE_vs_CFS_protein_as_fractions_shared_TPPquant.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P1
P2
dev.off()

plotS <- furrr::future_map(df_norm,function(x) try(plot_Splines(x,"P36507",df.temps,MD=TRUE,Filters=FALSE,fT=FALSE,show_results=FALSE,Peptide=TRUE,simulations=FALSE)))
check<-ggplot2::ggplot_build(plotS[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS,function(x) x$labels$title))
# plotS<-plotS[order(data)]
P1<-ggarrange(plotlist=plotS,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2 <- purrr::map(df_norm1,function(x) try(plot_Splines(x,"P36507",df.temps,MD=TRUE,Filters=FALSE,fT=TRUE,show_results=TRUE,Peptide=TRUE,simulations=FALSE)))
check<-ggplot2::ggplot_build(plotS2[[7]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P3<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
saveRDS(plotS2,"PEPTIDE_UNFILTERED_NO_IMAATSA.RDS")

pdf("Peptide_data_Shared_MEK2_as_replicates.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P2
dev.off()

plotS2 <- purrr::map(df_norm1,function(x) try(plot_Splines(x,"Q02750",df.temps,MD=TRUE,Filters=FALSE,fT=TRUE,show_results=FALSE,Peptide=TRUE,simulations=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[1]][[1]])
y<-get_legend(check$plot)

real<-purrr::map(plotS2,function(x) x[[1]])
sim<-purrr::map(plotS2,function(x) x[[2]])
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P2<-ggarrange(plotlist=real,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
P1<-ggarrange(plotlist=sim,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("PSM_MEK2_Cliff_curves_Before_and_after_int>70k_filter_and_SN>10_shared.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
P5
dev.off()
#plot volcano for unfiltered data
check<-purrr::map(plotS2,function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,benchmark=TRUE,Fhist=FALSE,labels=TRUE,type="other")))
check_<-ggplot2::ggplot_build(check[[1]])#fT is tied to benchmark
y<-get_legend(check_$plot)

P1<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

check2<-purrr::map(plotS2,function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,benchmark=TRUE,Fhist=TRUE,labels=TRUE,type="other")))
check_<-ggplot2::ggplot_build(check2[[1]])
y<-get_legend(check_$plot)
P2<-ggarrange(plotlist=check2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

check3<-purrr::map(plotS2,function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,benchmark=NA,Fhist=FALSE,labels=TRUE,type="other")))
pdf("PROTEIN_Napabucasin_Zebra_unique_volcano_splines_hist_targets.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
check
check2
check3
dev.off()
#Volcano plots using Z-scores
check2<-purrr::map(plotS2,function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,benchmark=FALSE,Fhist=FALSE,labels=FALSE,type="targets")))
check_<-ggplot2::ggplot_build(check2[[1]])
y<-get_legend(check_$plot)
P2<-ggarrange(plotlist=check2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("volcano_splines_PEPTIDE_filtered_zscores_targets.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P2
dev.off()
#plot volcano for gof filtered data (F statistic)
check<-purrr::map(plotS2[1],function(x) try(volcano_data(x,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=TRUE,fT=FALSE)))
check1<-ggplot2::ggplot_build(check[[1]])
y<-get_legend(check1$plot)
P2<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("volcano_splines_Peptide_filtered_IMAATSA_gof_panels.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P2
dev.off()

check<-TPPbenchmark_generic(f,volcano=FALSE,filters="HQ",Peptide=TRUE,filter=FALSE)

pdf("TPP_HQ_Peptide_filt_techreps_fractions.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
check1
dev.off()
#ot Number of curves
Check<-UpSet_curves(plotS2,Trilinear=FALSE,Splines=TRUE,Sigmoidal=FALSE,Peptide=FALSE,filter=FALSE)

pdf("Number_of_curves_upset_splines_PEPTIDE_filtered.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
Check
dev.off()
###################################################                                                                                                                                                                                                                                                                                                                            ##################################################
#Sigmoidal function with confidence intervals
###################################################

df_<-df_norm

plot_Sigmoidal<-function(df_,Protein,Peptide=FALSE){
  #remove duplicated intensity values
  df_<-df_ %>% distinct(.)
  df_$sample_name<-str_replace(df_$sample_name,"S","\u03A6")
  
  PlSig<-try(sigC(df_,Protein,Peptide=Peptide))
  
  ID<-unique(PlSig$uniqueID)
  
  sig<- try(sigfit(PlSig,Peptide=Peptide))
  return(sig)
}
plotS2<-ggplot()
plotS2<-purrr::map(df_,function(x) try(plot_Sigmoidal(x,"P36507",Peptide=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[1]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P2<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

plotS2<-ggplot()
plotS2<-purrr::map(df_,function(x) try(plot_Sigmoidal(x,"Q02750",Peptide=TRUE)))
check<-ggplot2::ggplot_build(plotS2[[3]])
y<-get_legend(check$plot)
# data<-unlist(lapply(plotS2,function(x) x$labels$title))
# plotS2<-plotS2[order(data)]
P1<-ggarrange(plotlist=plotS2,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)

pdf("Sigmoidal_curves_Peptide_MD_PD_HR_Settings_unfilt_2.pdf",encoding="CP1253.enc",compress=TRUE,width=12.13,height=7.93)
P1
P2
dev.off()




