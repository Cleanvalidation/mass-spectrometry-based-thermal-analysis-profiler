#import packages
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
theme_set(theme_bw())

#import data
PSMs <-  "~/CS7290/PSMs.xlsx"
PSMs<-readxl::read_excel(PSMs)
#Reaname column name according to Diffrent version #PSMs file:

#main header names change code
Changecolname <- function(non_edited_PSMs){
  editedPSMsfile <-non_edited_PSMs %>% 
    rename_at(vars(starts_with(c("12"))), dplyr::funs(str_replace(., "^", "Abundance.."))) %>% 
    rename_at(vars(starts_with(c("13"))), dplyr::funs(str_replace(., "^", "Abundance.."))) %>% 
    dplyr::rename_all(list(~make.names(.)))
  editedPSMsfile
}
editPSMs <- Changecolname(PSMs)
####

###############################################################
#Partition file names to generate other columns
Pattern <- function (editedPSMsfile){
  Annotation3<-data.frame(Run=(unique(editedPSMsfile$Spectrum.File)))
  
  RunID <- str_split(Annotation3$Run,pattern1,simplify = TRUE)[,1]
  #Annotation3<-data.frame(Annotation3)
  Annotation3$shortnames<- RunID 
  
  #check number of times data appears per run ID
  check<-Annotation3 %>% dplyr::group_by(shortnames) %>% dplyr::mutate(N = dplyr::n())
  check<-check %>% dplyr::filter(N>1) 
#filter out group sizes of 1
  Annotation3<-Annotation3 %>% subset(shortnames %in% unique(check$shortnames))
  pattern1<-stri_match_first_regex(Annotation3$Run,"[:upper:][:digit:][:digit:][:punct:]|[:upper:][:digit:][:punct:]")[,1]
  ##########
  RunID_uniqe <- unique(Annotation3$shortnames)
  RunID_uniqe <-data.frame(shortnames=RunID_uniqe )
  names(RunID_uniqe )<-"shortnames"
  size_RunID <-nrow(RunID_uniqe)
  
  Fra_UI<- vector() 
  for (i in 1:size_RunID) {
    print ( RunID_uniqe[i,1])
    Fra_UI[i]<-readline(prompt="What is the Sample ID from PD, please use numbers 0,1,2...")
    
    #F1="DMSO_1"
    #F2="DMSO_2"
    #F3="655_1"
    #F4="655_2"
  }
  Fra_UI <-  data.frame(shortnames=RunID_uniqe,Sampleid=(do.call(paste0, expand.grid(factor(c('F'),levels = c('F')),Fra_UI))))
  Annotation3 <-dplyr::left_join(Annotation3,Fra_UI,by="shortnames")
  Annotation3 <- Annotation3 %>% dplyr::mutate(Pattern1=pattern1)
  
  #generate Bioreplicate column
  Annotation3 <-Annotation3  %>% dplyr::mutate(BioReplicate=stri_match_first_regex(Annotation3$shortnames,"[:punct:][:digit:][:punct:]") )
  Annotation3 <-Annotation3  %>% dplyr::mutate(BioReplicate=stri_match_first_regex(Annotation3$BioReplicate,"[:digit:]") ) 
  Annotation3
}
Annotation3 <- suppressWarnings(Pattern(editPSMs))  


  ###################################################################### 
  #need edit?
  TMT6=c("126","127","128","129","130","131")
  TMT10= c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131")
  TMT11=c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131","131c")
 
  #order fractions
  Annotation3 <- Annotation3[order(Annotation3[,"Pattern1"]),]
  #split
  Annotation3<- Annotation3 %>% dplyr::group_split(shortnames,BioReplicate,Sampleid)
  #add fractions
  Annotation3<- lapply(Annotation3,function(x) x %>% dplyr::mutate(Fraction=row.names(.)))
  #repeat fractions 10 times because all fractions were labeled with TMT
  Annotation3<-lapply(Annotation3,function(x) x[rep(seq_len(nrow(x)), each = 10), ]) 
  #now group and label TMT
  Annotation3<-lapply(Annotation3,function(x) x %>% dplyr::mutate(Channel=rep(TMT10,nrow(x)/length(TMT10))))
  
  ###################################################################
  #Generating the Data frame consist of Annotation file columns
  Annotation3$Mixture <- ""
  Annotation3$TechRepMixture <- ""
  #Annotation3$Channel <- ""
  Annotation3$BioReplicate <- ""
  Annotation3$Condition<- "" 
  
  
  
  
  #generate Mixture information
  Annotation3<-Annotation3 %>% dplyr::mutate (Mixture="Empty")
  
  #generate Techrepmicture information
  Annotation3<- Annotation3 %>% dplyr::mutate (TechRepMixture=1)
  
  #rearrange column Run for easier generating the channel's column
  Annotation3 <- Annotation3[order(Annotation3[,1]),]
  
  #generate  TMT10channel(it needs edit)
  if(N==10){
    TMT=TMT10
  }else(N==6){
    TMT=TMT6
  }
  
  TMT=ifelse(N==6,TMT6,TMT11)
  TMT=ifelse(N==10,c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131"), ifelse(N==6,TMT6,c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131","131c")))  
  
  
  TMT6=data.frame("126","127","128","129","130","131")
  TMT10= c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131")
  TMT11=c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131","131c")
  N2<-NROW(unique(Annotation$Run))
  #Annotation <-Annotation %>% dplyr::mutate (Channel = rep(TMT10,N2))
  
  
  #generate condition column?
  Annotation3 <-Annotation3  %>% dplyr::mutate(Condition = ifelse(Channel == 126, "Norm",ifelse(Fraction == "F1" | Fraction == "F2",0,1)) )
  
  
  #########################################################################
  #clean data (it needs edit)
  nonID<-  Fra_userinput %>% filter_all(any_vars(. %in% "F0"))
  PSMsnew<-subset(editedPSMsfile, Spectrum.File==nonID)
  Annotationnew<-subset(Annotation3,Fraction==nonID)
  
  