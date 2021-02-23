library(tidyverse)
library("writexl")
library(dplyr)
library(MSstatsTMT)
library(stringr)
library(stringi)
theme_set(theme_bw())


#import data
f<- list.files(pattern='*PSMs.xlsx')
PSMs <-  "~/Cliff prot pep/PSMs.xlsx"
PSMs <- read_excel(f)

#main header names change code
Changecolname <- function(non_edited_PSMs){
  editedPSMsfile <-non_edited_PSMs %>% 
    rename_at(vars(starts_with(c("12"))), dplyr::funs(str_replace(., "^", "Abundance.."))) %>% 
    rename_at(vars(starts_with(c("13"))), dplyr::funs(str_replace(., "^", "Abundance.."))) %>% 
    dplyr::rename_all(list(~make.names(.)))
  names(editedPSMsfile)<-str_replace(names(editedPSMsfile)," ",".")
  editedPSMsfile
}
editPSMs <- Changecolname(PSMs)
####

###############################################################
#Partition file names to generate other columns
Pattern <- function (editedPSMsfile,TMT_n=10,solvent,treatment,sep){
  TMT6<-data.frame(TMT=NA)
  TMT10<-data.frame(TMT=NA)
  TMT11<-data.frame(TMT=NA)
  ###################################################################### 
  #get all possible TMT channels
  TMT6<- data.frame(TMT=c("126","127","128","129","130","131"))
  TMT10<- data.frame(TMT = c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131"))
  TMT11<- data.frame(TMT = c("126","127N","127C","128N","128C", "129N", "129C","130N", "130C", "131N","131C"))
  #allocate space for the data
  TMT<-data.frame(TMT=rep(as.character(NA),TMT_n))
  #get user input to assign TMT channels
  TMT<-ifelse(TMT_n==10,TMT %>% dplyr::mutate(TMT=TMT10$TMT),ifelse(TMT_n==11,TMT %>% dplyr::mutate(TMT=TMT11$TMT),TMT %>% dplyr::mutate(TMT=TMT6$TMT)))
  TMT<-data.frame(TMT=TMT[[1]])
  #Allocate memory in a data frame
  Annotation3<-data.frame(Run = NA,shortnames=NA,Mixture=NA,TechRepMixure=NA,BioReplicate=NA,Pattern1=NA,Channel=NA)
  #start filling out columns
  Annotation3<-data.frame(Run = unique(editedPSMsfile$Spectrum.File))
  if(is.null(treatment)&is.null(sep)){
  pattern1<-stri_match_first_regex(Annotation3$Run,"[:upper:][:digit:][:digit:][:punct:]|[:upper:][:digit:][:punct:]")[,1]
  RunID <- str_split(Annotation3$Run,pattern1,simplify = TRUE)[,1]
  }else{
    Control<-str_detect(Annotation3$Run,str_c(sep,solvent,sep))
    pattern1<-ifelse(Control==TRUE,solvent,treatment)
    pattern2<-stri_match_first_regex(Annotation3$Run,str_c(pattern1,"[:punct:][:digit:][:punct:]"))[,1]
    RunID <- pattern2#pattern2 has the name of the treatment (DMSO=none,655=inhibitor,1 and 2 are biological replicates)
  }
  Annotation3$shortnames<- RunID 
  
  #check number of times data appears per run ID
  check<-Annotation3 %>% dplyr::group_by(shortnames) %>% dplyr::mutate(N = dplyr::n())
  check<-check %>% dplyr::filter(N>1) #filter out group sizes of 1
  
  Annotation3<-Annotation3 %>% subset(shortnames %in% unique(check$shortnames))
  #get the letter patterns 
  pattern2<-stri_match_first_regex(Annotation3$Run,"[:upper:][:digit:][:digit:][:punct:]|[:upper:][:digit:][:punct:]")[,1]
 #pattern2 has the names of the fractions
  ##########
  RunID_uniqe <- unique(Annotation3$shortnames)
  RunID_uniqe <-data.frame(shortnames=RunID_uniqe )
  size_RunID <-nrow(RunID_uniqe)
  
  Fra_UI<- vector() 
  for (i in 1:size_RunID) {
    print ( RunID_uniqe[i,1])
    Fra_UI[i]<-readline(prompt="What is the Sample ID from PD, please use numbers 0,1,2...?")
    
    #F1="DMSO_1"
    #F2="DMSO_2"
    #F3="655_1"
    #F4="655_2"
  }
  Fra_UI <-  data.frame(shortnames=RunID_uniqe,sample_id=(do.call(paste0, expand.grid(factor(c('F'),levels = c('F')),Fra_UI))))
  
  Annotation3 <-dplyr::left_join(Annotation3,Fra_UI,by="shortnames")
  Annotation3 <- Annotation3 %>% dplyr::mutate(Pattern1=pattern1)
  
  #order fractions
 
 
  Annotation3<-dplyr::bind_rows(Annotation3)
  Frac<-stringi::stri_extract_last(Annotation3$Run,regex="[:punct:][:upper:][[:digit:]]+")
  Frac<-str_replace(Frac,"_","")
  Frac<-data.frame(Frac=as.factor(Frac))
  Annotation3<-cbind(Annotation3,Frac)
  Annotation3<- Annotation3 %>% dplyr::group_split(shortnames,sample_id)
  Annotation3<-lapply(Annotation3,function(x)x[order(x$Frac),])
  Fractions<-lapply(Annotation3,function(x)
                    data.frame(Fraction=seq(unique(x$Frac)),
                        Frac=unique(x$Frac)))
  
  Annotation3<-dplyr::bind_rows(Annotation3)
  Fractions<-Fractions[[1]]
  
  Annotation3<-Annotation3 %>% dplyr::right_join(Fractions,by="Frac") %>% dplyr::select(-Frac)
   #split
  Annotation3<- Annotation3 %>% dplyr::group_split(shortnames,sample_id)
  
  #generate BioReplicate information
  Annotation3<- lapply(Annotation3,function(x) x %>% dplyr::mutate(BioReplicate=ifelse(stringr::str_detect(x$shortnames,solvent),"vehicle","treated"),
                                                                   TechRepMixture=ifelse(stringr::str_detect(x$shortnames,'[:punct:]1[:punct:]'),1,2),
                                                                   ))
  Annotation3<-dplyr::bind_rows(Annotation3)
   #split
  Annotation3<- Annotation3 %>% dplyr::group_split(shortnames,BioReplicate,Fraction)
  #add fractions
  Annotation3<- lapply(Annotation3,function(x) x %>% dplyr::mutate(samplefractions=row.names(.)))
  #repeat fractions 10 times because all fractions were labeled with TMT
  Annotation3<-lapply(Annotation3,function(x) x[rep(seq_len(nrow(x)), each = 10), ]) 
  #now group and label TMT
  
  Annotation3<-lapply(Annotation3,function(x) cbind(x,TMT10))
  Annotation3<-lapply(Annotation3,function(x) x %>% 
                      dplyr::rename("Channel"="TMT"))
  #generate Condition nformation
  Annotation3 <-lapply(Annotation3,function(x) x  %>% dplyr::mutate(Condition = ifelse(Channel == "126", "Norm",ifelse(stringr::str_detect(x$shortnames,solvent),0,1))))
  #generate Mixture information
  Annotation3<-lapply(Annotation3,function(x) x %>% dplyr::mutate (Mixture=ifelse(stringr::str_detect(x$shortnames,solvent),paste(solvent,as.character(x$TechRepMixture)),paste(x$shortnames[!stringr::str_detect(x$shortnames,solvent)],as.character(x$TechRepMixture)))))
  #generate  final table
  Annotation3<-do.call(rbind.data.frame, Annotation3)
  return(Annotation3)
}


#call function
Annotation<- Pattern(editPSMs,10,"DMSO","655","_")#Enter 3 2 4 1 for the prompts
#Annotation$Mixture=Annotation$shortnames
#Annotation <- Annotation %>%dplyr::mutate(BioReplicate = ifelse(Channel == "126 "& Condition=="Norm","vehicle_Norm",ifelse(Channel == "126"& Condition=="1","treated_Norm",ifelse(Condition=="0","vehicle","treated"))))
#Annotation <- Annotation %>%dplyr::mutate(BioReplicate = ifelse(Fraction=="F1"| Fraction=="F2","vehicle","treated"))
#Annotation <- Annotation %>%dplyr::mutate(Bioreplicate=if(Bioreplicate=="vehicle"&Condition=="Norm","vehicle_Norm")

#remove empty run
editPSMs<-subset(editPSMs, Spectrum.File!="cpQEX190924_A549_DMSO_1_A10_eFT_30K_22min.raw")
#partition input data by solvents
#solvents <- unique(Annotation$Run[!str_detect(Annotation$Run,"DMSO")])
#editPSMs<-editPSMs %>% dplyr::filter(Spectrum.File %in% c(solvents))
#Annotation<-Annotation %>% dplyr::filter(Run %in% c(solvents))
########################################
#Masstats TMT
Norm_input<- PDtoMSstatsTMTFormat(
  editPSMs,
  Annotation,
  which.proteinid = "Protein.Accessions",
  useNumProteinsColumn = FALSE,
  useUniquePeptide = TRUE,
  rmPSM_withMissing_withinRun = TRUE,
  rmPSM_withfewMea_withinRun = TRUE,
  rmProtein_with1Feature = FALSE,
  summaryforMultipleRows = sum
)
#normalization with global normalization and imputation
finalNorm_global_imp<-proteinSummarization(
  Norm_input,
  method = "msstats",
  global_norm = TRUE,
  reference_norm = TRUE,
  remove_norm_channel = FALSE,
  remove_empty_channel = TRUE,
  MBimpute = TRUE,
  maxQuantileforCensored = NULL
)
#normalization with global normalization without imputation
finalNorm_global_noimp<-proteinSummarization(
  Norm_input,
  method = "msstats",
  global_norm = TRUE,
  reference_norm = TRUE,
  remove_norm_channel = FALSE,
  remove_empty_channel = TRUE,
  MBimpute = FALSE,
  maxQuantileforCensored = NULL
)
##################################
#normalization without global normalization without imputation ########Works
norm_no_global_no_imp<-proteinSummarization(
  Norm_input,
  method = "MedianPolish",
  global_norm = FALSE,
  reference_norm = TRUE,
  remove_norm_channel = FALSE,
  remove_empty_channel = TRUE,
  MBimpute = FALSE,
  maxQuantileforCensored = 0.6
)
#normalization without global normalization with imputation #########Works
norm_no_global_imp<-proteinSummarization(
  Norm_input,
  method = "MedianPolish",
  global_norm = FALSE,
  reference_norm = TRUE,
  remove_norm_channel = FALSE,
  remove_empty_channel = TRUE,
  MBimpute = TRUE,
  maxQuantileforCensored = 0.6
)
plot_Norm<-dataProcessPlotsTMT(
  data.peptide=Norm_input,
  data.summarization=finalNorm_noglobal_imp,
  type='QCPlot',
  ylimUp = FALSE,
  ylimDown = FALSE,
  x.axis.size = 10,
  y.axis.size = 10,
  text.size = 4,
  text.angle = 90,
  legend.size = 7,
  dot.size.profile = 2,
  ncol.guide = 5,
  width = 10,
  height = 10,
  which.Protein = "all",
  originalPlot = TRUE,
  summaryPlot = TRUE,
  address = ""
)

finalNorm_noglobal_imp<-proteinSummarization(
  Norm_input,
  method = "msstats",
  global_norm = FALSE,
  reference_norm = TRUE,
  remove_norm_channel = FALSE,
  remove_empty_channel = TRUE,
  MBimpute = TRUE,
  maxQuantileforCensored = NULL
)
write_xlsx(finalNorm_global_imp,"~/RMSStatsTMT/finalNorm_global_imp.xlsx") 
write_xlsx(finalNorm_global_noimp,"~/RMSStatsTMT/finalNorm_global_no_imp.xlsx") 
write_xlsx(finalNorm_noglobal_imp,"~/RMSStatsTMT/MSSTATSTMT_logS_no_global_imp.xlsx") ##This worked
write_xlsx(norm_no_global_no_imp,"~/RMSStatsTMT/norm_no_global_no_imp.xlsx")
