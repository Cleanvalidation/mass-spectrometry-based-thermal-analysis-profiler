#this simulation renames data from TPP to MSStatsTMT format
MSstats_format_TO_TPP<-function(Sim,df.temps,CARRIER=TRUE){
  x<-Sim
  if(isTRUE(CARRIER)){
    if(any(stringr::str_detect(df.temps$temp_ref,"131C"))){
      df.temps<-df.temps[!stringr::str_detect(df.temps$temp_ref,"131C"),]
    }
    
  }
  df.temps<-df.temps %>% dplyr::filter(df.temps$temp_ref %in% unique(x$Channel))
  x<-x %>% dplyr::filter(x$Channel %in% unique(df.temps$temp_ref))
  #rename_TPP
  if(any(names(x)=="Protein")&!any(names(x)=="uniqueID")){
    x$uniqueID<-x$Protein
  }
  if(any(names(x)=="uniqueID")&!any(names(x)=="gene_name")){
    x$gene_name<-x$uniqueID
  }
  if(any(names(x)=="Protein")&!any(names(x)=="gene_name")){
    x$gene_name<-x$Protein
  }
  if(any(names(x)=="treatment")&!any(names(x)=="Condition")){
    x$Condition<-x$treatment
  }
  if(any(names(x)=="dataset")&!any(names(x)=="Condition")){
    x$Condition<-x$dataset
  }
  if(!any(names(x)=="uniqueID")&any(names(x)=="Protein")){
    x$uniqueID<-x$Protein
  }
  if(any(names(x)=="Channel")&!any(names(x)=="temp_ref")){
    x$temp_ref<-x$Channel
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x$I<-x$value
  }
  if(any(names(x)=="Abundance")&!any(names(x)=="I")){
    x$I<-x$Abundance
  }
  if(any(stringr::str_detect(x$Condition,"_"))){
    x$Condition<-stringr::str_extract(stringr::str_to_lower(x$Condition),"[[:lower:]]+")
  }
  if(any(names(x)=="Subject")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Subject
  }
  if(any(names(x)=="Mixture")&!any(names(x)=="Experiment")){
    x$Experiment<-x$Mixture
  }
  if(any(names(x)=="Mixture")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Mixture
  }
  if(any(names(x)=="Experiment")){
    TechReps<-data.frame(Experiment=as.factor(unique(x$Experiment)),TechRepMixture=rep(c(1,2),2))
    x$Experiment<-as.factor(x$Experiment)
    if(any(names(x)=="TechRepMixture")){
      x <-x %>% dplyr::select(-TechRepMixture)
    }
    x<-x %>% dplyr::inner_join(TechReps,by="Experiment")
  }else if(any(names(x)=="Subject")){
    TechReps<-data.frame(Experiment=as.factor(unique(x$Subject)),TechRepMixture=rep(c(1,2),2))
    x<-x %>% dplyr::rename(Experiment=Subject)
    if(any(names(x)=="TechRepMixture")){
      x <-x %>% dplyr::select(-TechRepMixture)
    }
    x<-x %>% dplyr::inner_join(TechReps,by="Experiment")
  }else if(any(names(x)=="Condition")){
    if(any(stringr::str_detect(x$Condition,"_"))){
      x$Experiment<-stringr::str_extract(stringr::str_to_lower(x$Condition),"[[:lower:]]+")
    }
  }


#stopifnot(all(stringr::str_detect(names(x),c("uniqueID","gene_name","sample_id","TechRepMixture","I","temp_ref"))))
TPP_Cliff<-dplyr::bind_rows(x) %>% dplyr::select(uniqueID,gene_name,sample_id,TechRepMixture,I,temp_ref,Experiment) %>% dplyr::filter(!is.na(I))
TPP_Cliff_pivot<-pivot_wider(
  TPP_Cliff,
  id_cols = NULL,
  names_from = temp_ref,
  names_prefix = "rel_fc_",
  names_sep = "_",
  names_repair = "minimal",
  values_from = I,
  values_fn=unique
  
)

TPP_Cliff<-TPP_Cliff_pivot %>% distinct(.)
check<-names(TPP_Cliff)
check1<-check[stringr::str_detect(check,"[:digit:][:upper:]")]
#column numbers that have reporter ion data
data2<-which(check %in% check1)
#replace C or N with L and H
check1<-stringr::str_replace(check1,"C","H")
check1<-stringr::str_replace(check1,"N","L")
#replace names
check[data2]<-check1
names(TPP_Cliff)<-check
x1<-x %>% dplyr::select(uniqueID,gene_name,Condition,sample_id,TechRepMixture)
TPP_Cliff<-TPP_Cliff %>% dplyr::inner_join(x1)
TPP_Cliff$treatment<-stringr::str_extract(str_to_lower(TPP_Cliff$Condition),"[[:lower:]]+")
TPP_Cliff$Condition<-ifelse(TPP_Cliff$treatment=="vehicle","Vehicle","Treatment")
TPP_Cliff$Experiment<-ifelse(TPP_Cliff$Condition=="Treatment",
                             paste0("Treatment_",TPP_Cliff$TechRepMixture),paste0(TPP_Cliff$Condition,"_",TPP_Cliff$TechRepMixture))

TPP_Cliff$ComparisonVT1<-NA
TPP_Cliff$ComparisonVT2<-NA

TPP_Cliff$ComparisonVT1<-ifelse(TPP_Cliff$TechRepMixture==1,"x","")
TPP_Cliff$ComparisonVT2<-ifelse(TPP_Cliff$TechRepMixture==2,"x","")

check1<-check[stringr::str_detect(check,"rel_fc_[[:digit:]]+|rel_fc_[[:digit:]]+[:upper:]")]
#column numbers that have reporter ion data
data2<-which(check %in% check1)

config<-dplyr::bind_rows(TPP_Cliff)
names(config)<-stringr::str_replace(names(config),"rel_fc_","")
check<-c(config %>% dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2),config[data2])
temp_ref<-stringr::str_replace(df.temps$temp_ref,"C","H")
temp_ref<-stringr::str_replace(temp_ref,"N","L")

temps<-df.temps %>% dplyr::mutate(temp_ref=temp_ref) %>% distinct(.)
temps<-pivot_wider(temps,names_from=temp_ref,values_from=temperature)
temps<-purrr::map_dfr(seq_len(nrow(TPP_Cliff)), ~temps)
temp_ref<-stringr::str_replace(names(temps),"C","H")
temp_ref<-stringr::str_replace(temp_ref,"N","L")

names(temps)<-temp_ref
TPP_Cliff<-cbind(TPP_Cliff,temps)
#keep two replicates
TPP_Cliff<-TPP_Cliff %>% dplyr::filter(ComparisonVT1=="x" | ComparisonVT2=="x")
TPP_Cliff$qssm<-as.integer(5)
TPP_Cliff$qupm<-as.integer(10)

temp_ref<-stringr::str_replace(df.temps$temp_ref,"C","H")
temp_ref<-stringr::str_replace(temp_ref,"N","L")

TPPconfig<-TPP_Cliff %>%
  dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,tidyr::starts_with("1")) %>%
  distinct(.)%>%
  dplyr::mutate(Experiment=as.character(Experiment)) %>% 
  dplyr::arrange(Experiment)  

TPPdata<-TPP_Cliff %>% 
  dplyr::select(uniqueID,gene_name,Experiment,qssm,qupm,tidyr::starts_with("rel_fc")) %>%
  distinct(.) %>% 
  dplyr::filter(!is.na(Experiment)) %>%
  dplyr::arrange(Experiment) %>% 
  dplyr::mutate(Experiment=as.character(Experiment)) %>% 
  dplyr::group_split(gene_name)
TPPdata<-TPPdata %>% purrr::keep(function(x) length(unique(x$Experiment))==4)

TPPdata<-dplyr::bind_rows(TPPdata) %>% dplyr::group_split(Experiment)
names(TPPdata)<-unique(TPPconfig$Experiment)
resultPath<-file.path(getwd())

TPPdata<-purrr::map(TPPdata,function(x) x %>% dplyr::filter(!is.na(rel_fc_126)))
TPPdata1<-purrr::map(TPPdata,function(x) x %>% dplyr::filter(!is.na(gene_name)) %>% reshape2::melt())
check<-purrr::map2(TPPdata,TPPdata1,function(x,y)y %>% dplyr::right_join(x))
# 
# TRreqs<-tpptrDefaultNormReqs()
# TRreqs$fcRequirements$thresholdLower<-c(-1,-1,-1)
# TRreqs$fcRequirements$thresholdUpper<-c(Inf,Inf,Inf)
# TRreqs$fcRequirements<-TRreqs$fcRequirements[1:2,]
# 
check<-lapply(check,function(x) x %>% dplyr::mutate(uniqueID=as.factor(uniqueID),
                                                    gene_name=uniqueID))
out<-list(TPPconfig,check)
names(out)<-c("TPPconfig","TPPdata")
return(out)

}

runTPP_splineFtest<-function(normData,DF=future::availableCores()*0.7){
  
  tpptrData<-TPP::tpptrImport(configTable = normData$TPPconfig,data=normData$TPPdata)
  
  
  check_longTable <- TPP::tpptrTidyUpESets(tpptrData)
  #TPP fit splines
  
  SimSplineFits <- TPP::tpptrFitSplines(data = check_longTable, 
                                        factorsH1 = c("condition"),
                                        splineDF = DF,
                                        nCores = 6)
  
  SimFtest<-tpptrFTest(SimSplineFits)
  return(SimFtest)
}
runTPP_sigmoidtest<-function(normData,DF=future::availableCores()*0.7){

  tpptrData<-TPP::tpptrImport(configTable = normData$TPPconfig,data=normData$TPPdata)
  
  
  check_longTable <- TPP::tpptrTidyUpESets(tpptrData)

  #TPP fit sigmoids
  Fit_sig<-TPP::tpptrCurveFit(check_longTable,doPlot=FALSE,maxAttempts=50,nCores=DF)
 
  
pval_results<-tpptrAnalyzeMeltingCurves(Fit_sig)
  return(pval_results)
}
#Rename output from PDtoMSStatsTMTformat()

Rename2_MSStatsTMT_function<-function(PSM_result,Run){
  PSM_result$BioReplicate=NA
  PSM_result$Mixture=NA
  if(any(names(PSM_result)=="Run")){
    if(any(is.na(PSM_result$Mixture))|!any(names(PSM_result)=="Mixture")){
      PSM_result$Mixture<-paste0(ifelse(stringr::str_detect(PSM_result$Run,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(PSM_result$Run,"carrier")==TRUE,"C",NA)),'_',
                                 ifelse(stringr::str_detect(PSM_result$Run,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(PSM_result$Run,"r_FAIMS")==TRUE,"F",NA)),'_',
                                 ifelse(stringr::str_detect(PSM_result$Run,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(PSM_result$Run,"S_Phi")==TRUE,"S",NA)))
    }
  }else if(any(names(PSM_result)=="Spectrum_File")){
    
    PSM_result$Mixture<-paste0(ifelse(stringr::str_detect(PSM_result$Spectrum_File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(PSM_result$Spectrum_File,"carrier")==TRUE,"C",NA)),'_',
                               ifelse(stringr::str_detect(PSM_result$Spectrum_File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(PSM_result$Spectrum_File,"r_FAIMS")==TRUE,"F",NA)),'_',
                               ifelse(stringr::str_detect(PSM_result$Spectrum_File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(PSM_result$Spectrum_File,"S_Phi")==TRUE,"S",NA)))
    
  }else{
    PSM_result$Mixture<-PSM_result$Subject
    PSM_result$Run<-Run
    PSM_result$File_ID<-PSM_result$Mixture
  }
  if(any(is.na(PSM_result$BioReplicate))&any(names(PSM_result)=="File_ID")){
    PSM_result$BioReplicate<-PSM_result$File_ID
    PSM_result$Mixture<-PSM_result$File_ID
  }
  if(any(is.na(PSM_result$TechRepMixture))&any(names(PSM_result)=="Run")){
    if(any(stringr::str_detect(PSM_result$Run,"."))){
      PSM_result$TechRepMixture<-stringr::str_extract(PSM_result$Run,"[[:digit:]]+.raw")
      PSM_result$TechRepMixture<-stringr::str_extract(PSM_result$TechRepMixture,"[[:digit:]]+")
    }else{
      PSM_result$TechRepMixture<-stringr::str_extract(PSM_result$Run,"[[:digit:]]+r")
      PSM_result$TechRepMixture<-stringr::str_extract(PSM_result$TechRepMixture,"[[:digit:]]+")
    }
  }else if(any(names(PSM_result)=="Subject")){#if this is a simulation
    len<-length(unique(stringr::str_extract(PSM_result$Protein,"_[[:digit:]]+")))
    TechReps<-data.frame(Experiment=as.factor(unique(PSM_result$Subject)),TechRepMixture=rep(c(1,2),2))
    PSM_result$Experiment<-NA
    PSM_result$Experiment<-PSM_result$Mixture
    if(any(names(PSM_result)=="TechRepMixture")){
      PSM_result<-PSM_result%>% dplyr::select(-TechRepMixture)
    }
    PSM_result<-PSM_result %>% dplyr::inner_join(TechReps,by="Experiment")
    PSM_result$BioReplicate<-NA
    
    PSM_result$Mixture<-stringr::str_remove(PSM_result$Protein,"_Sim_[[:digit:]]+")
    PSM_result$BioReplicate<-paste0(PSM_result$Mixture,"_",PSM_result$Channel)
  }else{
    
    TechReps<-data.frame(Experiment=as.factor(unique(x$Mixture)),TechRepMixture=rep(c(1,2),2))
    PSM_result$Experiment<-NA
    PSM_result$Experiment<-PSM_result$Mixture
    if(any(names(PSM_result)=="TechRepMixture")){
      PSM_result<-PSM_result%>% dplyr::select(-TechRepMixture)
    }
    PSM_result<-PSM_result %>% dplyr::inner_join(TechReps,by="Experiment")
    PSM_result$BioReplicate<-NA
    
    PSM_result$Mixture<-stringr::str_remove(PSM_result$Protein,"_Sim_[[:digit:]]+")
    PSM_result$BioReplicate<-paste0(PSM_result$Mixture,"_",PSM_result$Channel)
  }
  
  
  if(any(is.na(PSM_result$Condition))&any(names(PSM_result)=="Run")){
    PSM_result$Condition<-ifelse(stringr::str_detect(PSM_result$Run,"DMSO"),paste0(PSM_result$Channel,"_","vehicle"),paste0(PSM_result$Channel,"_","treated"))
    PSM_result$Condition<-ifelse(PSM_result$Channel=="126","Norm",PSM_result$Condition)
  }else{
    PSM_result$Condition<-as.factor(PSM_result$Group)
    PSM_result$sample_id<-PSM_result$File_ID
    
  }
  if(!any(names(PSM_result)=="Fraction")&any(names(PSM_result)=="File_ID")){
    if(any(stringr::str_detect(PSM_result$File_ID,"F[[:digit:]]+.[[:digit:]]+"))){
      PSM_result$Fraction<-stringr::str_remove(PSM_result$File_ID,"F[[:digit:]]+.")
      PSM_result$Fraction<-stringr::str_extract(PSM_result$Fraction,"[[:digit:]]+")
    }else{
      PSM_result$Fraction<-1
    }
  }
  x<-PSM_result %>% as.data.frame()
  
  return(x)
}

TPP_NPARC_calc<-function(Sim,method="NPARC",DF=5,df.temps,CARRIER=TRUE){  
  if(method=="NPARC"){
    #rename data to MSStatsTMT format
    start=proc.time()
    Sim_All_Data<-list(ProteinLevelData=suppressWarnings(MSstats_format_TO_TPP(Sim,df.temps,CARRIER=CARRIER)))
    end=proc.time()
    print(paste0("Renamed data to TPP format in  ",as.numeric(signif((end-start)[1],2))))
    
    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for spline DF=5
    start=proc.time()
    Sim_NPARC<-runTPP_splineFtest(Sim_All_Data$ProteinLevelData,DF=DF)
    end=proc.time()
    print(paste0("TPP results calculated in ",as.numeric(signif((end-start)[1],2))))
    Sim_NPARC$unmoderatedFp_val<-1-pf(Sim_NPARC$F_statistic, Sim_NPARC$df1,Sim_NPARC$df2) 
    
    return(Sim_NPARC)
  }else if(method=="TPP"){
    #rename data to MSStatsTMT format
    start=proc.time()
    Sim_All_Data<-list(ProteinLevelData=suppressWarnings(MSstats_format_TO_TPP(Sim,df.temps,CARRIER=TRUE)))
    end=proc.time()
    print(paste0("Renamed data to TPP format in  ",as.numeric(signif((end-start)[1],2))))
    
    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for spline DF=5
    start=proc.time()
    Sim_TPP<-runTPP_sigmoidtest(Sim_All_Data$ProteinLevelData,DF=DF)
    end=proc.time()
    print(paste0("TPP results calculated in ",as.numeric(signif((end-start)[1],2))))
    
    
    return(Sim_TPP)
  }
}


KL_Calc<-function(Sim_NPARC,Corr=NA){
    #remove simulation names and numbers from the data 
    Sim_NPARC$Protein<-stringr::str_remove(as.character(Sim_NPARC$uniqueID),"_Sim_[[:digit:]]+")

    Corr$Protein<-as.character(Corr$Protein)
    #bin the RE correlation from the data, include zeros on the bin
    Corr<-Corr %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,0.7,by=0.1),0.81),right=FALSE))
    #From NPARC results rename uniqueID to Protein and remove Simulation names and numbers from the data
    Sim_NPARC<-Sim_NPARC %>% dplyr::select(Protein,p_NPARC) 
    Sim_NPARC$Protein<-stringr::str_remove(as.character(Sim_NPARC$Protein),"_Sim_[[:digit:]]+")
    Sim_Cor<-Sim_NPARC %>% dplyr::inner_join(Corr,by="Protein")
    #join NPARC results to correlation matrix
    Sim_Cor<-Sim_Cor %>% dplyr::filter(!is.na(RE_corr))
    Sim_Cor<-Sim_Cor %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,0.7,by=0.1),0.81),right=FALSE))
    
    testcount<-Sim_Cor %>%
      dplyr::group_by(binCorr) %>%
      dplyr::summarise(nProteins=length(unique(Protein)),
                       TypeI_error=mean(p_NPARC<0.05,na.rm=TRUE)) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(prop_Prot=nProteins/sum(nProteins)) %>% 
      distinct() %>% 
      dplyr::distinct()
    
    Sim_Sig<-Sim_Cor %>% dplyr::filter(p_NPARC<0.05) %>%
      dplyr::group_by(binCorr) %>%
      dplyr::mutate(TypeI_error=mean(p_NPARC<0.05,na.rm=TRUE))%>%
      dplyr::ungroup()
    
    #identify the significant ones (False positives)
    test_sig<-Sim_Sig %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,0.7,by=0.1),0.81),right=FALSE))
    
    testcount_sig<-test_sig %>%
      dplyr::group_by(binCorr) %>%
      dplyr::summarise(nProteins=length(unique(Protein)))%>% 
      dplyr::ungroup() %>%
      dplyr::mutate(prop_Prot=nProteins/sum(nProteins)) %>% 
      distinct() %>% 
      dplyr::distinct()
    
    kL_bin<-testcount_sig$prop_Prot*log(testcount_sig$prop_Prot/testcount$prop_Prot)
    kL<-sum(kL_bin)
    xout<-data.frame(KL=kL)
    
    xout<-list(xout,testcount,testcount_sig)
    names(xout)<-c("KL_result","bins_all","bins_sig")
    return(xout)
    
}
  
KL_Calc_unmod<-function(Sim_NPARC,Corr=NA){
  #remove simulation names and numbers from the data 
  Sim_NPARC$Protein<-stringr::str_remove(as.character(Sim_NPARC$uniqueID),"_Sim_[[:digit:]]+")
  
  Corr$Protein<-as.character(Corr$Protein)
  #bin the RE correlation from the data, include zeros on the bin
  Corr<-Corr %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,0.7,by=0.1),0.81),right=FALSE))
  #From NPARC results rename uniqueID to Protein and remove Simulation names and numbers from the data
  Sim_NPARC<-Sim_NPARC %>% dplyr::select(Protein,unmoderatedFp_val) 
  Sim_NPARC$Protein<-stringr::str_remove(as.character(Sim_NPARC$Protein),"_Sim_[[:digit:]]+")
  Sim_Cor<-Sim_NPARC %>% separate(Protein,into=c("x","RE_corr"),sep="_") %>% dplyr::select(-x)
  #join NPARC results to correlation matrix
  Sim_Cor<-Sim_Cor %>% dplyr::filter(!is.na(RE_corr)) %>% dplyr::mutate(RE_corr=as.numeric(RE_corr))
  Sim_Cor<-Sim_Cor %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,max(Sim_Cor$RE_corr),by=0.1),max(Sim_Cor$RE_corr)+0.1),right=FALSE))
  
  testcount<-Sim_Cor %>%
    dplyr::group_by(binCorr) %>%
    dplyr::summarise(nProteins=length(unique(Protein)),
                     TypeI_error=mean(unmoderatedFp_val<0.05,na.rm=TRUE)) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(prop_Prot=nProteins/sum(nProteins)) %>% 
    distinct() %>% 
    dplyr::distinct()
  
  Sim_Sig<-Sim_Cor %>% dplyr::filter(unmoderatedFp_val<0.05) %>%
    dplyr::group_by(binCorr) %>%
    dplyr::mutate(TypeI_error=mean(unmoderatedFp_val<0.05,na.rm=TRUE))%>%
    dplyr::ungroup()
  
  #identify the significant ones (False positives)
  test_sig<-Sim_Sig %>% dplyr::mutate(binCorr=cut(RE_corr,breaks=c(seq(0,0.7,by=0.1),0.81),right=FALSE))
  
  testcount_sig<-test_sig %>%
    dplyr::group_by(binCorr) %>%
    dplyr::summarise(nProteins=length(unique(Protein)))%>% 
    dplyr::ungroup() %>%
    dplyr::mutate(prop_Prot=nProteins/sum(nProteins)) %>% 
    distinct() %>% 
    dplyr::distinct()
  
  kL_bin<-testcount_sig$prop_Prot*log(testcount_sig$prop_Prot/testcount$prop_Prot)
  kL<-sum(kL_bin)
  xout<-data.frame(KL=kL)
  
  xout<-list(xout,testcount,testcount_sig)
  names(xout)<-c("KL_result","bins_all","bins_sig")
  return(xout)
  
}
getTypeIerr<-function(Result){
  #remove simulation names and numbers from the data 
  Sim_NPARC$Protein<-stringr::str_remove(as.character(Result$uniqueID),"_Sim_[[:digit:]]+")
  
    Sim_Cor<-Sim_NPARC %>%
      tidyr::separate(Protein,sep="_",into = c("x","Corr")) %>%
      dplyr::select(-x) %>% 
      dplyr::mutate(Corr=as.numeric(Corr))
    
    #visualize intra-subject correlation
    Sim_Cor<-Sim_Cor %>% dplyr::filter(!is.na(Corr))
    
    Sim_Sig<-Sim_Cor %>% dplyr::filter(p_NPARC<0.05)
    
    testcount<-Sim_Cor %>% dplyr::group_by(Corr) %>%
      dplyr::summarise(TypeI_error=mean(p_NPARC<0.05,na.rm=TRUE)) %>% 
      dplyr::ungroup() 
    
    testcount_sig<-Sim_Sig %>%
      dplyr::group_by(Corr) %>% dplyr::count() %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(prop_Prot=n/sum(n))

    
    kL_bin<-testcount_sig$prop_Prot*log(testcount_sig$prop_Prot/testcount$prop_Prot)
    
    xout<-list(testcount,testcount_sig,tibble(overall_typeIerr=mean(Sim_Cor$p_NPARC<0.05,na.rm=TRUE),
                                              KL=round(sum(kL_bin,na.rm=TRUE),3)))
    names(xout)<-c("BinCorData","SigBinCorData","TypeI_error")
    
    return(xout)
  
}

getTypeIerr_unmod<-function(Result){
  #remove simulation names and numbers from the data 
  Sim_NPARC$Protein<-stringr::str_remove(as.character(Result$uniqueID),"_Sim_[[:digit:]]+")
  
  Sim_Cor<-Sim_NPARC %>%
    tidyr::separate(Protein,sep="_",into = c("x","Corr")) %>%
    dplyr::select(-x) %>% 
    dplyr::mutate(Corr=as.numeric(Corr))
  
  #visualize intra-subject correlation
  Sim_Cor<-Sim_Cor %>% dplyr::filter(!is.na(Corr))
  
  Sim_Sig<-Sim_Cor %>% dplyr::filter(unmoderatedFp_val<0.05)
  
  testcount<-Sim_Cor %>% dplyr::group_by(Corr) %>%
    dplyr::summarise(TypeI_error=mean(unmoderatedFp_val<0.05,na.rm=TRUE)) %>% 
    dplyr::ungroup() 
  
  testcount_sig<-Sim_Sig %>%
    dplyr::group_by(Corr) %>% dplyr::count() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(prop_Prot=n/sum(n))
  
  
  kL_bin<-testcount_sig$prop_Prot*log(testcount_sig$prop_Prot/testcount$prop_Prot)
  
  xout<-list(testcount,testcount_sig,tibble(overall_typeIerr=mean(Sim_Cor$p_NPARC<0.05,na.rm=TRUE),
                                            KL=round(sum(kL_bin,na.rm=TRUE),3)))
  names(xout)<-c("BinCorData","SigBinCorData","TypeI_error")
  
  return(xout)
  
}
