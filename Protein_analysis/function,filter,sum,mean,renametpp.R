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
      data<-parallel::mclapply(data,replicate_labels,mc.cores=availableCores())
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
      data<-parallel::mclapply(data,replicate_labels,mc.cores=availableCores())
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
                            nCores = availableCores(),
                            resultPath = resultPath, 
                            plotCurves = FALSE,
                            normalize = FALSE)
  return(TRresults)
}
filter_Peptides<-function(df_,S_N,PEP,XCor,Is_Int,Missed_C,Mods,Charg,DeltaMppm,Occupancy,filter_rank=FALSE,keep_shared_proteins=FALSE,CFS=TRUE,Frac=FALSE){
  #remove shared proteins
  if(!any(names(df_)=="sample_name")){
    df_$sample_name<-df_$Spectrum.File[1]
  }
  if(any(stringr::str_detect(names(df_),"Percolator.PEP"))){
    pep<-names(df_)[stringr::str_detect(names(df_),"Percolator.PEP")]
    df_<-df_ %>% dplyr::rename("Percolator_PEP"=pep)
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
  if(any("Average_Reporter_S_N" %in% names(df_))){
    df_<-df_ %>% dplyr::filter("Average_Reporter_S_N">S_N,Percolator_PEP<PEP,Charge<Charg,MissedCleavages<Missed_C,abs(DeltaM)<DeltaMppm)
    df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","S_N"="Average_Reporter_S_N","PEP"="Percolator_PEP",
                              "IonInjTime"="Ion_Inject_Timems_",
                              "I_Interference"="Isolation_Interference_")
    if(any(names(df_)=="Channel_Occupancy_")){
      df_<-df_ %>% dplyr::filter(Channel_Occupancy_>Occupancy)
    }
    #rank by the highest intensity channel
    if(isTRUE(Frac)){
      rank<-df_%>% dplyr::filter(temp_ref=="126") %>% 
        dplyr::mutate(rank=dplyr::ntile(.$S_N,3)) %>% dplyr::select(uniqueID,Annotated_Sequence,treatment,Modifications,Spectrum.File,Annotated_Sequence,rank,S_N,PEP,MissedCleavages,sample_id,sample_name,Fraction,Positions_in_Master_Proteins)
    }else{
      rank<-df_%>% dplyr::filter(temp_ref=="126") %>% 
        dplyr::mutate(rank=dplyr::ntile(.$S_N,3)) %>% dplyr::select(uniqueID,Annotated_Sequence,treatment,Modifications,Spectrum.File,Annotated_Sequence,rank,S_N,PEP,MissedCleavages,sample_id,sample_name,replicate,Positions_in_Master_Proteins)
    }
    #remove na values in uniqueID's
    if(isTRUE(Frac)){
      rank<-rank %>%dplyr::filter(!is.na(rank),!is.na(uniqueID)) %>%  dplyr::select(uniqueID,Annotated_Sequence,treatment,Modifications,Spectrum.File,sample_id,sample_name,rank,Fraction,S_N,Positions_in_Master_Proteins) %>% distinct(.)
    }else{
      rank<-rank %>%dplyr::filter(!is.na(rank),!is.na(uniqueID)) %>%  dplyr::select(uniqueID,Annotated_Sequence,treatment,Modifications,Spectrum.File,sample_id,sample_name,rank,replicate,S_N,Positions_in_Master_Proteins) %>% distinct(.)
    }
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
    
    
    
  }else{
    
    #df_<-df_ %>% dplyr::rename("Missed_Cleavages"="#_Missed_Cleavages") 
    df_<-df_ %>% dplyr::filter(Percolator_PEP<PEP,Charge<Charg,MissedCleavages<Missed_C,abs(DeltaM)<DeltaMppm)
    df_<-df_%>% dplyr::rename("uniqueID"="Accession","I"="value","PEP"="Percolator_PEP")
    if(length(XCor)==2){
      df_<-df_ %>% dplyr::filter(XCorr>XCor[1],XCorr<XCor[2])
    }else{
      df_<-df_ %>% dplyr::mutate(XCor_l=ifelse(Charge==2 & XCorr > 1.8,TRUE,ifelse(Charge>2 & XCorr > XCor,TRUE,FALSE)))
    }
    df_<-df_ %>% distinct(.)
    #remove the carrier channel 
    rank<-df_ %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(I,3)) %>% dplyr::select(-temp_ref,-I)
    #remove na values in uniqueID's
    rank<-rank %>%dplyr::filter(!is.na(rank)) %>%  dplyr::select(Annotated_Sequence,uniqueID,sample_id,sample_name,treatment,rank) %>% distinct(.)
    #convert to data.table
    rank<-data.table::data.table(rank)
    df_<-data.table::data.table(df_)
    #get unique PSMs
    rank<-unique(rank)
    
    #Use I column to distinguish duplicated rank values/PSM
    col_n<-dplyr::intersect(names(rank),names(df_))
    #set join columns
    data.table::setkeyv(rank,cols=col_n)
    data.table::setkeyv(df_,cols=col_n)
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