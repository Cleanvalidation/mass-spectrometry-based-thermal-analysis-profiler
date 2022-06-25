library(scam)
library(mgcv)
library(marginaleffects)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(ggplot2)

library(testthat)

library(doParallel)
library(foreach)


filter_good_data = function(data){
  
  good_data = data |> 
    group_by(Accession,sample_name) |> filter(sample_name=="C_F_E",
                                              length(unique(treatment))==2,
                                              length(I)>= 30,
                                              length(I[treatment=="vehicle"])>=20,
                                              length(I[treatment=="treated"])>=20) |>dplyr::group_split()
  good_data<-good_data|>purrr::keep(function(y) sum(is.na(y$relAbundance))/nrow(y)<0.3)|> dplyr::bind_rows()
  #max(I) <=1.5) |> ungroup()
  return(good_data)
}

humandata = readRDS("All_Proteins_Human_dataset.RDS")
humandata_good = humandata|>dplyr::bind_rows() |> filter_good_data() 
# ct = humandata_good |> group_by(Accession) |> summarise(n=n(),ntmt=length(unique(treatment)))
# min(ct$n) #check filters
# min(ct$ntmt)


ground_truth_data = readRDS("Ground_truth_Human_dataset.rds")
ground_truth_data_good = ground_truth_data |> filter_good_data()



bothdata_good<-dplyr::bind_rows(humandata_good,ground_truth_data_good)

bothdata<-dplyr::bind_rows(humandata,ground_truth_data)

#a = ground_truth_data |> filter(Accession=="P36507") #Bad: O60330,O43805,O43829 #Good: P3607, Q02750
# a = bothdata_good |> filter(Accession=="O60306")
# asplit = a |> group_split(temperature)
# 
# ground_truth_split = ground_truth_data_good |> group_split(Accession,temperature)
# 


fit_scam_nlm = function(accession_data){
  mod = scam(I ~ s(temperature,by=treatment,bs="mpd",k=5) + treatment,data=accession_data,
             optimizer="nlm")
  return(mod)
}

fit_scam_default = function(accession_data){
  mod = scam(I ~ s(temperature,by=treatment,bs="mpd",k=5) + treatment,data=accession_data)
  return(mod)
}


fit_scam = function(accession_data){
  ##Uses efs optimizer -- use this function
  mod = scam(I ~ s(temperature,by=treatment,bs="mpd",k=5) + treatment,data=accession_data,
             optimizer="efs")
  return(mod)
}

poss_fit_scam = possibly(.f=fit_scam,otherwise="Error")
#get Tm values from marginal effectgs

grid_search_Tm = function(model,start=34,end=67,by=0.5){
  #initial grid search is wide, then zooms in 
  #Finds the location where derivative dY/dtemp is most negative (hence which.min)
  
  #first stage 
  initial_temps = seq(start,end,by)
  initial_grid_treatment = marginaleffects(model,datagrid(temperature=initial_temps,treatment="treated"),var="temperature")
  initial_grid_vehicle = marginaleffects(model,datagrid(temperature=initial_temps,treatment="vehicle"),var="temperature")
  
  initial_Tm_treatment = initial_grid_treatment$temperature[which.min(initial_grid_treatment$dydx)]
  initial_Tm_vehicle = initial_grid_vehicle$temperature[which.min(initial_grid_vehicle$dydx)]
  
  #zoomed in stage 
  zoomed_temps_treatment = seq(0.9*initial_Tm_treatment,1.1*initial_Tm_treatment,0.1)
  zoomed_temps_vehicle = seq(0.9*initial_Tm_vehicle,1.1*initial_Tm_vehicle,0.1)
  zoomed_grid_treatment = marginaleffects(model,datagrid(temperature=zoomed_temps_treatment,treatment="treated"),var="temperature")
  zoomed_grid_vehicle = marginaleffects(model,datagrid(temperature=zoomed_temps_vehicle,treatment="vehicle"),var="temperature")
  
  Tm_treatment = zoomed_grid_treatment$temperature[which.min(zoomed_grid_treatment$dydx)]
  Tm_vehicle = zoomed_grid_vehicle$temperature[which.min(zoomed_grid_vehicle$dydx)]
  
  result = tibble(Tm_treatment=Tm_treatment,Tm_vehicle=Tm_vehicle)
  return(result)
}
fit_scam_marginal_ATE = function(accession_data){
  stopifnot(length(unique(accession_data$Accession)) == 1)
  #model = scam(I ~ s(temperature,by=treatment,bs="mpd",k=5) + treatment,data=accession_data)
  model = poss_fit_scam(accession_data)
  if("character" %in% class(model)){
    error_df = tibble(p.value=NA,Accession=unique(accession_data$Accession))
    return(error_df)
  }
  adj_r2 = summary(model)$r.sq
  ATE = comparisons(model,variables="treatment") |> tidy() #average treatment effect
  ATE$adj_r2 = adj_r2
  ATE$Accession = unique(accession_data$Accession)
  ATE = cbind(ATE,grid_search_Tm(model)) #get Tm's 
  return(ATE)
}
fit_scam_marginal_DATE = function(accession_data){
  stopifnot(length(unique(accession_data$Accession)) == 1)
  #model = scam(I ~ s(temperature,by=treatment,bs="mpd",k=5) + treatment,data=accession_data)
  model = poss_fit_scam(accession_data)
  if("character" %in% class(model)){
    error_df = tibble(p.value=NA,Accession=unique(accession_data$Accession),adj_r2=NA,slope.pval=NA)
    return(error_df)
  }
  adj_r2 = summary(model)$r.sq
  ATE = marginaleffects(model,
                        var="temperature",
                        newdata=datagrid(treatment=c("treated","vehicle"),temperature = seq(37,67,length.out=nrow(predict(model)))),
                        hypothesis=c(rep(1/nrow(predict(model)),nrow(predict(model))),rep(-1/nrow(predict(model)),nrow(predict(model))))) |> tidy() #DATE
  ATE$adj_r2 = adj_r2
  ATE$Accession = unique(accession_data$Accession)
  ATE$slope.pval=ATE$p.value
  ATE<-ATE %>% dplyr::select(-p.value)|>dplyr::bind_rows()
  return(ATE)
}
permute_treatment = function(accession_temp_data){
  labels = accession_temp_data$treatment
  
  permuted_labels = sample(labels,size=length(labels),replace=FALSE)
  
  permuted_data = accession_temp_data |> mutate(treatment = permuted_labels)
  return(permuted_data)
}

permute_treatment_within_group = function(split_group_data){
  
  permuted_within_group_data = map(split_group_data,.f=permute_treatment) |> bind_rows()
  return(permuted_within_group_data)
}

unittest_permutation_within_group = function(accession_data){
  count_orig = accession_data |> group_by(temperature,treatment,sample_name) |> count()
  split_accession_data = accession_data |> group_split(temperature)
  permuted_data = permute_treatment_within_group(split_accession_data)
  count_permute  = permuted_data |> group_by(temperature,treatment,sample_name) |> count()
  expect_equal(count_orig,count_permute)
  print("Test Passed")
}



unittest_permutation_within_group(dplyr::bind_rows(bothdata_good))



## Computes p values and returns df-- for original p values

compute_pvalues_DATE = function(fulldata,workers=4){
  
  data_gdf_accession = fulldata |> group_split(Accession,sample_name)
  plan(multisession,workers = workers)
  original_result = future_map(.f=fit_scam_marginal_DATE,data_gdf_accession) |> bind_rows() 
  
  DATE_result = future_map2(data_gdf_accession,as.list(1:length(data_gdf_accession)),
                            function(x,y){message(y)
                              z<-fit_scam_marginal_DATE(x)
                              return(z)}) |>
    relocate(Accession,.before=type)
  DATE_result<-cbind(DATE_result,original_result)#get Tm values from original_result
  DATE_result$dTm<-DATE_result$Tm_treatment-DATE_result$Tm_vehicle
  return(DATE_result)
}

compute_pvalues_ATE = function(fulldata,workers=4){
  
  data_gdf_accession = fulldata |> group_split(Accession,sample_name)
  plan(multisession,workers = workers)
  original_result = future_map(.f=fit_scam_marginal_ATE,data_gdf_accession) |> bind_rows() 
  
  original_result$dTm<-original_result$Tm_treatment-original_result$Tm_vehicle
  return(original_result)
}
## Computes p values but returns vector -- meant for pertmutation p values bc no df accession label info needed

compute_pvalues_vec = function(fulldata,workers=8){
  data_gdf_accession = fulldata |> group_split(Accession,sample_name)
  plan(multisession,workers = workers)
  original_result = future_map(.f=fit_scam_marginal,data_gdf_accession) |> bind_rows() |>
    pull(p.value)
  DATE_result = future_map(.f=fit_scam_marginal,data_gdf_accession) |> bind_rows() |>
    pull(p.value)
  return(DATE_result)
}


## 

compute_permutation_null_dist = function(fulldata,workers=8,runs=10){
  fulldata_split = fulldata |> group_split(Accession,temperature,sample_name)
  
  p_value_list = vector(mode="list",length = runs)
  for(i in 1:runs){
    fulldata_permuted = permute_treatment_within_group(fulldata_split)
    p_value_list[[i]] = compute_pvalues_vec(fulldata_permuted,workers=workers)
  }
  p_values = unlist(p_value_list)
  return(p_values)
}



#system.time(gpermuted <- permute_treatment_within_group(ground_truth_split))

accessions = unique(bothdata_good$Accession)
# 
# partialbothdata_good = bothdata_good |> filter(Accession %in% accessions[1:300]) #just to test initially 
# 
# accessions[290:295] #this had the bad range


start=proc.time()
og_pvals = compute_pvalues_ATE(bothdata_good)
end=proc.time()
print(end-start)

#benchmark # of fitted proteins for CFE
ATE_all<-nrow(og_pvals)
ATE_pval_0.05<-nrow(og_pvals %>% dplyr::filter(p.value<0.05))
ATE_pval_R2<-nrow(og_pvals %>% dplyr::filter(p.value<0.05,adj_r2>0.8))
ATE_pval_R2_GT<-

start=proc.time()
operm = compute_permutation_null_dist(bothdata_good,runs=10) #1 run to test for now
end=proc.time()
print(end-start)

Get_FDR<-function(og_pvals,operm,alpha,B){
  pj<-og_pvals
  pj_b<-operm
  #get parameters for FDR for a two-sample t-test
  R = sum(nrow(pj[pj$p.value<alpha,]),na.rm=TRUE) #alternative
  V_est = sum(length(pj_b[pj_b<alpha]),na.rm=TRUE)/B #null hypothesis
  #compute plug-in FDR for a two-sample t-test
  FDR = V_est/R
  FDR=round(FDR,3)
  df = data.frame(FDR=FDR,alpha=alpha,B=B)
  return(df)
}
FDR_0_05<-Get_FDR(og_pvals,operm,0.05,10)
FDR_0_01<-Get_FDR(og_pvals,operm,0.01,10)

FDR_df<-rbind(FDR_0_01,FDR_0_05)


possible_shifters<-og_pvals[og_pvals$p.value<0.05,][og_pvals[og_pvals$p.value<0.05,]$Accession %in% ground_truth_data$Accession,]
#high and low confidence
# "O43252" "O95352" "P02686" "P02768" "P04406" "P06730" "P10398" "P11142" "P11233" "P20936" "P21359" "P23458"
# "P27348" "P27361" "P28482" "P31785" "P36507" "P40763" "P42336" "P42345" "P42574" "P49841" "P51114" "P60709"
# "P60763" "P61981" "P62879" "P62987" "Q00653" "Q02750" "Q13153" "Q13541" "Q7LG56" "Q9UHA4" "Q9Y2Q5" "Q9Y6R4"
#only high confidence
possible_shifters_hc<-og_pvals[og_pvals$p.value<0.05,][og_pvals[og_pvals$p.value<0.05,]$Accession %in% ground_truth_data[which(ground_truth_data$confidence=="high"&ground_truth_data$Protein.FDR.Confidence.Combined=="High"),]$Accession,] 
# "P10398" "P21359" "P23458" "P27348" "P27361" "P28482" "P31785" #"P36507" "P40763" "P42336" "P49841" "P60709"
#"P61981" #"Q02750" #"Q13153" "Q7LG56" "Q9UHA4" "Q9Y2Q5" "Q9Y6R4"
a = dplyr::bind_rows(bothdata_good) |> filter(Accession=="O75165")
moda_nlm = fit_scam_nlm(a)
moda_default = fit_scam_default(a)
moda_efs = fit_scam(a)
a$predI_nlm = predict(moda_nlm)
a$predI_default = predict(moda_default)
a$predI_efs = predict(moda_efs)

ggplot(a,aes(x=temperature,y=I,col=treatment)) + geom_point()+ geom_line(aes(x=temperature,y=predI_efs))+ggtitle(a$Accession)

#marginal effects

tidy((marginaleffects(moda_nlm)))
tidy((marginaleffects(moda_default)))
tidy((marginaleffects(moda_efs)))

plot_cap(moda_nlm,effect="temperature",condition="treatment")
plot_cap(moda_default,effect="temperature",condition="treatment")
plot_cap(moda_efs,effect="temperature",condition="treatment")+ggtitle(a$Accession)

#from the list of candidates, calculate marginal effects on these
possible_shifters_hc<-dplyr::bind_rows(possible_shifters_hc) %>%
  dplyr::group_by(Accession) %>%
  dplyr::group_split()

bothdata_good<-dplyr::bind_rows(bothdata_good) %>%
  dplyr::filter(Accession %in% dplyr::bind_rows(possible_shifters_hc)$Accession) %>% 
  dplyr::group_by(Accession) %>% 
  dplyr::group_split() 
Get_me_targets<-function(a,ground_truth){
  ME <- a %>%
    dplyr::mutate(stderr_great_than_estimate=ifelse((abs(a$conf.high)-abs(a$conf.low))<abs(a$estimate)&a$adj_r2>0.8,TRUE,FALSE))

  moda_efs <-fit_scam(ground_truth)
  ground_truth1<-ground_truth %>%
    dplyr::mutate(Tm=list(grid_search_Tm(moda_efs))) %>% 
    tidyr::unnest(Tm) %>% 
    dplyr::mutate(dTm=Tm_treatment-Tm_vehicle) %>% 
    dplyr::select(Accession,Tm_treatment,Tm_vehicle,dTm,Protein.FDR.Confidence.Combined,confidence)
  
  df<-ground_truth1 %>% dplyr::right_join(ME) %>% distinct(.)
  return(df)
  
}
hi<-purrr::map2(possible_shifters_hc,bothdata_good,function(x,y)Get_me_targets(x,y)) %>% dplyr::bind_rows()

ordered_by_estimates_CFE_<-hi[order(hi$estimate),]
Target_list_CFE_<-hi[hi$stderr_great_than_estimate==TRUE,]
Target_list_CFE_<-Target_list[order(Target_list$dTm,decreasing=TRUE),]


#get Venn Diagrams with limma
#make sure the results are present in ground truth 
ATE_DF_GT<-original_result%>% dplyr::filter(Accession %in% ground_truth_data_good$Accession)
DATE_DF_GT<-DATE_result%>% dplyr::filter(Accession %in% ground_truth_data_good$Accession)
#make sure the results are filtered by p-value based on the respective test
ATE_DF_GT_filt<-ATE_DF_GT%>% dplyr::filter(Accession %in% ground_truth_data_good$Accession)
DATE_DF_GT_filt<-DATE_DF_GT%>% dplyr::filter(Accession %in% ground_truth_data_good$Accession)

ATE_DF$intersection<-ifelse(ATE_DF$Accession%in%DATE_DF$Accession,1,0)
ME_DF<-ATE_DF[,(ncol(ATE_DF)-1):ncol(ATE_DF)]
a<-vennCounts(ME_DF)
vennDiagram(a)

