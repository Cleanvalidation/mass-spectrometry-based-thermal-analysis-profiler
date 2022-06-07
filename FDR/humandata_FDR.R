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
  good_data = data |> group_by(Accession) |> filter(length(unique(treatment))==2,
                                                    length(I)>= 30) |> ungroup()
                                                    #max(I) <=1.5) |> ungroup()
  return(good_data)
}

humandata = readRDS("All_Proteins_Human_dataset.rds")
humandata_good = humandata |> filter_good_data()
ct = humandata_good |> group_by(Accession) |> summarise(n=n(),ntmt=length(unique(treatment)))
min(ct$n) #check filters
min(ct$ntmt)


ground_truth_data = readRDS("Ground_truth_Human_dataset.rds")
ground_truth_data_good = ground_truth_data |> filter_good_data()

bothdata_good = bind_rows(humandata_good,ground_truth_data_good)
rm(humandata_good,humandata)


#a = ground_truth_data |> filter(Accession=="P36507") #Bad: O60330,O43805,O43829 #Good: P3607, Q02750
a = bothdata_good |> filter(Accession=="O60306")
asplit = a |> group_split(temperature)

ground_truth_split = ground_truth_data_good |> group_split(Accession,temperature)



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


moda_nlm = fit_scam_nlm(a)
moda_default = fit_scam_default(a)
moda_efs = fit_scam(a)
a$predI_nlm = predict(moda_nlm)
a$predI_default = predict(moda_default)
a$predI_efs = predict(moda_efs)

ggplot(a,aes(x=temperature,y=I,col=treatment)) + geom_point()+ geom_line(aes(x=temperature,y=predI_efs))
tidy((marginaleffects(moda_nlm)))
tidy((marginaleffects(moda_default)))
tidy((marginaleffects(moda_efs)))

plot_cme(moda_nlm,effect="temperature",condition="treatment")
plot_cme(moda_default,effect="temperature",condition="treatment")
plot_cme(moda_efs,effect="temperature",condition="treatment")



fit_scam_marginal = function(accession_data){
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
  count_orig = accession_data |> group_by(temperature,treatment) |> count()
  split_accession_data = accession_data |> group_split(temperature)
  permuted_data = permute_treatment_within_group(split_accession_data)
  count_permute  = permuted_data |> group_by(temperature,treatment) |> count()
  expect_equal(count_orig,count_permute)
  print("Test Passed")
}


unittest_permutation_within_group(ground_truth_data)



## Computes p values and returns df-- for original p values

compute_pvalues_df = function(fulldata,workers=4){
  data_gdf_accession = fulldata |> group_split(Accession)
  plan(multisession,workers = workers)
  original_result = future_map(.f=fit_scam_marginal,data_gdf_accession) |> bind_rows() |>
    relocate(Accession,.before=type)
  return(original_result)
}

## Computes p values but returns vector -- meant for pertmutation p values bc no df accession label info needed

compute_pvalues_vec = function(fulldata,workers=8){
  data_gdf_accession = fulldata |> group_split(Accession)
  plan(multisession,workers = workers)
  original_result = future_map(.f=fit_scam_marginal,data_gdf_accession) |> bind_rows() |>
    pull(p.value)
  return(original_result)
}


## 

compute_permutation_null_dist = function(fulldata,workers=8,runs=10){
  fulldata_split = fulldata |> group_split(Accession,temperature)
  
  p_value_list = vector(mode="list",length = runs)
  for(i in 1:runs){
    fulldata_permuted = permute_treatment_within_group(fulldata_split)
    p_value_list[[i]] = compute_pvalues_vec(fulldata_permuted,workers=workers)
  }
  p_values = unlist(p_value_list)
  return(p_values)
}



system.time(gpermuted <- permute_treatment_within_group(ground_truth_split))

accessions = unique(bothdata_good$Accession)

partialbothdata_good = bothdata_good |> filter(Accession %in% accessions[1:300]) #just to test initially 

accessions[290:295] #this had the bad range

start=proc.time()
og_pvals = compute_pvalues_df(bothdata_good)
end=proc.time()
print(end-start)

og_pvals |> filter(is.na(p.value))


start=proc.time()
operm = compute_permutation_null_dist(bothdata_good,runs=1) #1 run to test for now
end=proc.time()
print(end-start)


## TO DO: original NA p value proteins should also not be considered in the permutation, and should be filtered out first. 