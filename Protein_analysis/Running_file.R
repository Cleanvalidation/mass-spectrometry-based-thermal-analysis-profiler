source('Libraries.R')
source('Rcpp_quantile.R')
source('Datapreprocessing.R')
source('Normalize_cetsa.R')
source('statistics_cetsa.R')
memory.limit(175921900000)#set for 16 GB RAM
plan(multicore,workers=availableCores())
options(future.globals.maxSize = 8000 * 1024^2)

df.t <- function(n,temperatures,protein_path,sample_mapping_name=NA){
  if(!is.logical(sample_mapping_name)){
    TMT<-read_xlsx(sample_mapping_name) %>% 
      dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample"="Sample","sample_name"="MS_sample_number","time_point"="Time_point")  
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


df.temps<-df.t(11,temperatures=c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67,68))
df_raw <- read_cetsa("C:/Users/tushi/OneDrive/Desktop/Project","C:/Users/tushi/OneDrive/Desktop/Project","_Proteins",Peptide="PG",Frac=TRUE,CFS=FALSE,solvent="Control",CARRIER=FALSE,rank=TRUE,sub=50,temperatures=df.temps,baseline="min",NORM="QUANTILE")

