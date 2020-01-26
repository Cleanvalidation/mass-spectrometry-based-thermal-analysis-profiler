setwd("~/cetsa_processing-master/cetsa_01-master/R")
# # df<- read_excel("eFT_30K_all5samples_PROTEINS.xlsx")
df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4','F5'), sample_name = c('MEK_1','MEK_2', 'MEK_3','DMSO_1','DMSO_3'), stringsAsFactors = FALSE)

f<-"~/cetsa_processing-master/cetsa_01-master/R/test/proteins.xlsx"
df_raw <- read_cetsa(f)
df_clean <- clean_cetsa(df_raw, temperatures = df.temps, samples = df.samples)
df_norm <- normalize_cetsa(df_clean , df.temps$temperature)
df_normC <- df_norm %>% select(-value,-correction)
df_normC <- df_normC %>% dplyr::rename(sample,c("dataset"))
colnames(df_normC)<-c("uniqueID","dataset","C","I")
tppData<-df_normC

#Test Cliff's two proteins
Data1<-tppData %>% dplyr::filter(uniqueID == "P36507",dataset%in%c("F1","F2","F3"))
Data0<-tppData %>% dplyr::filter(uniqueID == "P36507",dataset%in%c("F4","F5"))

Data1<-tppData %>% dplyr::filter(uniqueID == "Q02750",dataset%in%c("F1","F2","F3"))
Data0<-tppData %>% dplyr::filter(uniqueID == "Q02750",dataset%in%c("F4","F5"))
