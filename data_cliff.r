library(multidplyr)
library(nplr)
library(tidyverse)
library(broom)
library(knitr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
library(readxl)
library(nls2)
library(rlist)
library(stringr)
library(plyr)
library(dplyr)
library(minpack.lm)
#readCliffdataset
# # df<- read_excel("eFT_30K_all5samples_PROTEINS.xlsx")
df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4','F5'), sample_name = c('MEK_1','MEK_2', 'MEK_3','DMSO_1','DMSO_3'), stringsAsFactors = FALSE)
f<-"~/cetsa_processing-master/cetsa_01-master/R/test/proteins.xlsx"
df_raw <- read_cetsa(f)

df_clean <- clean_cetsa(df_raw, temperatures= df.temps, samples= df.samples)
df_norm <- normalize_cetsa(df_clean, temperatures = df.temps$temperature)
df_normC <- df_norm %>% select(-value,-correction)
df_normC <- df_normC %>% separate(sample,c("dataset","replicate"))
colnames(df_normC)<-c("uniqueID","dataset","replicate","temperature","relAbundance")
tppData<-df_normC

# #Read Childs dataset
# setwd("~/test_01")
# tppData <- readRDS("tppData.Rds")
# tppData %>% head %>% kable(digits = 2)
#generate a summary of the table contents
# tppData %>%
#   mutate(compoundConcentration = factor(compoundConcentration),
#          replicate = factor(replicate),
#          dataset = factor(dataset)) %>%
#   summary()
#remove decoy proteins which contain ### prefix
# tppData <- tppData %>% filter(!grepl(":###[[:alnum:]]*###",uniqueID))
# #remove proteins that do not have at least one unique peptide
# tppData <- filter(tppData, uniquePeptideMatches >=1)
# #Remove proteins that only contain missing values
# tppData <- tppData %>% filter(!is.na(relAbundance))
#remove all proteins which were not reproducibly observed with full melt
#curves in both treated groups.  #ten temperatures were recorded, so there
#should be 10 measurements
tppData <- tppData %>%
  dplyr::group_by(dataset, uniqueID) %>%
  mutate(n=n()) %>%
  dplyr::group_by(dataset) %>%
  mutate(max_n = max(n)) %>%
  filter(n == max_n) %>%
  dplyr::select(-n,-max_n) %>%
  ungroup
#count the remaining proteins in the dataset
tppData %>%
  distinct(dataset,uniqueID) %>%
  distinct %>%
  dplyr::group_by(dataset) %>%
  tally %>%
  kable()
#two concentrations of dasatinib : 0.5 uM and 5 uM
#combine both datasets into one with 4 vehicle and two treated repl
tppData <- tppData %>%
  mutate(replicate = ifelse(dataset == "Dasatinib 5",
                            yes = replicate + 2,
                            no = replicate)) %>%
  mutate(dataset = gsub("0.5| 5","",dataset))
#check result
tppData %>%
  distinct(dataset,replicate, compoundConcentration) %>%
  filter(compoundConcentration>0) %>%
  dplyr::rename('drugconcentration (treatment groups)'=compoundConcentration) %>%
  kable()
#save preprocessed
saveRDS(tppData,"tppData_preprocessed.Rds")
tppData_preprocessed <- readRDS("~/test_01/tppData_preprocessed.Rds")
tppData<-tppData_preprocessed
#filter the data for the best compound
#########################################CAUTION: I assigned this to the same variable
#To look up a compound, please edit dataset (drug) and uniqueID
stk4 <- filter(tppData, dataset == "Staurosporine", uniqueID=="STK4_IPI00011488")
stk4 <- filter(tppData, dataset == "Panobinostat", uniqueID=="HDAC10")
stk4 <- filter(tppData, dataset == "Panobinostat", uniqueID=="HDAC8")
stk4 <- filter(tppData, dataset == "ATP", uniqueID=="ABL1_IPI00221171")
stk4 <- subset(tppData, uniqueID=="P36507")
stk4 <- subset(tppData, uniqueID=="P46940")
stk4<-subset(tppData, uniqueID=="Q9Y6I9")
stk4<-subset(tppData, uniqueID=="Q9Y6I9")
stk4<-subset(tppData, uniqueID=="P61457")
stk4<-subset(tppData, uniqueID=="Q9HBU6")
###Pick one stk4 from above, ignore the rest


stk4 %>% filter(compoundConcentration == 20, replicate ==1) %>%
  dplyr::select(-dataset) %>% kable(digits=2)

#Generate a plot of the measurements

stk4_plot <- ggplot(stk4,aes(x=temperature, y = relAbundance))+
  geom_point(aes(shape = factor(replicate),color = factor(dataset)),
             size = 2)+theme_bw()+ggtitle(stk4$uniqueID) + scale_color_manual("molar drug concentration",
                                                                       values = c("#808080","#da7f2d"))
print(stk4_plot)


###Model fitting:create a function where null hypothesis fits a sigmoid curve over all data

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

##function created
#Execute null functions
stk4
nullFit <- fitSingleSigmoid(x = stk4$temperature, y = stk4$relAbundance)
summary(nullFit)

#append predictions and residuals in the same table

nullPredictions <- broom::augment(nullFit)

#view
nullPredictions %>% filter(x %in% c(46,49)) %>% kable()

#append and show fit in the plot

stk4$nullPrediction <- nullPredictions$.fitted
stk4$nullResiduals <- nullPredictions$.resid

stk4_plot <- stk4_plot + geom_line(data = stk4,aes(y=nullPrediction))

#view plot
print(stk4_plot)

#alternative hypothesis: separate melting curves per exp group

alternativePredictions <- stk4 %>%
  dplyr::group_by(dataset) %>%
  do({
    fit = fitSingleSigmoid(x = .$temperature,
                           y= .$relAbundance)
    broom::augment(fit)
  }) %>%
  ungroup %>%
  #rename columns to merge into a data frame
  dplyr::rename(alternativePredictions = .fitted,
                alternativeResiduals = .resid,
                temperature = x,
                relAbundance = y)
stk4 <- stk4 %>%
  left_join(alternativePredictions,
            by = c("relAbundance","temperature","dataset")) %>%
  distinct()
#add the curves to the plot
stk4_plot <- stk4_plot + geom_line(data = distinct(stk4,temperature,alternativePredictions,dataset),
                                   aes(y=alternativePredictions,color = factor(dataset)))
print(stk4_plot)


#compute RSS values to quantify GOF relative to null

rssPerModel <- stk4 %>%
  summarise(rssNull = sum(nullResiduals^2),
            rssAlternative = sum(alternativeResiduals^2))
kable(rssPerModel, digits = 4)
###############################################################Amanda's script
#Curve fitting
#splituniqueID and dataset
TppData <-tppData %>% split(list(tppData$uniqueID,tppData$dataset))
d<-list.count(TppData)
#preallocate new_a

new_a<-data.frame(uniqueID = character(d), RSS= numeric(d),AUC= numeric(d),INFL= numeric(d),gof= numeric(d))
#iterate over proteins to get outputs
pb <- txtProgressBar(min = 1, max = d, style = 3)
for (i in 1:d) {
  options(warn=2)
  setTxtProgressBar(pb, value = i)
  x <- try(nplr(x=TppData[[i]]$temperature, y =TppData[[i]]$relAbundance, npars=5,silent =TRUE))
  df<-length(TppData[[i]]$temperature)
  if ('nplr' %in% class(x)) {
    new_a$RSS[i]<-sqrt(as.numeric(getStdErr(x)[2])*(df-2))
    new_a$AUC[i]<-as.numeric(getAUC(x)[2])
    new_a$INFL[i]<-10^as.numeric(getInflexion(x)[1])
    new_a$gof[i]<-as.numeric(getGoodness(x)[2])
  } else {
    new_a[i, c('RSS', 'AUC', 'INFL', 'gof')] <- NA
  }

}
close(pb)
new_a$uniqueID <- data.frame(ls(TppData)[1:d])
saveRDS(new_a,'Calt_pp.Rds')

new_A<-new_a
#split the uniqueID from the dataset
test<-separate(new_a$uniqueID,ls.TppData..1.d.,'.')
colnames(test)<-"uniqueID"
new_a<-data.frame(test,new_A$dataset,new_a$compoundConcentration,new_a$replicate,new_a$RSS,new_a$AUC,new_a$INFL,new_a$gof)
colnames(new_a)<-c("uniqueID","dataset","compoundConcentration","replicate","RSS","AUC","INFL","gof")

new_aRSS<- ddply(new_a,c("uniqueID","dataset","compoundConcentration"),summarize,RSS =sum(RSS)^2,AUC=mean(AUC),INFL=mean(INFL),gof=mean(gof))
new_a<-new_aRSS %>% drop_na()


saveRDS(new_a,'Calt.Rds')


TppData <-tppData %>% split(list(tppData$uniqueID))
d<-list.count(TppData)
new<-data.frame(uniqueID = character(d), RSS= numeric(d),AUC= numeric(d),INFL= numeric(d),gof= numeric(d))
df<-length(TppData[[1]]$temperature)
head(new)
#iterate over proteins to get outputs
pb <- txtProgressBar(min = 1, max = d, style = 3)
for (i in 1:d) {
  setTxtProgressBar(pb, value = i)
  options(warn=2)
  x <- try(nplr(x=TppData[[i]]$temperature, y =TppData[[i]]$relAbundance, npars=5))
  df<-length(TppData[[i]]$temperature)
  if ('nplr' %in% class(x)) {
    new$RSS[i]<-as.numeric(getStdErr(x)[2]*(df-2))
    new$AUC[i]<-as.numeric(getAUC(x)[2])
    new$INFL[i]<-10^as.numeric(getInflexion(x)[1])
    new$gof[i]<-as.numeric(getGoodness(x)[2])
  } else {
    new[i, c('RSS', 'AUC', 'INFL', 'gof')] <- NA
  }

}
close(pb)
new$uniqueID <- data.frame(ls(TppData)[1:d])
test<-separate(new$uniqueID,ls.TppData..1.d.,'.')
colnames(test)<-"uniqueID"
new<-data.frame(test,new$compoundConcentration,new$replicate,new$dataset,new$RSS,new$AUC,new$INFL,new$gof)
colnames(new)<-c("uniqueID","compoundConcentration","replicate","dataset","RSS","AUC","INFL","gof")
new_RSS<- ddply(new,c("uniqueID","dataset"),summarize,RSS =sum(RSS)^2,AUC=mean(AUC),INFL=mean(INFL),gof=mean(gof))
new<-new_RSS %>% drop_na()

head(new)
saveRDS(new,'CNULL.Rds')

# new_a<-new_a%>% subset(.$RSS<2) %>% subset(.$RSS>-2)
#
# new<-new %>% subset(.$RSS<2) %>% subset(.$RSS>-2)

#create table with results
r<-data.frame(new$uniqueID,new$RSS-(new_a$RSS),new$AUC-new_a$AUC,new$INFL-new_a$INFL)
colnames(r)<-c("uniqueID","RSSdiff","AUCdiff","TmDIFF")
##########################################
#save output
d <-length(TppData)
d <- data.frame(uniqueID = character(d), RSSdiff= numeric(d),AUCdiff= numeric(d),INFLdiff= numeric(d))

d$RSSdiff <- d$RSS.x-d$RSS.y
d$AUCdiff<-d$AUC.x-d$AUC.y
d$INFLdiff<-d$INFL.x-d$INFL.y
results<-data.frame(d$uniqueID,d$dataset,d$RSSdiff,d$AUCdiff,d$INFLdiff)

hist(results$RSSdiff)
rss5PL<-kable(results,digits = 4)


DataC<-tppData %>% subset(uniqueID=="P36507")
DataC<-tppData %>% subset(uniqueID=="P46940")
DataC<-tppData %>% subset(uniqueID=="O15037")

Data0<-tppData %>% subset(tppData$uniqueID %in% DataC$uniqueID) %>% filter(.$dataset==unique(DataC$dataset)[1])
#if this is 2d CETSA
#%>% filter(compoundConcentration == 0 )

Data1<-tppData %>% subset(tppData$uniqueID %in% DataC$uniqueID) %>% filter(.$dataset==unique(DataC$dataset)[2])

#if this is 2d CETSA
#%>% filter(compoundConcentration > 0 )
#Null hypothesis
x<-nplr(x=DataC$temperature,y=DataC$relAbundance,silent=FALSE,npars=5)
xn <-as.numeric(10^getXcurve(x))
yn <-as.numeric(getYcurve(x))
t <- unique(DataC$temperature)
#get fitted temperature values
veh <-data.frame(xn,yn)
datan<-veh %>% subset(t %in% xn)


#alt Hypothesis
x<-nplr(x=Data0$temperature,y=Data0$relAbundance,silent=FALSE,npars=5)
x0<-as.numeric(10^getXcurve(x))
y0 <-as.numeric(getYcurve(x))
veh <-data.frame(x0,y0)
data0<-veh %>% subset(t %in% x0)

x<-nplr(x=Data1$temperature,y=Data1$relAbundance,silent=FALSE,npars=5)
x20 <-as.numeric(10^getXcurve(x))
y20 <-as.numeric(getYcurve(x))
veh <-data.frame(x20,y20)
data20<-veh %>% subset(t %in% x20)
#plot the fitted values
plot(data0$x0,data0$y0,main='vehicle fit')
plot(data20$x20,data20$y20,main='treated fit')

#plot the original values
orig0<-Data0
plot(orig0$temperature,orig0$relAbundance,main=Data0$uniqueID[1])

origt<-Data1
plot(origt$temperature,origt$relAbundance,main=Data1$uniqueID[1])
data<-rbind(orig0,origt)
#Null hypothesis

colnames(datan)<-c('temperature','relAbundance')
colnames(data0)<-c('temperature','relAbundance')
colnames(data20)<-c('temperature','relAbundance')

#For vehicle##########
#fit has 200 values, I only need 20 x values

control<-ggplot(data,aes(x=temperature, y = relAbundance))+geom_line(data = datan,color='Black')
print(control)
control<-control+geom_point(aes(shape = factor(replicate),color = factor(dataset)))
print(control)
control<-control+geom_line(data = data20,color='Orange')+geom_line(data = data0,color='Gray')+ggtitle(Data0$uniqueID)
print(control)
