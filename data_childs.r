library(multidplyr)
library(dplyr)
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
library(MSstats) #Load MSstats package
library(readr)
library(tidyr) #Required package for normalization
library(dplyr) #Required package for normalization
library(ggplot2)
library(minpack.lm)

#Read Childs dataset
setwd("~/test_01")
tppData <- readRDS("tppData.Rds")
tppData %>% head %>% kable(digits = 2)
#generate a summary of the table contents
tppData %>%
  mutate(compoundConcentration = factor(compoundConcentration),
         replicate = factor(replicate),
         dataset = factor(dataset)) %>%
  summary()
#remove decoy proteins which contain ### prefix
tppData <- tppData %>% filter(!grepl(":###[[:alnum:]]*###",uniqueID))
#remove proteins that do not have at least one unique peptide
tppData <- filter(tppData, uniquePeptideMatches >=1)
#Remove proteins that only contain missing values
tppData <- tppData %>% filter(!is.na(relAbundance))
#remove all proteins which were not reproducibly observed with full melt
#curves in both treated groups.  #ten temperatures were recorded, so there
#should be 10 measurements
tppData <- tppData %>%
  dplyr::group_by(dataset, uniqueID) %>%
  dplyr::mutate(n=n()) %>%
  group_by(dataset) %>%
  dplyr::mutate(max_n = max(n)) %>%
  dplyr::filter(n == max_n) %>%
  dplyr::select(-n,-max_n) %>%
  ungroup
#count the remaining proteins in the dataset
tppData %>%
  dplyr::distinct(dataset,uniqueID) %>%
  dplyr::group_by(dataset) %>%
  dplyr::tally() %>%
  knitr::kable()
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
tppData <- readRDS("~/CS7290-/external data/tppData_preprocessed.Rds")
#########################################CAUTION: I assigned this to the same variable
#To look up a compound, please edit dataset (drug) and uniqueID
stk4 <- filter(tppData, dataset == "Staurosporine", uniqueID=="STK4_IPI00011488")
stk4 <- filter(tppData, dataset == "Panobinostat", uniqueID=="HDAC10")
stk4 <- filter(tppData, dataset == "Panobinostat", uniqueID=="RAP2B")


stk4 <- filter(tppData,uniqueID=="HDAC6_IPI00005711",dataset == "Dasatinib")
stk4 <- filter(tppData, dataset == "ATP", uniqueID=="C6ORF106_IPI00013269")
stk4 <- subset(tppData, uniqueID=="P36507")
stk4 <- subset(tppData, uniqueID=="P46940")
stk4<-subset(tppData, uniqueID=="Q9Y6I9")
stk4<-subset(tppData, uniqueID=="Q9Y6I9")
stk4<-subset(tppData, uniqueID=="P61457")
stk4<-subset(tppData, uniqueID=="Q9HBU6")
###Pick one stk4 from above, ignore the rest

##select data
#select desired protein

stk4 %>% filter(compoundConcentration == 20, replicate ==1) %>%
  dplyr::select(-dataset) %>% kable(digits=2)

#Generate a plot of the measurements

stk4_plot <- ggplot(stk4,aes(x=temperature, y = relAbundance))+
  geom_point(aes(shape = factor(replicate),color = factor(compoundConcentration)),
             size = 2)+theme_bw()+ggtitle(stk4$uniqueID) 
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
#Execute function

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
  group_by(compoundConcentration) %>%
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
            by = c("relAbundance","temperature","compoundConcentration")) %>%
  distinct()
#add the curves to the plot
stk4_plot <- stk4_plot + geom_line(data = distinct(stk4,temperature, compoundConcentration,alternativePredictions),
                                   aes(y=alternativePredictions,color = factor(compoundConcentration)))+ggtitle(stk4$uniqueID)
print(stk4_plot)


#compute RSS values to quantify GOF relative to null

rssPerModel <- stk4 %>%
  summarise(rssNull = sum(nullResiduals^2),
            rssAlternative = sum(alternativeResiduals^2)) %>% data.frame()
kable(rssPerModel, digits = 4)

saveRDS(tppData,"tppData_preprocessed.Rds")
tppData_preprocessed <- readRDS("~/test_01/tppData_preprocessed.Rds")
tppData<-tppData_preprocessed
#filter the data for the best compound


###############################################################Amanda's script
#Curve fitting
#splituniqueID and dataset while removing data accessions with one replicate and with no values
TppData <-tppData %>% split(list(tppData$uniqueID,tppData$dataset,tppData$compoundConcentration)) %>%  discard(function(x) nrow(x) == 0) %>% discard(function(x) nrow(x) <= 10)
d<-list.count(TppData)
#preallocate new_a

new_a<-data.frame(uniqueID = character(d), RSS= numeric(d),AUC= numeric(d),INFL= numeric(d),gof= numeric(d))
#iterate over proteins to get outputs
options(warn=2) #set warnings as errors

pb <- txtProgressBar(min = 1, max = d, style = 3)
for (i in 1:d) {
  setTxtProgressBar(pb, value = i)
  x <- try(nplr(x=TppData[[i]]$temperature, y =TppData[[i]]$relAbundance, npars=5,silent = TRUE))
  df<-length(TppData[[i]]$temperature)
  if ('nplr' %in% class(x)) {
    new_a$compoundConcentration[i]<-paste(TppData[[i]]$compoundConcentration[1])
    new_a$replicate[i]<-paste(TppData[[i]]$replicate[1])
    new_a$dataset[i]<-paste(TppData[[i]]$dataset[1])
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
saveRDS(new_a,'Chalt_pp.Rds')

new_A<-new_a
new_a<-new_A
#split the uniqueID from the dataset
test<-tidyr::separate(new_a$uniqueID,ls.TppData..1.d.,'.')
colnames(test)<-"uniqueID"
new_a<-data.frame(test,new_A$dataset,new_a$compoundConcentration,new_a$replicate,new_a$RSS,new_a$AUC,new_a$INFL,new_a$gof)
colnames(new_a)<-c("uniqueID","dataset","compoundConcentration","replicate","RSS","AUC","INFL","gof")


new_aRSS<- ddply(new_a,c("uniqueID","dataset","compoundConcentration"),summarize,RSS =sum(RSS,na.rm=TRUE)^2,AUC = mean(AUC,na.rm=TRUE),INFL = mean(INFL,na.rm=TRUE),gof=mean(gof,na.rm=TRUE))
new_a$AUCdiff <- numcolwise(.fun = function(x) {x - x[1]})(new_a[5])
new_a$INFLdiff <- numcolwise(.fun = function(x) {x - x[1]})(new_a[6])


new_a<-new_aRSS

saveRDS(new_a,'CHalt.Rds')


TppData <-tppData %>% split(list(tppData$uniqueID,tppData$dataset)) %>%  discard(function(x) nrow(x) == 0)
d<-list.count(TppData)
new<-data.frame(uniqueID = character(d), RSS= numeric(d),AUC= numeric(d),INFL= numeric(d),gof= numeric(d))
df<-length(TppData[[1]]$temperature)

head(new)
options(warn=2)
#iterate over proteins to get outputs
pb <- txtProgressBar(min = 1, max = d, style = 3)
for (i in 1:d) {
  setTxtProgressBar(pb, value = i)

  x <-try(nplr(x=TppData[[i]]$temperature, y =TppData[[i]]$relAbundance, npars=5))

  df<-length(TppData[[i]]$temperature)
  if ('nplr' %in% class(x)) {
    new$compoundConcentration[i]<-paste(TppData[[i]]$compoundConcentration[1])
    new$replicate[i]<-paste(TppData[[i]]$replicate[1])
    new$dataset[i]<-paste(TppData[[i]]$dataset[1])
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
saveRDS(new,'Chnull_pp.Rds')

New<-new
new<-New
test<-separate(new$uniqueID,ls.TppData..1.d.,'.')
colnames(test)<-"uniqueID"
new<-data.frame(test,new$compoundConcentration,new$replicate,new$dataset,new$RSS,new$AUC,new$INFL,new$gof)
colnames(new)<-c("uniqueID","compoundConcentration","replicate","dataset","RSS","AUC","INFL","gof")
new_RSS<- ddply(new,c("uniqueID","dataset"),summarize,RSS =sum(RSS)^2,AUC=diff(AUC),INFL=mean(INFL),gof=mean(gof))
new<-new_RSS %>% drop_na()

saveRDS(new,'CHNULL.Rds')

new<-new %>% subset(new$uniqueID %in% new_a$uniqueID)
new_a<-new_a %>% subset(new_a$uniqueID %in% new$uniqueID)

#join filtered data to test
r<-left_join(new,new_a,"uniqueID","dataset")
r<-r
d<-nrow(r)
#create table with results
new_RSS<- ddply(R,c("uniqueID","dataset"),summarize,RSS =r$RSS.x-r$RSS.y,AUC_d = r$AUC.x-r$AUC.y, Tm_d=r$INFL.x-r$INFL.y,gof_n=r$gof.x,
                gof_a=r$gof.y)

new_RSS$uniqueID <- r$uniqueID
new_RSS$dataset<-r$dataset.x
R<-new_RSS
colnames(R)<-c("uniqueID",'dataset',"RSSdiff","AUCdiff","TmDIFF","gof_n","gof_a")


#Filter results
test<-R %>% filter(gof_n > 0.9991) %>% filter(gof_a>0.9991) %>% filter(RSSdiff<1)%>% filter(RSSdiff> -1)
testID<-tppData %>% filter(uniqueID %in% test$uniqueID)
data<-testID%>%group_by(dataset, uniqueID) %>%
  tally() %>% subset(n>30)
#subset proteins from optimization model
test1<-test %>% subset(test$uniqueID %in% data$uniqueID)
saveRDS(test1,'filteredProteins.rds')
# t <-tppData$uniqueID[test1$uniqueID %in% test1$uniqueID]
# #Change the names according to MSStats requirement:
#

df_out<-data.frame(INTENSITY = tppData$relAbundance, CONCENTRATION = tppData$temperature,NAME = tppData$uniqueID,REPLICATE = tppData$replicate,dataset=tppData$dataset,CC=tppData$compoundConcentration)
spike <- df_out %>% filter(NAME == "STK4") %>% filter(CC== 0)
blankdata<-df_out %>% filter(NAME=="NQO2",CONCENTRATION>0)
blankdata$NAME<-"STK4"
blankdata$CONCENTRATION <-0#set all blanks to concentration zero to prevent errors
blankdata$INTENSITY<-blankdata$INTENSITY-1#Subtract one from concentration values
colnames(blankdata)[1]<-"INTENSITY"
colnames(blankdata)[2]<-"CONCENTRATION"
datain<-rbind(spike,blankdata)
#
# #FINALLY test
test<-nonlinear_quantlim(datain,alpha = 0.05,Npoints = 100,Nbootstrap=2000)
plot_quantlim(spikeindata = tes,quantlim_out = test, dir_output = "~/test_01")





DataC<-tppData %>% subset(uniqueID=="15 KDA PROTEIN._IPI00879051")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="MED21_IPI00013677")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="SUPT3H_IPI00025792")%>% subset(dataset=="ATP")
DataC<-tppData %>% subset(uniqueID=="PRKAB2_IPI00013905")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="POLR3H_IPI00289667")%>% subset(dataset=="ATP")
DataC<-tppData %>% subset(uniqueID=="HDAC10_IPI00012439")%>% subset(dataset=="Dasatinib")
DataC<-tppData %>% subset(uniqueID=="HDAC10_IPI00012439")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="HDAC10")%>% subset(dataset=="Panobinostat")#filtered out
DataC<-tppData %>% subset(uniqueID=="HDAC1_IPI00013774")%>% subset(dataset=="Dasatinib")


DataC<-tppData %>% subset(uniqueID=="RNF181")%>% subset(dataset=="Panobinostat")
DataC<-tppData %>% subset(uniqueID=="HNRNPM_IPI00171903")%>% subset(dataset=="Dasatinib")
DataC<-tppData %>% subset(uniqueID=="CDK4_IPI00007811")%>%  subset(dataset=="ATP")
DataC<-tppData %>% subset(uniqueID=="NUCB2_IPI00009123") %>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="SUZ12_IPI00299526")%>% subset(dataset=="Staurosporine")

DataC<-tppData %>% subset(uniqueID=="SF3B6")%>% subset(dataset=="Panobinostat")#
DataC<-tppData %>% subset(uniqueID=="YEATS4_IPI00008536")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="TPD52_IPI00873344")%>% subset(dataset=="Staurosporine")
DataC<-tppData %>% subset(uniqueID=="DNM1L_IPI00235412")%>% subset(dataset=="Dasatinib")
DataC<-tppData %>% subset(uniqueID=="RPS20")%>% subset(dataset=="Panobinostat")
DataC<-tppData %>% subset(uniqueID=="RPS20_IPI00794659")%>% subset(dataset=="ATP")
DataC<-tppData %>% subset(uniqueID=="SMG8")%>% subset(dataset=="Panobinostat")
DataC<-tppData %>% subset(uniqueID=="ESYT1")
DataC<-tppData %>% subset(uniqueID=="CKAP5")
DataC<-tppData %>% subset(uniqueID=="C2CD5")
DataC<-tppData %>% subset(uniqueID=="USB1")
DataC<-tppData %>% subset(uniqueID=="POLR3GL")%>% subset(dataset=="Panobinostat")
DataC<-tppData %>% subset(uniqueID=="POLR3GL_IPI00031117")%>% subset(dataset=="Staurosporine")

#  case_when(nrow(DataC)==0,


nrow(DataC)

Data0<-tppData %>% subset(tppData$uniqueID %in% DataC$uniqueID) %>% filter(.$dataset==DataC$dataset[1]) %>% filter(compoundConcentration == 0 )


Data1<-tppData %>% subset(tppData$uniqueID %in% DataC$uniqueID) %>% filter(.$dataset==DataC$dataset[1]) %>% filter(compoundConcentration >0 )


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
#plot(data0$x0,data0$y0,main='vehicle fit')
#plot(data20$x20,data20$y20,main='treated fit')

#plot the original values
orig0<-Data0
#plot(orig0$temperature,orig0$relAbundance,main=Data0$uniqueID[1])

origt<-Data1
#plot(origt$temperature,origt$relAbundance,main=Data1$uniqueID[1])
data<-rbind(orig0,origt)
#Null hypothesis

colnames(datan)<-c('temperature','relAbundance')
colnames(data0)<-c('temperature','relAbundance')
colnames(data20)<-c('temperature','relAbundance')

#For vehicle##########
#fit has 200 values, I only need 20 x values

control<-ggplot(data,aes(x=temperature, y = relAbundance))+geom_line(data = datan,color='Black')
control<-control+geom_point(aes(shape = factor(replicate),color = factor(compoundConcentration)))
control<-control+geom_line(data = data20,color='Orange')+geom_line(data = data0,color='Gray')+ggtitle(Data0$uniqueID)
print(control)

