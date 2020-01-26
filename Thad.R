---title: "R Notebook"
output: html_notebook
---
  ```{R}
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
library(MSstats) # Load MSstats package
library(readr)
library(tidyr) # Required package for normalization
library(dplyr) # Required package for normalization
library(ggplot2)
library(RCurl)
library(gridExtra)
library(minpack.lm)
library(pspline)
library(xts)
theme_set(theme_bw())


```
Start with the processed data
```{R Read Preprocessed Data}
tppData <- readRDS("tppData_preprocessed.Rds")
```
Visualize it
```{R Visualize}
tppData %>% head %>% kable(digits = 2)
```
Summarize it
```{R Summarize}
tppData %>%
  mutate(
    compoundConcentration = factor(compoundConcentration),
    replicate = factor(replicate),
    dataset = factor(dataset)) %>%
  summary()
```
Make a few minor adjustments
```{R Clean Data}
# remove decoy proteins which contain # # # prefix
tppData <- tppData %>%
  filter(!grepl(":###[[:alnum:]]*###", uniqueID))
# remove proteins that do not have at least one unique peptide
tppData <- filter(tppData, uniquePeptideMatches >= 1)
# Remove proteins that only contain missing values
tppData <- tppData %>%
  filter(!is.na(relAbundance))
# remove all proteins which were not reproducibly observed with full melt
# curves in both treated groups.  # ten temperatures were recorded, so there
# should be 10 measurements
# tppData <- tppData %>%
#   group_by(dataset, uniqueID) %>%
#   mutate(n=n()) %>%
#   group_by(dataset) %>%
#   mutate(max_n = max(n)) %>%
#   filter(n == max_n) %>%
#   dplyr::select(-n,-max_n) %>%
#   ungroup
# count the remaining proteins in the dataset
tppData %>%
  distinct(dataset,uniqueID) %>%
  distinct %>%
  group_by(dataset) %>%
  tally %>%
  kable()
```
```{R Combine Datasets}
# two concentrations of dasatinib : 0.5 uM and 5 uM
# combine both datasets into one with 4 vehicle and two treated repl
tppData <- tppData %>%
  mutate(
    replicate = ifelse(
      dataset == "Dasatinib 5",
      yes = replicate + 2,
      no = replicate)) %>%
  mutate(dataset = gsub("0.5| 5", "", dataset))
# check result
tppData %>%
  distinct(dataset,replicate, compoundConcentration) %>%
  filter(compoundConcentration > 0) %>%
  dplyr::rename('drugconcentration (treatment groups)' = compoundConcentration) %>%
  kable()
```
```{R Filter by one protein}
stk4 <- filter(tppData, dataset == "Staurosporine", uniqueID == "STK4_IPI00011488")
# Pick a background protein (find the needle within tppData
# # select data
# select desired protein
stk4 %>% filter(compoundConcentration == 20, replicate == 1) %>%
  dplyr::select(-dataset) %>%
  kable(digits = 2)
# Generate a plot of the measurements
stk4_plot <- ggplot(stk4, aes(x = temperature, y = relAbundance)) +
  geom_point(
    aes(shape = factor(replicate),
        color = factor(compoundConcentration)),
    size = 2
  ) +
  ggtitle(stk4$uniqueID) +
  scale_color_manual(
    print(stk4$dataset[1]), values = c("#808080","#DA7F2D")
  ) +
  theme_bw()
stk4_plot
```
```{R Null Fit}
# Model fitting: create a function where null hypothesis fits a sigmoid curve over all data
fitSingleSigmoid <- function(x , y, start =c(Pl=0, a = 550, b = 10))
{
  try(
    nls(
      formula = y ~ (1-Pl)/(1+exp((b-a/x)))+Pl,
      start = start,
      data = list(x=x,y=y),
      na.action = na.exclude,
      algorithm = "port",
      lower = c(0.0,1e-5,1e-5),
      upper = c(1.5,15000,250),
      control = nls.control(maxiter = 50)),
    silent = TRUE
  )
}
nullFit <- fitSingleSigmoid(x = stk4$temperature, y = stk4$relAbundance)
summary(nullFit)
```
```{R View Null Predictions}
# append predictions and residuals in the same table
nullPredictions <- broom::augment(nullFit)
# view
nullPredictions %>%
  filter(x %in% c(46,49)) %>%
  kable()
```
```{R Plot Null Preditions}
# append and show fit in the plot
stk4$nullPrediction <- nullPredictions$.fitted
stk4$nullResiduals <- nullPredictions$.resid
# add null line
stk4_plot <- stk4_plot +
  geom_line(data = stk4, aes(y = nullPrediction))
stk4_plot
```
```{R Alternative}
# alternative hypothesis: separate melting curves per exp group
alternativePredictions <- stk4 %>%
  group_by(compoundConcentration) %>%
  do(
    {
      fit = fitSingleSigmoid(x = .$temperature, y= .$relAbundance)
      broom::augment(fit)
    }
  ) %>%
  ungroup %>%
  # rename columns to merge into a data frame
  dplyr::rename(
    alternativePredictions = .fitted,
    alternativeResiduals = .resid,
    temperature = x,
    relAbundance = y
  )
stk4 <- stk4 %>%
  left_join(
    alternativePredictions,
    by = c("relAbundance", "temperature", "compoundConcentration")
  ) %>%
  distinct()
# add the curves to the plot
stk4_plot <- stk4_plot +
  geom_line(
    data = distinct(
      stk4, temperature, compoundConcentration, alternativePredictions
    ),
    aes(
      y = alternativePredictions, color = factor(compoundConcentration)
    )
  ) +
  ggtitle(stk4$uniqueID)
stk4_plot
```
```{R Compute RSS}
# compute RSS values to quantify GOF relative to null
rssPerModel <- stk4 %>%
  summarise(
    rssNull = sum(nullResiduals^2),
    rssAlternative = sum(alternativeResiduals^2)) %>%
  data.frame()
kable(rssPerModel, digits = 4)
```
```{R}
# filter the data for the best compound
# get the original protein
DataC <- filter(tppData, dataset == "Staurosporine", uniqueID == "STK4_IPI00011488")
# get the treated and vehicle filtered data
Data0 <- tppData %>%
  subset(tppData$uniqueID %in% DataC$uniqueID) %>%
  filter(.$dataset==DataC$dataset[1]) %>%
  filter(compoundConcentration == 0 )
Data1 <- tppData %>%
  subset(tppData$uniqueID %in% DataC$uniqueID) %>%
  filter(.$dataset == DataC$dataset[1]) %>%
  filter(compoundConcentration > 0)
# Null hypothesis
z <- nplr(x = DataC$temperature, y = DataC$relAbundance, silent = FALSE, npars = 5,useLog=FALSE) 


xn <- (getXcurve(z))
yn <- as.numeric(getYcurve(z))



t <- unique(DataC$temperature)
# get fitted temperature values
veh <- data.frame(xn, yn)
datan <- veh %>% subset(t %in% xn)
# alt Hypothesis
x <- nplr(x = Data0$temperature, y = Data0$relAbundance, silent = FALSE, npars = 5)
x0 <- as.numeric(10^getXcurve(x))
y0 <- as.numeric(getYcurve(x))
veh <- data.frame(x0, y0)
data0 <- veh %>%
  subset(t %in% x0)
y <- nplr(x = Data1$temperature, y = Data1$relAbundance, silent = FALSE, npars = 5)
x20 <- as.numeric(10^getXcurve(y))
y20 <- as.numeric(getYcurve(y))
veh <- data.frame(x20, y20)
data20 <- veh %>% subset(t %in% x20)
# plot the fitted values
# plot the original values
orig0 <- Data0
# plot(orig0$temperature,orig0$relAbundance,main=Data0$uniqueID[1])
origt <- Data1
# plot(origt$temperature,origt$relAbundance,main=Data1$uniqueID[1])
data <- rbind(orig0, origt)
# Null hypothesis
colnames(datan) <- c('temperature','relAbundance')
colnames(data0) <- c('temperature','relAbundance')
colnames(data20)<- c('temperature','relAbundance')
# plot results
control <- ggplot(data, aes(x = temperature, y = relAbundance)) +
  geom_line(data = datan, color = 'Black')
control <- control +
  geom_point(aes(shape = factor(replicate),color = factor(compoundConcentration)))
control <- control +
  geom_line(data = data20, color = 'Orange') +
  geom_line(data = data0, color = 'Gray') +
  ggtitle(Data0$uniqueID)
print(control)

c<-NA
tr<-NA
c<-nplr(x = Data0$temperature, y = Data0$relAbundance, silent = FALSE, npars = 5,useLog=FALSE)
tr<-nplr(x = Data1$temperature, y = Data1$relAbundance, silent = FALSE, npars = 5,useLog=FALSE)


theme_set(theme_bw())
plot(c, pcol="grey40", lcol="gray1",  showInfl=TRUE,xlab= "" ,ylab="",ylim=c(0,1.1))
par(new=TRUE,mfrow = c(1, 1))
plot(tr, pcol="grey40", lcol="chocolate1", showInfl=TRUE, cex.main=0.5,xlab="Temperature",ylab="Relative Intensity",ylim=c(0,1.1))
box(which = "plot", lty = "solid")

#control params for prediction intervals
B<-c@pars$bottom
T_<-c@pars$top
xmid<-c@pars$xmid
b<-c@pars$scal
s<-c@pars$s

#bootstrap data

names(Data0)[names(Data0) == "temperature"] <- "C"
names(Data0)[names(Data0) == "relAbundance"] <- "I"

#bootstrap to generate random $Ivalues with variance on Intensity
BS<-NA
BA<-NA
BSvar<-list(NA)
BP<-list(NA)
temps<-NA
BC<-NA
DF<-NA
for (i in 1:100){
  BS<-sample_n(Data0,20,replace=TRUE)
  #generate mean intensities per $C value
  BA<-BS %>% group_by(C) %>% dplyr::summarise(n=n(),mean=mean(I,na.rm=TRUE),var=var(I,na.rm=TRUE))
  BA<-na.omit(BA)
  #generate intensity values
  BSvar[[i]]<-rnorm(length(BA$var),BA$mean, sqrt(BA$var))
  
}
BSvar<-unlist(BSvar)
names(BSvar)<-"I"
LOW<-NA
HI<-NA
n<-NA
T_1<-NA
n<-length(BSvar)
T_1 <- qt(.975, n-2)
stdErr<-as.numeric(c@stdErr[1])

#ypred


newy<-BSvar
#mean y pred
ybar<-mean(BSvar)
#generate CIs
CI<- T_1*(stdErr)*((1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
PI <- T_1*stdErr*sqrt((1+1/n+(newy - ybar)^2/sum((newy - ybar)^2)))

PLOW<-newy-PI
PHIGH<-newy+PI

CLOW<-newy-CI
CHIGH<-newy+CI
#Use the 5parameter equation to obtain concentrations


Data<-data.frame(as.numeric(BSvar),PLOW,PHIGH,CLOW,CHIGH)
names(Data)[1]<-"I"
BC<-Data %>% data.frame() %>% 
  dplyr::mutate(C=xmid - 1/b*
                  log10(((T_ - B)/(Data$I - B))^(1/s)-1)) %>% na.omit(.)

names(BC)[1]<-"I"
BC$group<-"vehicle"


#plot prediction intervals for treated 

#treated params for prediction intervals
B<-tr@pars$bottom
T_<-tr@pars$top
xmid<-tr@pars$xmid
b<-tr@pars$scal
s<-tr@pars$s

#bootstrap data

names(Data1)[names(Data1) == "temperature"] <- "C"
names(Data1)[names(Data1) == "relAbundance"] <- "I"

#bootstrap to generate random $Ivalues with variance on Intensity
BS<-NA
BA<-NA
BSvar<-list(NA)
BP<-list(NA)
temps<-NA
BC1<-NA
DF<-NA
for (i in 1:100){
  BS<-sample_n(Data1,20,replace=TRUE)
  #generate mean intensities per $C value
  BA<-BS %>% group_by(C) %>% dplyr::summarise(n=n(),mean=mean(I,na.rm=TRUE),var=var(I,na.rm=TRUE))
  BA<-na.omit(BA)
  #generate intensity values
  BSvar[[i]]<-rnorm(length(BA$var),BA$mean, sqrt(BA$var))
  
}
BSvar<-unlist(BSvar)
names(BSvar)<-"I"
LOW<-NA
HI<-NA
n<-NA
T_2<-NA
n<-length(BSvar)
T_2 <- qt(.975, n-2)
stdErr<-as.numeric(tr@stdErr[1])

#ypred
newy<-BSvar


#mean y pred
ybar<-mean(BSvar)
CI<- T_1*(stdErr)*((1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
PI <- T_1*stdErr*sqrt((1+1/n+(newy - ybar)^2/sum((newy - ybar)^2)))

PLOW<-newy-PI
PHIGH<-newy+PI

CLOW<-newy-CI
CHIGH<-newy+CI
#Use the 5parameter equation to obtain concentrations


Data1<-data.frame(as.numeric(BSvar),PLOW,PHIGH,CLOW,CHIGH)



names(Data1)[1]<-"I"
BC1<-NA
BC1<-Data1 %>% data.frame() %>% 
  dplyr::mutate(C=xmid - 1/b*
                  log10(((T_ - B)/(Data1$I - B))^(1/s)-1)) %>% na.omit(.)

names(BC1)[1]<-"I"
BC1$group<-"treated"

PIS<-NA
PIS<-rbind(BC,BC1)
#plot upper PI
theme_set(theme_bw())

PredPI<-ggplot(data = PIS,x=C,y=I)+geom_point(aes(x=C,y=I))+geom_line(aes(x=C,y=I,group=group)) +geom_ribbon(aes(x=C,ymin=PLOW,ymax=PHIGH,colour=group),alpha=0.2)+xlab("Temperature")+ylab("Relative Intensity")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+scale_x_continuous(limits=c(40,65))+scale_y_continuous(limits=c(0,1.2))
print(PredPI)
CPI<-NA
CPI<-PredPI+geom_ribbon(aes(x=C,ymin=CLOW,ymax=CHIGH,colour=group))
print(CPI)
#predict with 5 pl
pred1<-predict(c, interval="predict")

#summarize RSS values null for 5 PL
SEnull<-z@stdErr[[1]]
RSSnull<-SEnull*(length(DataC$temperature)-2)
#RSS values alt control
SEc<-c@stdErr[[1]]
RSScl<-SEc*(length(Data0$I)-2)
#RSS values alt treated
SEtr<-tr@stdErr[[1]]
RSStr<-SEtr*(length(Data0$I)-2)

results5<-data.frame(RSSnull,RSScl+RStr)
names(results5)<-c("RSSnull","RSSalt")

results5$RSSdiff<-results5$RSSnull-results5$RSSalt

#sigmoidal results

rssPerModel$RSSdiff<-rssPerModel$rssNull-rssPerModel$rssAlternative


#Hypothesis testing H0: 5 parameter model fits better than sigmoidal fit
dfa<-nrow(Data0)*3-5
dfn<-nrow(Data0)*3-4
F5<-dfa*(results5$RSSdiff-rssPerModel$RSSdiff)/dfn*rssPerModel$RSSdiff



```
