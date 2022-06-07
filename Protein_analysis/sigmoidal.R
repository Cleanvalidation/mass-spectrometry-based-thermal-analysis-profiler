##################################################
#Sigmoidal function with confidence intervals
###################################################


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

repeatFits <- function(x,y, start= c(Pl = 0, a = 550, b=10),
                       seed = NULL, alwaysPermute = FALSE, maxAttempts = 100){
  i <-0
  doFit <- TRUE
  doVaryPars <-alwaysPermute
  
  if(!is.null(seed)){
    set.seed(seed)
  }  
  while (doFit){
    startTmp <-start * (1+doVaryPars*runif(1, -0.5,0.5))
    m <- fitSingleSigmoid(x=x, y=y, start = startTmp)
    i <- i +1
    doFit <- inherits(m,"try-error") & i < maxAttempts
    doVaryPars <- TRUE
  }
  return(m)
}

computeRSS <- function(x,y,start = c(Pl = 0, a = 550, b=10),seed = NULL,
                       alwaysPermute = FALSE,maxAttempts = 50){
  #start model fitting
  fit <- repeatFits(x=x, y=y,start=start,seed=seed,alwaysPermute=alwaysPermute,maxAttempts = maxAttempts)
  if(!inherits(fit,"try-error")){
    #if model fit converged, calculate RSS and parameters
    resid <-residuals(fit)
    rss <- sum(resid^2,na.rm = TRUE)
    fittedValues <- sum(!is.na(resid))
    params <- coefficients(fit)
    #Compute Tm for comparison (not needed for HT)
    tm <- params[["a"]]/(params[["b"]] - log(1-params[["Pl"]])/(1/2 - params[["Pl"]]-1))
  } else {
    #if model did not converge, return default vals
    rss <- NA
    fittedValues <-0
    params <- c(Pl=NA,a=NA,b=NA)
    tm <- NA
  }
  out <- tibble(rss=rss, fittedValues = fittedValues, tm= tm,
                a = params[["a"]],b = params[["b"]], Pl = params[["Pl"]])
  return(out)
  
  
}
# #GoF by RSS comparison for each protein
# #RSS difference
# 
computeRSSdiff <- function(x,y,treatment,maxAttempts = 50, repeatsIfNeg = 100){
  rssDiff <- -1
  repeats <- 0
  alwaysPermute <- FALSE
  
  start0 = start1 <- c(Pl=0,a=550,b=10)
  
  while((is.na(rssDiff) | rssDiff<0) & repeats <= repeatsIfNeg){
    
    nullResults <- computeRSS(x=x, y=y, start = start0, seed=repeats,
                              maxAttempts = maxAttempts,
                              alwaysPermute = alwaysPermute)
    
    altResults <- tibble(x,y,treatment) %>%
      group_by(treatment) %>%
      purrr::map_dfr({
        fit = computeRSS(x=.$x, y = .$y,start = start1, seed=repeats,
                         maxAttempts = maxAttempts,
                         alwaysPermute = alwaysPermute)
        
      }) 
    rss0 <- nullResults$rss
    rss1 <-sum(altResults$rss)#combined alternative RSS values
    rssDiff <- rss0-rss1#difference between null and combined alternative RSS values
    
    if(is.na(rssDiff) | rssDiff <0){
      repeats <- repeats +1
      alwaysPermute <- TRUE
      start1 <- c(Pl = nullResults[["Pl"]],a = nullResults[["a"]],b=nullResults[["b"]])
      
    }
  }
  
  n0 <-nullResults$fittedValues
  n1 <-sum(altResults$fittedValues)
  
  tm <- altResults %>%
    mutate(key = paste0("tm_",treatment)) %>%
    dplyr::select(key,tm) %>% spread(key,tm)
  
  out <- tibble(rss0,rss1,rssDiff,n0,n1,repeats) %>%
    cbind(tm)
  return(out)
}
sigCI <- function(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...)
{
  method <- match.arg(method)
  
  format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
          "%")
  
  ## Taken from confint.nls
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- seq_along(pnames)
  if (is.numeric(parm)) 
    parm <- pnames[parm]
  
  ## Taken from confint.default and modified slightly to use t-distribution
  asCI <- function(object, parm, level)
  {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    #        pct <- stats:::format.perc(a, 3)
    pct <- format.perc(a, 3)
    fac <- qt(a, df.residual(object))
    
    parmInd <- match(parm, pnames)
    ci <- array(NA, dim = c(length(parmInd), 2), dimnames = list(parm, pct))
    ses <- sqrt(diag(vcov(object)))[parmInd]
    ci[] <- cf[parmInd] + ses %o% fac
    ci
  }
  
  ## Taken from confint.nls
  asProf <- function(object, parm, level)
  {
    
    utils::flush.console()
    object <- profile(object, which = parm, alphamax = (1 - level)/4)
    confint(object, parm = parm, level = level, ...)    
  }
  
  switch(method, asymptotic = asCI(object, parm, level), profile = asProf(object, parm, level))
}

sigC<-function(df_,Protein,Peptide=FALSE,stats=FALSE){
  Protein<-as.character(Protein)
  df_$C<-as.numeric(as.character(df_$C))
  if(isTRUE(stats)){
    DFN<-df_
    df_1<-df_%>% dplyr::filter(treatment=="treated")
    df_<-df_%>% dplyr::filter(treatment=="vehicle")
  }else{
    DFN<-df_%>% dplyr::filter(uniqueID %in% Protein)
    df_1<-df_%>% dplyr::filter(uniqueID%in%Protein,treatment=="treated")
    df_<-df_%>% dplyr::filter(uniqueID%in%Protein,treatment=="vehicle")
  }
  
  if(isTRUE(Peptide)){
    df_<-df_ %>% dplyr::rename("value"="I3","temperature"="C")
    df_1<-df_1 %>% dplyr::rename("value"="I3","temperature"="C")
    DFN<-DFN %>% dplyr::rename("value"="I3","temperature"="C")
  }else{
    df_<-df_ %>% dplyr::rename("value"="I","temperature"="C")
    df_1<-df_1 %>% dplyr::rename("value"="I","temperature"="C")
    DFN<-DFN %>% dplyr::rename("value"="I","temperature"="C")
  }
  
  if(isTRUE(stats)){
    nlm2<-function(x){
      nlm2<-x %>% 
        dplyr::mutate(fit=list(cetsa_fit2(.,norm=FALSE)),
                      Tm=as.numeric(coef(fit[[1]][[1]])[3]))
      
      return(nlm2) 
    }
    df_<-df_ %>% dplyr::group_split(uniqueID)
    df_1<-df_1 %>% dplyr::group_split(uniqueID)
    DFN<-DFN %>% dplyr::group_split(uniqueID)
    if (.Platform$OS.type=="windows"){
      nlm2<-parallel::mclapply(df_,nlm2)
    }else{
      nlm2<-parallel::mclapply(df_,nlm2,mc.cores=availableCores())
    }
  }else{
    nlm2<-df_  %>%   
      dplyr::mutate(fit=list(cetsa_fit2(.,norm=FALSE)))
    nlm2$Tm<-NA
    nlm2$Tm<-as.numeric(coef(nlm2$fit[[1]][[1]])[3])
  }
  if(class(nlm2$fit[[1]][[1]])=='try-error'){
    warning("the function could not fit the data")
  }
  
  
  
  #ready confidence intervals for plot
  result<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(.$fit[[1]][[2]]$lwr.conf),
                                                        HI = list(.$fit[[1]][[2]]$upr.conf),
                                                        CP=list(.$fit[[1]][[2]]$temperature),
                                                        IP=list(.$fit[[1]][[2]]$value),
                                                        dof=summary(.$fit[[1]][[1]])$df[2])# lower CI
  result<-result %>% dplyr::rename("I"="value","C"="temperature")
  
  result <-result %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]][[1]]))
  
  
  nlm1<-list()
  CT<-list()
  dfc<-list()
  
  #sigmoidal fit for treated
  nlm2<-df_1  %>%   
    dplyr::mutate(fit=list(try(cetsa_fit2(.,norm=FALSE))))
  nlm2$Tm<-NA
  nlm2$Tm<-as.numeric(coef(nlm2$fit[[1]][[1]])[3])
  
  if(class(nlm2$fit[[1]][[1]])=='try-error'){
    warning("the function could not fit the data")
  }
  
  #ready confidence intervals for plot
  result1<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(.$fit[[1]][[2]]$lwr.conf),
                                                         HI = list(.$fit[[1]][[2]]$upr.conf),
                                                         CP=list(.$fit[[1]][[2]]$temperature),
                                                         IP=list(.$fit[[1]][[2]]$value),
                                                         dof=summary(.$fit[[1]][[1]])$df[2])# lower CI
  result1<-result1 %>% dplyr::rename("I"="value","C"="temperature")
  
  result1 <-result1 %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]][[1]]))
  
  
  
  return(rbind(Pred,Pred1))
}
sigfit<-function(SigF,Peptide=FALSE){
  
  #
  if(isTRUE(Peptide)){
    
    Pred<-SigF %>%
      subset(treatment=="vehicle") %>% 
      dplyr::select(uniqueID,Annotated_Sequence,treatment ,C,I,CP,IP,CC,Tm,rss,
                    LOW,HI,sample_name,sample_id,missing_pct)
    
    Pred1<-SigF%>%
      subset(treatment=="treated") %>% 
      dplyr::select(uniqueID,Annotated_Sequence,treatment ,C,I,CP,IP,CC,Tm,rss,
                    LOW,HI,sample_name,sample_id,missing_pct)
    
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    
    #unnest data
    Pred<-Pred %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred1<-Pred1 %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred$AUC<-round(pracma::trapz(Pred$IP),2)
    Pred1$AUC<-round(pracma::trapz(Pred1$IP),2)
    Pred1$dAUC<-abs(Pred1$AUC[1]-Pred$AUC[1])
    
    Pred1$dAUC<-pracma::trapz(Pred1$IP[(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))-1):(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))+1)])-pracma::trapz(Pred$IP[(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))-1):(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))+1)])
    Pred1$dAUC<-abs(round(Pred1$dAUC[1],3))
    
    #max replicates
    roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
      if(length(x) != 1) stop("'x' must be of length 1")
      10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
    }
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    
    #append missing value data
    if(isTRUE(Peptide)){
      
      Pred<-Pred %>% distinct(.) %>% dplyr::group_by(uniqueID,treatment,sample_id,Annotated_Sequence) %>% dplyr::group_split(.)
      Pred<-purrr::map(Pred,function(x) x %>% dplyr::mutate(missing_v=100*(roundUpNice(length(unique(x$I)))-length(unique(x$I)))/roundUpNice(length(unique(x$I)))))
      Pred<-dplyr::bind_rows(Pred)
      Pred<-Pred %>% dplyr::group_split(uniqueID,treatment)
      Pred<-purrr::map(Pred,function(x) x %>% dplyr::mutate(missing_v=getmode(x$missing_v)))
      
      Pred1<-Pred1 %>% distinct(.) %>% dplyr::group_by(uniqueID,treatment,sample_id,Annotated_Sequence) %>% dplyr::group_split(.)
      Pred1<-purrr::map(Pred1,function(x) x %>% dplyr::mutate(missing_t=100*(roundUpNice(length(unique(x$I)))-length(unique(x$I)))/roundUpNice(length(unique(x$I)))))
      Pred1<-dplyr::bind_rows(Pred1)
      Pred1<-Pred1 %>% dplyr::group_split(uniqueID,treatment)
      Pred1<-purrr::map(Pred1,function(x) x %>% dplyr::mutate(missing_t=getmode(x$missing_t)))
      
      
      Pred<-dplyr::bind_rows(Pred)
      Pred1<-dplyr::bind_rows(Pred1)
    }else{
      Pred$missing_v<-round(Pred$missing_v[1],0)
      Pred1$missing_t<-round(Pred1$missing_t[1],0)
    }
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x=C,y=IP,color=treatment))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred1,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      ggplot2::annotate("text", x=45, y=-0.65, label= paste("missing: ",round(Pred$missing_v[1],0),"%"),colour="#00BFC4")+ 
      ggplot2::annotate("text", x=45, y=-0.75, label= paste("missing: ",round(Pred1$missing_t[1],0),"%"),colour="#F8766D")+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = 37, xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +
      ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = 37, xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2)+ylim(-0.8,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
    print(P1)
  }else{
    Pred<-SigF %>%
      subset(treatment=="vehicle") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,Cell_Component,Bio_Process,Coverage,MW_kDa,
                    LOW,HI,sample_name,sample_id,missing_pct,rank)
    Pred1<-SigF%>%
      subset(treatment=="treated") %>% 
      dplyr::select(uniqueID,treatment ,C,I,CP,IP,CC,Tm,rss,Cell_Component,Bio_Process,Coverage,MW_kDa,
                    LOW,HI,sample_name,sample_id,missing_pct,rank)
    
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #unnest data
    Pred<-Pred %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    Pred1<-Pred1 %>% tidyr::unnest(cols=c(CP,IP,LOW,HI))
    
    Pred$AUC<-round(pracma::trapz(Pred$IP),2)
    Pred1$AUC<-round(pracma::trapz(Pred1$IP),2)
    Pred1$dAUC<-abs(Pred1$AUC[1]-Pred$AUC[1])
    
    
    Pred1$dAUC<-pracma::trapz(Pred1$IP[(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))-1):(which(abs(Pred1$IP-0.5)==min(abs(Pred1$IP-0.5)))+1)])-pracma::trapz(Pred$IP[(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))-1):(which(abs(Pred$IP-0.5)==min(abs(Pred$IP-0.5)))+1)])
    Pred1$dAUC<-abs(round(Pred1$dAUC[1],3))
    
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x=C,y=IP,color=treatment))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred1,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      ggplot2::annotate("text", x=45, y=-0.35, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=45, y=-0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=45, y=-0.55, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1],Pred1$sample_name[1]))+
      ggplot2::annotate("text", x=45, y=-0.65, label= paste("missing: ",round(Pred$missing_v[1],0),"%"),colour="#00BFC4")+ 
      ggplot2::annotate("text", x=45, y=-0.75, label= paste("missing: ",round(Pred1$missing_t[1],0),"%"),colour="#F8766D")+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = 37, xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +
      ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = treatment))+
      ggplot2::geom_line(Pred,mapping=ggplot2::aes(x=CP,y=IP),alpha = 0.2 ) +
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(x=CP,y=IP,ymin = LOW, ymax = HI ,fill=treatment),alpha = 0.2 ) +
      annotate("text",
               x = round(Pred$Tm[1],1),
               y = -0.05,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = 37, xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2)+ylim(-0.8,1.5)+xlim(37,68)+theme(legend.position="bottom")
    
    
    print(P1)
    
    
  }
  
}
