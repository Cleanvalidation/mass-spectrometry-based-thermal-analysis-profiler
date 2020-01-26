library(tidyverse)
library(minpack.lm)
library(dplyr)
library(rlist)
library(ggplot2)
library(data.table)
library(tidyverse)
library(knitr)
library(ggthemes)
library(gridExtra)
library(grid)
library(readxl)
library(nls2)
library(tidyr)
library(minpack.lm)
library(nlstools)
library(pkgcond)
library(rlist)
theme_set(theme_bw())

#' Read in data
#'
#' Read in Excel file and apply minimal pre-processing
#'
#' @param  f.  Path of Excel file output from Proteome Discoverer
#' @return a dataframe containing extracted information
#'
#' @importFrom readxl read_excel
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom stringr str_extract
#' @export
read_cetsa <- function(f) {
  df.raw <- read_excel(f)
  
  df <- df.raw %>%
    dplyr::select(Accession, dplyr::starts_with('Abundance')) %>%
    dplyr::select(-contains('(Grouped)')) %>%
    gather('id', 'value', -Accession) %>%
    dplyr::mutate(sample = str_extract(id, '(?<=Abundance: ).*(?=:)')) %>%
    dplyr::mutate(temp_ref = str_extract(id, '(?<=: )([0-9]+[N|C]?)'))
  
  return(df)
}


#' Clean data
#'
#' Clean CETSA data
#'
#' @param df.  Data frame returned by read_cetsa
#' @param temperatures.  Data frame of temperatures related to TMT tags - columns =
#'     temp_ref, temperature
#' @param samples.  Optional data frame of sample names - columns = sample, sample_name
#' @param separator.  Character used to separate parts of sample name.  Sample name is constructed as
#'     sampleroot_A_B where sampleroot is the name of the sample, A is the number of the biological
#'     replicate and B is the number of the technical replicate.
#'
#' @return a data frame of clean data
#'
#' @import dplyr
#'
#' @export
clean_cetsa <- function(df, temperatures = NULL, samples = NULL) {
  if (is.null(temperatures)) warning('No temperature data')
  if (!is.null(samples)) {
    df <- df %>%
      dplyr::left_join(samples, by = c('sample' = 'sample_id'))
  }
  if (!is.null(temperatures)) {
    df <- df %>%
      dplyr::left_join(temperatures, by = 'temp_ref')
    
  } else {
    df <- df %>%
      dplyr::rename(temperature = temp_ref)
  }
  df <- df %>%
    dplyr::select(Accession, sample, temperature,value) %>%
    dplyr::filter(!is.na(temperature),!is.na(value)) %>%
    dplyr::group_by(Accession,sample) %>%
    dplyr::mutate(value = value / value[temperature == min(temperature)]) %>% ungroup()
  
  return(df)
}


#' normalize CETSA data
#'
#' Normalize data according to Pelago paper.
#' Normalizaion works by establishing a subset of curves which pass a set of criteria
#' (see below) and then calculating the median fold change through this set.  A curve is
#' then fit through these points and used to generate correction factors.
#'
#' Curve selection criteria:
#'   \itemize{
#'   \item All 10 temperature points present
#'   \item T7/T1 between 0.4 and 0.6
#'   \item T9/T1 < 0.3
#'   \item T10/T1 < 0.2
#' }
#'
#' @param df. A cleaned CETSA data frame
#' @param temperatures. A list of experimental temperatures
#'
#' @return normalized CETSA data frame
#'
#' @import dplyr
#' @import tidyr
#' @import nls2
#' @export
normalize_cetsa <- function(df, temperatures) {
  temperatures <- sort(temperatures)
  df.jointP <- df %>%
    dplyr::group_by(Accession, sample) %>%
    dplyr::filter(n() == 10) %>%  ## all 10 temperatures are present
    dplyr::summarise(T7 = value[temperature == temperatures[[7]]] / value[temperature == temperatures[[1]]],
                     T9 = value[temperature == temperatures[[9]]] / value[temperature == temperatures[[1]]],
                     T10 = value[temperature == temperatures[[10]]] / value[temperature == temperatures[[1]]]) %>% filter(T7 >= 0.4 & T7 <= 0.6 & T9 < 0.3 & T10 < 0.2)%>% ungroup()
  
  ## split by sample group and filter
  l.bytype <- split.data.frame(df.jointP, df.jointP$sample)
  
  ## determine which group contains the greatest number of curves and use this for normalization
  n.filter <- lapply(l.bytype, nrow)
  df.normP <- l.bytype[[which.max(n.filter)]]
  norm.accessions <- df.normP$Accession
  
  ## calculate median for each sample group
  
  df.mynormset <- df %>% base::subset(Accession %in% norm.accessions)
  
  df.median <- df %>%
    group_by(sample,temperature) %>%
    dplyr::summarise(value = median(value))
  
  
  ## fit curves to the median data
  df.fit <- df.median %>%
    group_by(sample) %>% 
    dplyr::do(fit = cetsa_fit(d = ., norm = FALSE))
  ## calculate the fitted values
  d<-list.count(df.fit$fit)
  df.fittedVals<-0
  for (i in 1:d){
    if (is.na(df.fit$fit[i])) {
      df.fittedVals[i][[1]] <-NA_real_
    } else {
      df.fittedVals[i] <- as.data.frame(predict(df.fit$fit[[i]]))#interesting, plot(df.fit$fit[[1]])
    } 
    
  }
  
  df.fittedVals<- df.fittedVals %>% as.data.frame()
  names(df.fittedVals) <- df.fit$sample
  
  df.fittedVals <- df.fittedVals %>% gather()
  colnames(df.fittedVals)<-c('sample','fitted_values')
  
  df.fittedValst<-df.fittedVals %>% group_by(sample)
  
  ## calculate ratios between the fitted curves and the median values
  df.out <- df.median %>%
    data.frame(dplyr::full_join(df.fittedValst,df.median,'sample')) %>%
    dplyr::mutate(correction = df.fittedValst$fitted_values / value) %>% select('sample','temperature','value','fitted_values','correction')
  
  ## apply normalization factor to data
  df <- df %>% dplyr::left_join(df.out %>% select(sample,temperature,correction), by = c('sample', 'temperature')) %>%
    dplyr::mutate(norm_value = value * correction)
  
  return(df)
}

#
#' fit curves to CETSA data
#'
#' fit curves to normalized or non-normalized CETSA data
#'
#' @param df. data frame containing CETSA data
#' @param normalized_data  Boolean.  If true then use normalized data otherwise use non-normalized
#' @param n_cores.  Number of cores for parallel processing
#' @param separator.  Character used to separate sample name from replicate.  If NULL then
#'     replicates are treated as individual samples.
#'
#' @return list
#'
#' @import nls2
#' @import dplyr
#' @importFrom tidyr separate
#' @import multidplyr
#' @importFrom tibble rowid_to_column
#'
#' @export
curves_cetsa <- function(df, normalized_data = TRUE, n_cores = 1, separator = NULL) {
  cluster <- create_cluster(cores = n_cores)
  
  if (!is.null(separator)) {
    df <- df %>%
      tidyr::separate(sample, c('sample', 'replicate'), sep = separator, convert = TRUE)
  }
  
  df.curve <- df %>%
    dplyr::group_by(sample, Accession) %>%
    multidplyr::partition(sample, Accession, cluster=cluster) %>%
    multidplyr::cluster_assign_each('normalized_data', normalized_data) %>%
    multidplyr::cluster_assign_each('fit.cetsa', fit.cetsa) %>%
    multidplyr::cluster_assign_each('cetsa_fit', cetsa_fit) %>%
    multidplyr::cluster_library('nls2') %>%
    dplyr::do(fit = cetsa_fit(d = ., norm = normalized_data)) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    tibble::rowid_to_column('ref') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Tm = Tm(fit)) %>%
    dplyr::select(ref, Accession, sample, fit, Tm)
  return(df.curve)
}


#' calculate CETSA statistics
#'
#' calculate CETSA statistics - top and bottom plateau gradient and goodness of fit
#'
#' @param df.  data frame containing curve-fitted CETSA data
#' @param plateau_temps.  Vector containing plateau values (T1A, T1B, T2A, T2B)
#'
#' @return data frame containing curves with stats
#'
#' @import dplyr
#' @importFrom tibble as.tibble
#'
#' @export
stats_cetsa <- function(df, plateau_temps = c(37,40,64,67)) {
  
  df <- df%>%
    mutate(gof = gof(fit))
  
  responses <- apply(df, 1, function(x) {
    tryCatch({
      predict(x['fit'][[1]], list(t = plateau_temps))
    }, error = function(e) {
      rep(NA, length(plateau_temps))
    })
  })
  df.responses <- as.data.frame(as.matrix(t(responses)))
  df.responses$ref <- as.numeric(row.names(df.responses))
  df.responses$slopeStart <- apply(df.responses, 1, function(x) abs(x[2] - x[1]) / (plateau_temps[2] - plateau_temps[1]))
  df.responses$slopeEnd <- apply(df.responses, 1, function(x) abs(x[4] - x[3]) / (plateau_temps[4] - plateau_temps[3]))
  
  df <- cbind(df, slope_start=df.responses$slopeStart, slope_end=df.responses$slopeEnd)
  return(as.tibble(df))
}


#' tag CETSA curves
#'
#' Apply filters and tag curves
#'
#' @param df.  data frame containing curve-fitted CETSA data
#' @param tag_fit.  Boolean.  Tag if poor curve fit
#' @param tag_range.  Boolean.  Tag if Tm is outside experimental temperature range
#' @param tag_gof.  Boolean.  Tag if poor goodness-of-fit
#' @param tag_plateau.  Boolean.  Tag if plateau(s) not flat
#' @param plateau_filter.  Filter value for tag_plateau (default = 0.01)
#' @param gof_filter.  Filter value for tag_gof (default = 0.085)
#'
#' @return data frame containing tagged results
#'
#' @import dplyr
#' @importFrom purrr map
#'
#' @export
tag_cetsa <- function(df, tag_fit = TRUE, tag_range = TRUE, tag_gof = TRUE, tag_plateau = TRUE, plateau_filter = 0.01, gof_filter = 0.085) {
  
  df <- df %>% ungroup()
  
  if (tag_fit) {
    df <- df %>%
      mutate(tag_fit = if_else(is.na(fit), 'bad fit', NA_character_))
  }
  
  if (tag_range) {
    df <- df %>%
      mutate(temps = map(df$fit, c('data', 't'))) %>%
      rowwise() %>%
      mutate(t_min = min(temps, na.rm=T)) %>%
      mutate(t_max = max(temps, na.rm=T)) %>%
      ungroup() %>%
      mutate(tag_range = if_else(is.na(fit), 'no data',
                                 if_else(Tm < t_min | Tm > t_max, 'outside experimental range', NA_character_))) %>%
      select (-temps, -t_min, -t_max)
  }
  
  if (tag_gof) {
    df <- df %>%
      mutate(tag_gof = if_else(gof >= gof_filter, 'poor gof', NA_character_))
  }
  
  if (tag_plateau) {
    df <- df %>%
      rowwise() %>%
      mutate(tag_plateau = if_else(min(slope_start, slope_end) >= plateau_filter, 'poor plateau - low and high',
                                   if_else(slope_start >= plateau_filter, 'poor plateau - low',
                                           if_else(slope_end >= plateau_filter, 'poor plateau - high', NA_character_))))
  }
  
  ## add tag_summary column
  use_cols <- names(df)[which(startsWith(names(df), 'tag_'))]
  if (length(use_cols) > 0) {
    df$tag_summary <- apply(df, 1, function(x) {
      cols <- sapply(use_cols, function(y) x[y])
      cols <- cols[!is.na(cols)]
      paste0(cols, collapse = '; ')
    })
    
    df$tagged <- nchar(df$tag_summary) > 0
  }
  
  return(df)
}


#' extract curve fit parameters from a CETSA data frame
#'
#' extract curve fit parameters from a CETSA data frame
#'
#' @param df.  Data frame containing curve fit column
#' @param fit.  Curve fit column name (default = "fit")
#' @param remove_fit.  Boolean.  If true then remove the fit column
#'
#' @import dplyr
#'
#' @export
extract_parameters <- function(df, fit_col = 'fit', remove_fit = TRUE) {
  fit_col_ref <- which(names(df) == fit_col)
  if (length(fit_col_ref) == 0) {
    message('fit column not found')
    return(NULL)
  }
  
  # get variable names
  fit_types <- sapply(df[[fit_col]], class)
  first_match <- which(!fit_types == 'logical')[1]
  param_names <- names(df[[fit_col]][[first_match]]$m$getPars())
  
  df <- df %>%
    rowwise() %>%
    mutate(params = curve_parameters(!!as.name(fit_col)))
  
  # rearrange column order
  if (remove_fit) {
    df <- df[, c(1:(fit_col_ref-1), ncol(df), (fit_col_ref+1):(ncol(df)-1))]
  } else {
    df <- df[, c(1:fit_col_ref, ncol(df), (fit_col_ref+1):(ncol(df)-1))]
  }
  
  return(df)
}


#' convert to wide format
#'
#' convert to wide format
#'
#' @param df. Data frame containing Accession, sample and params columns
#'
#' @return a wide data frame
#'
#' @import dplyr
#' @importFrom tidyr spread
#'
#' @export
long_to_wide_cetsa <- function(df) {
  df %>%
    select(Accession, sample, params) %>%
    spread(sample, params)
}


#' curve fitting equation
#'
#' CETSA curve fitting equation
#' Curve equation \eqn{y = \frac{(1 - p)}{(1 + e^{(-k(\frac{1}{t} - \frac{1}{m}))})} + p}
#'
#' @param p.  Curve parameter
#' @param k.  Curve parameter
#' @param m.  Curve parameter
#' @param t.  Curve variable
#'
#' @return position on curve at t
#'
fit.cetsa <- function(p, k, m, t) {
  (1 - p)/(1 + exp(-k*(1/t - 1/m))) + p
}


#' curve fitting function
#'
#' CETSA curve fitting function
#' Fit a curve to a set of data
#'
#' @param d.  data frame containing temperature and value or norm_value columns
#' @param norm.  Boolean.  If true then use norm_value column otherwise use value column
#'
#' @return nls2 curve model
#'
#' @import nls2
#'
## fit curve using mstherm equation
cetsa_fit <- function(d, norm = FALSE) {
  if (sum(!is.na(d$value)) < 2 | sum(!is.na(d$temperature)) < 2) return(NA)
  result = tryCatch({
    if (!norm) {
      myData <- list(t = d$temperature, y = d$value)
    } else {
      myData <- list(t = d$temperature, y = d$norm_value)
    }
    fine_start <- expand.grid(p=c(0,0.3),k=seq(0,8000,by=200),m=seq(30,80,by=10))
    new_start <- nls2(y ~ fit.cetsa(p, k, m, t),
                      data = myData,
                      start = fine_start,
                      algorithm = "brute-force",
                      control = nls.control(warnOnly=T,maxiter=5000))
    nls2(y ~ fit.cetsa(p, k, m, t),
         data = myData,
         start = new_start,
         control = nls.control(warnOnly=F),
         algorithm = "port",
         lower = c(0, 1, 10),
         upper = c(0.4, 100000, 100))
  }, error = function(err) {
    return(NA)
  })
  return(result)
}


#' return Tm from CETSA curve
#'
#' return Tm from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return Tm.
#'
Tm <- function(f) {
  if (length(f) == 1) return(NA)
  f$m$getPars()[['m']]
}


#' return gof from CETSA curve
#'
#' return gof from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return gof.
#'
gof <- function(f) {
  if (length(f) == 1) return(NA)
  sigma(f)
}


#' return parameters from CETSA curve
#'
#' return parameters from CETSA curve
#'
#' @param f. fitted curve
#'
#' @return json formatted string of parameters
#'
#' @importFrom jsonlite toJSON
#'
curve_parameters <- function(f) {
  if (length(f) == 1) return(NA)
  as.character(toJSON(c(f$data, params=list(as.list(f$m$getPars())))))
}


#' plot a series of CETSA curves
#'
#' plot a series of CETSA curves
#'
#' @param r. Subset of rows from CETSA table
#'
#' @import RColorBrewer
#' @importFrom jsonlite fromJSON
#'
#' @export
plot_cetsa <- function(r) {
  ## plot a series of curves
  cols <- brewer.pal(length(r)-1, 'Dark2')
  legend_inc <- c()
  firstplot <- TRUE
  for (c in 2:length(r)) {
    if (!is.na(r[[c]])) {
      params <- fromJSON(r[[c]])
      
      df <- data.frame(t = params$t, y = params$y) %>%
        dplyr::group_by(t) %>%
        dplyr::summarise(av = mean(y, na.rm = T), sd = ifelse(n() == 1, 0, sd(y, na.rm = T)))
      
      ## plot experimental points
      if (firstplot) {
        plot(x = df$t, y = df$av, col = cols[c-1], main = r[[1]], xlab="Temperature", ylab="Normalized Response", cex.axis = 0.8, cex.main = 0.8, cex = 0.5, ylim = c(0, 1.2))
        firstplot <- FALSE
      } else {
        points(x = params$t, y = params$y, col = cols[c-1], cex = 0.5)
      }
      
      ## error bars
      if(!all(df$sd == 0)) {
        arrows(df$t, df$av-df$sd, df$t, df$av+df$sd, col = cols[c-1], length=0.05, angle=90, code=3)
      }
      
      ## plot fitted curve
      t_range <- seq(range(params$t)[1], range(params$t)[2], by = 0.2)
      y_pred <- fit.cetsa(p=params$params$p, k=params$params$k, m=params$params$m, t_range)
      lines(t_range, y_pred, col = cols[c-1])
      
      ## add line through Tm
      abline(v = params$params$m, col = cols[c-1], lty = 2, lwd = 0.5)
      
      ## add to legend
      legend_inc <- c(legend_inc, c)
    }
  }
  ## include legend
  if (!firstplot) {
    legend("bottomleft",
           legend = sapply(names(r)[legend_inc], function(x) {paste0(x, ': ', round(params$params$m, 2), '°C')}),
           pch = rep(1, length(legend_inc)),
           col = cols[legend_inc - 1],
           cex = 0.4)
  }
}
#find patterns in data
find_pat = function(pat, x) 
{
  ff = function(.pat, .x, acc = if(length(.pat)) seq_along(.x) else integer(0L)) {
    if(!length(.pat)) return(acc)
    
    if(is.na(.pat[[1L]])) 
      Recall(.pat[-1L], .x, acc[which(is.na(.x[acc]))] + 1L)
    else 
      Recall(.pat[-1L], .x, acc[which(.pat[[1L]] == .x[acc])] + 1L)
  }
  
  return(ff(pat, x) - length(pat))
}  

##################################

# read in the the STK4 data
d <- readRDS("tppData_preprocessed.Rds")
#standardize names
d<-dplyr::rename(d,C=temperature,I=relAbundance,CC=compoundConcentration)
#prepare a list of datasets
datalist<-unique(d$dataset)
#prepare a list of concentration values
conclist<-unique(d$CC)
#prepare a list of proteins
protlist<-unique(d$uniqueID)

###selection of sample compound
d<-d %>%
  subset(uniqueID == "STK4_IPI00011488") %>%
  dplyr::filter(dataset=="Staurosporine" )
#subset control
Data0<-d %>%
  subset(uniqueID == "STK4_IPI00011488") %>%
  dplyr::filter(CC == 0,dataset=="Staurosporine" )

#subset treated
Data1<-d %>%
  subset(uniqueID == "STK4_IPI00011488") %>%
  dplyr::filter(CC > 0,dataset=="Staurosporine" )
#null hypothesis
DF<-d
DF<-DF[order(DF$C),]
#Alternative hypothesis
d<-Data0
d1<-Data1

#boostrap function
#returns bootstrapped intensities
#n is the sampling window
#N is the number of iterations
BStrap<-function(Data0,n,N){
  
  BS<-NA
 # BS1<-NA
  BA<-NA
 # BA1<-NA
  BSvar<-list(NA)
  #BSvar1<-list(NA)
  #Bootstrap
  for (i in 1:N){
    BS<-dplyr::sample_n(data.frame(Data0),n,replace=TRUE)
  #  BS1<-sample_n(data.frame(Data1),n,replace=TRUE)
    #generate mean intensities per $C value
    BA<-BS %>% dplyr::group_by(C) %>% dplyr::summarise(n=n(),mean=mean(I,na.rm=TRUE),var=var(I,na.rm=TRUE))
   # BA1<-BS1 %>% group_by(C) %>% dplyr::summarise(n=n(),mean=mean(I,na.rm=TRUE),var=var(I,na.rm=TRUE))
    
    #generate intensity values
    BSvar[[i]]<-rnorm(length(BA$var),BA$mean, sqrt(BA$var))
   # BSvar1[[i]]<-rnorm(length(BA1$var),BA1$mean, sqrt(BA1$var))
    
  }
  
  BSvar<-unlist(BSvar)
  #BSvar1<-unlist(BSvar1)
  
  BSvar<-na.omit(BSvar)
  
}

BSvar<-suppress_warnings(BStrap(Data0,20,100))
BSvar1<-suppress_warnings(BStrap(Data1,20,100))
#set up line region function
#returns lineregion column
#d is the alternative data
#Data0 is the same alternative data which will be modified for line region
DLR<-function(d,Data0){
  
  #df_1<-NA
  df <- d %>% dplyr::group_by(C) %>% dplyr::summarise(n=n(),I=mean(I))
 # df_1 <- d1 %>% group_by(C) %>% dplyr::summarise(n=n(),I=mean(I))
  
  df<-df[order(df$C),]
  #df_1<-df_1[order(df_1$C),]
  df$DeltaIntensity <- rep(0,nrow(df))
 # df_1$DeltaIntensity <- rep(0,nrow(df_1))
  # iterate through and find the change from the last value
  for(i in 2:nrow(df)){
    df$DeltaIntensity[i] <- abs(df$I[i] - df$I[i-1])
  }
  # for(i in 2:nrow(df_1)){
  #   df_1$DeltaIntensity[i] <- abs(df_1$I[i] - df_1$I[i-1])
  # }
  # get the scores of the changes
  df$DeltaIntensityZ <- scale(df$DeltaIntensity)
 # df_1$DeltaIntensityZ <- scale(df_1$DeltaIntensity)
  getFirstSignChange <- function(df, reverse = FALSE)
  {
    ### supporting function
    # sub-function to get the sign
    getSign <- function(x)
    {
      if(x < 0){
        "Negative"
      } else if(x > 0){
        "Positive"
      } else {
        0
      }
    }
    ### main logic of the function
    # start at 2; there is no change for the first one
    start <- 2
    end <- length(df)
    # if reverse flag it true, start at end and count backward
    if(isTRUE(reverse)){
      start <- length(df)
      end <- 2
    }
    # iterate though the vector...
    for(i in start:end) {
      # ...capture current and previous signs
      currentSign <- getSign(df[i])
      previousSign <- getSign(df[i-1])
      if(currentSign != previousSign) {
        # return the index just before the change
        return(i-1)
      }
    }
  }
  #get sign change indexes
  upperIndex <- getFirstSignChange(df$DeltaIntensityZ, reverse = FALSE)
  lowerIndex <- getFirstSignChange(df$DeltaIntensityZ, reverse = TRUE)
  UC<-NA
  LC<-NA
  
  UC<-as.numeric(df[upperIndex,"I"])
  LC<-as.numeric(df[lowerIndex,"I"])
  
  L1<-data.frame(matrix(rep(1,each=nrow(df)), ncol=1, byrow=TRUE))
  names(L1)<-"LineRegion"
  
  df<-data.frame(df,L1)
  # second region, in the line
  df[upperIndex+1:lowerIndex,"LineRegion"] <- 2#small bugs were fixed here
  # the second hockey stick
  df[lowerIndex+1:nrow(df),"LineRegion"] <- 3#small bugs were fixed here
  df<-df %>% dplyr::select(I,C,LineRegion)
  df<-na.omit(df)
  
  df_ <- Data0
  
  df_<-df_[order(df_$C),]
  
  L1<-data.frame(matrix(rep(1,each=nrow(df_)), ncol=1, byrow=TRUE))
  names(L1)<-"LineRegion"
  
  df_<-data.frame(df_,L1)
  
  df_$LineRegion[df_$I >= UC]<-1
  
  df_$LineRegion[df_$I < LC]<-3
  
  df_$LineRegion[df_$I < UC & df_$I >= LC]<-2
  
  df<-df_
  
  ########################
  
  #here we would like to split the data by regions to calculate linear models
  df1<-df %>% dplyr::filter(LineRegion==1)
  df2<-df %>% dplyr::filter(LineRegion==2)
  df3<-df %>% dplyr::filter(LineRegion==3)
  
  #Line region check starts here
  # #Split the data into line regions
  # if(any(!is.na(tppData))){
  #   df1<-df %>% dplyr::filter(LineRegion==1&dataset==c("F4","F5"))
  #   df2<-df %>% dplyr::filter(LineRegion==2&dataset==c("F4","F5"))
  #   df3<-df %>% dplyr::filter(LineRegion==3&dataset==c("F4","F5"))
  # }
  
  #Inital guess:define the changepoint as the last points of regions 1 and 2
  change1<-tail(df1$C,1)
  # 
  change2<-tail(df2$C,1)# 
  # 
  
  #Set mean and variance of the top plateau
  Eq1<-mean(df1$I)
  vEq1<-var(df1$I)
  
  #define the number of samples in the blank
  nblank<-nrow(df1)
  
  #calculate t statistic
  alpha = 0.05
  tstat<-qt(1-alpha,nblank-1)*sqrt(vEq1*nblank+(vEq1*nblank)/(nblank-1))
  
  #define confidence intervals for the blank
  CI_1H<-Eq1+tstat*sqrt(vEq1)
  CI_1L<-Eq1-tstat*sqrt(vEq1)
  
  #define the start of region 2 to overlay CI control
  df2$LineRegion<-ifelse(df2$I<CI_1L,2,1)
  
  
  #if the intensity value for LR 2 is below CI, keep the LineRegion
  
  #Now repeat for change 2
  #Set mean and variance of the bottom plateau control
  Eq1<-mean(df3$I)
  vEq1<-var(df3$I)
  #define the number of samples in the bottom plateau
  nblank<-nrow(df3)
  
  #calculate t statistic
  alpha = 0.05
  tstat<- qt(1-alpha,nblank-1)*sqrt(vEq1*nblank+(vEq1*nblank)/(nblank-1))
  
  
  #define confidence intervals for the blank
  CI_1H<-Eq1+tstat*sqrt(vEq1)
  
  CI_1L<-Eq1-tstat*sqrt(vEq1)
  
  
  #If Eq11 is nan set CI_1H1 as the max Intensity in region 3
  #define the start of region 2 to overlay CI
  ph1<-df2$I[df2$I %in% tail(df2$I)]
  ph1<-ifelse(ph1>=CI_1H,2,3)
  df2$LineRegion[df2$I %in% tail(df2$I)]<-ph1
  
  #if the intensity values of df2 fall above the 95% intervals
  #LR 2 is conserved, otherwise, set as LR 3
  
  df<-rbind(df1,df2,df3)
  
  #df1<-df1 %>% filter(CC>0)
  df$LineRegion<-as.numeric(df$LineRegion)
  
  #changepoint check control:
  ctest<-df %>% dplyr::group_by(C,LineRegion) %>% dplyr::summarise(n=n()) %>% dplyr::ungroup()
  ctest<-na.omit(ctest)
  if(any(ctest$n<max(ctest$n))){
    #get the data points that belong to two regions
    split<-ctest %>% dplyr::filter(n<max(ctest$n))%>% dplyr::select(C,LineRegion)
    split$LineRegion<-as.numeric(split$LineRegion)
    
    #find duplicates and subtract 1 from LineRegion
    #find duplicates and get the first one
    dap<-split
    dap<-plyr::ddply(dap,"C",summarize,LineRegion = min(LineRegion))
   # dap<-dap[1,]
    
    for(i in 1:length(unique(dap$C))){
      rows<-find_pat(dap$C[i],df$C)
      df[rows,]$LineRegion<-dap$LineRegion[i]
    }
    df<-na.omit(df)
    
    # dap<-split
    # dap<-plyr::ddply(dap,"C",summarize,LineRegion = max(LineRegion))
    # dap<-dap[2,]
    # 
    # for(i in 1:length(unique(dap$C))){
    #   rows<-find_pat(dap$C[i],df$C)
    #   df[rows,]$LineRegion<-dap$LineRegion[i]
    # }
    # df<-na.omit(df)
    
  }
  #run again to correct LineRgions on the second changepoint
  # ctest<-df %>% dplyr::group_by(C,LineRegion) %>% dplyr::summarise(n=n()) 
  # if(any(ctest$n < max(ctest$n))){
  #   #get the data points that belong to two regions
  #   split<-ctest %>% dplyr::filter(n<max(ctest$n))%>% dplyr::select(C,LineRegion)
  #   split$LineRegion<-as.numeric(split$LineRegion)
  #   
  #   #find duplicates and subtract 1 from LineRegion
  #   #find duplicates and get the first one
  #   dap<-split
  #   dap<-dap %>% dplyr::group_by(C) %>% dplyr::summarise(LineRegion=min(.$LineRegion)) %>% dplyr::ungroup()
  #   
  #   for(i in 1:length(unique(dap$C))){
  #     rows<-find_pat(dap$C[i],df$C)
  #     df[rows,]$LineRegion<-dap$LineRegion[i]
  #   }
  #   df<-na.omit(df)
  #   
  # }
   
  
}

df_<-DLR(d,Data0)
df_1<-DLR(d1,Data1)
#calculate trilinear model stats 
#DF is the null data (both conditions)
#df is the control
#df1 is the treated
tlstat<-function(DF,df,df1){

change1<- tail(df$C[df$LineRegion==1],1)
change2<- tail(df$C[df$LineRegion==2],1)
change3<- tail(df1$C[df1$LineRegion==1],1)
change4<- tail(df1$C[df1$LineRegion==2],1)

#Calculate linear models per line region
mean1<-df  %>% 
  dplyr::group_by(LineRegion) %>% 
  dplyr::do( M1 = (lm(.$I ~ .$C, . ) ) ) 
#define linear models with outputs
for (i in 1:list.count(mean1$M1)){
  mean1$formula<-formula(mean1$M1[[i]])
  mean1$slope[i]<-coef(mean1$M1[[i]])[2]
  mean1$intercept[i]<-coef(mean1$M1[[i]])[1]
  mean1$rss[i]<-deviance(mean1$M1[[i]])
  mean1$CI[i]<-(confint(mean1$M1[[i]])[1]+confint(mean1$M1[[i]])[2])/2
  mean1$dataset<-"control"
  
}
mean1_1<-df1  %>% 
  dplyr::group_by(LineRegion) %>% 
  dplyr::do( M1 = (lm(.$I ~ .$C, . ) ) ) 
#define linear models with outputs
for (i in 1:list.count(mean1_1$M1)){
  mean1_1$formula<-formula(mean1_1$M1[[i]])
  mean1_1$slope[i]<-coef(mean1_1$M1[[i]])[2]
  mean1_1$intercept[i]<-coef(mean1_1$M1[[i]])[1]
  mean1_1$rss[i]<-deviance(mean1_1$M1[[i]])
  mean1_1$CI[i]<-(confint(mean1_1$M1[[i]])[1]+confint(mean1_1$M1[[i]])[2])/2
  mean1_1$dataset<-"treated"
}

#call data for the null hypothesis
DF$C<-as.numeric(DF$C)
Mup<-mean(c(change1,change3))#change this
Mdo<- mean(c(change2,change4))#change this        
DF1<-DF %>% base::subset(DF$C<= Mup) %>% dplyr::mutate(LineRegion=1)
DF2<-DF %>% base::subset(DF$C<Mdo & DF$C > Mup) %>% dplyr::mutate(LineRegion=2)
DF3<-DF %>% base::subset(DF$C>= Mdo) %>% dplyr::mutate(LineRegion=3)
DF<-rbind(DF1,DF2,DF3)
mean3<-DF  %>% 
  dplyr::group_by(LineRegion) %>% 
  dplyr::do( M1 = (lm(.$I ~ .$C, . ) ) ) 
#define nls linear models with outputs
for (i in 1:list.count(mean3$M1)){
  mean3$formula<-formula(mean3$M1[[i]])
  mean3$slope[i]<-coef(mean3$M1[[i]])[2]
  mean3$intercept[i]<-coef(mean3$M1[[i]])[1]
  mean3$rss[i]<-deviance(mean3$M1[[i]])
  mean3$CI[i]<-(confint(mean3$M1[[i]])[1]+confint(mean3$M1[[i]])[2])/2
  mean3$dataset<-"null"
  
}


results<-list(mean1,mean1_1,mean3)
}

tlresults<-tlstat(DF,df_,df_1)

#plot confidence intervals
#tl results are the model fit results
#DF is the null data 
#df_ is the control 
#df_1 is the treated 
tlCI<-function(tlresults,DF,df_,df_1,overlay=TRUE){

  pred1<-predict(tlresults[[3]]$M1[[1]], interval="confidence")
  pred2<-predict(tlresults[[3]]$M1[[2]], interval="confidence")
  pred3<-predict(tlresults[[3]]$M1[[3]], interval="confidence")
  Pred<-NA
  if (any(!is.na(pred3)) & any(!is.na(pred2))){
    FIT<-list(pred1[,1],pred2[,1],pred3[,1])
    LOW<-list(pred1[,2],pred2[,2],pred3[,2])
    HI<-list(pred1[,3],pred2[,3],pred3[,3])
  }else if(any(is.na(pred2)) & any(!is.na(pred3))){
    FIT<-list(pred1[,1],pred3[,1])
    LOW<-list(pred1[,2],pred3[,2])
    HI<-list(pred1[,3],pred3[,3])
  }else{
    FIT<-list(pred1[,1],pred2[,1])
    LOW<-list(pred1[,2],pred2[,2])
    HI<-list(pred1[,3],pred2[,3])
  }
  Pred<-data.frame(unlist(FIT),unlist(LOW),unlist(HI))
  Pred<-Pred[1:length(DF$C),]
  Pred<-data.frame(Pred,DF$C,DF$I)
  names(Pred)<-c("fit","lower","upper","C","I")
  Pred$Treatment<-tlresults[[3]]$dataset[1]
  
  PLN<-ggplot(Pred, aes(x = C,y = fit,color=Treatment)) +geom_point(aes(x=C,y=I))+geom_ribbon(data=Pred,aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ xlab("Temperature (\u00B0C)")+ylab("Relative Intensity")+ annotate("text", x=62, y=1, label= paste("RSS= ",round(sum(tlresults[[3]]$rss),3)))
  
  pred1<-predict(tlresults[[2]]$M1[[1]], interval="confidence")
  pred2<-predict(tlresults[[2]]$M1[[2]], interval="confidence")
  pred3<-predict(tlresults[[2]]$M1[[3]], interval="confidence")
  Pred1<-NA
  if (any(!is.na(pred3)) & any(!is.na(pred2))){
    FIT<-list(pred1[,1],pred2[,1],pred3[,1])
    LOW<-list(pred1[,2],pred2[,2],pred3[,2])
    HI<-list(pred1[,3],pred2[,3],pred3[,3])
  }else if(any(is.na(pred2)) & any(!is.na(pred3))){
    FIT<-list(pred1[,1],pred3[,1])
    LOW<-list(pred1[,2],pred3[,2])
    HI<-list(pred1[,3],pred3[,3])
  }else{
    FIT<-list(pred1[,1],pred2[,1])
    LOW<-list(pred1[,2],pred2[,2])
    HI<-list(pred1[,3],pred2[,3])
  }
  Pred1<-data.frame(unlist(FIT),unlist(LOW),unlist(HI))
  Pred1<-Pred1[1:length(df_1$C),]
  Pred1<-data.frame(Pred1,df_1$C,df_1$I)
  names(Pred1)<-c("fit","lower","upper","C","I")
  Pred1$Treatment<-tlresults[[2]]$dataset[1]
  
  PLR_P1<-ggplot(Pred1, aes(x = C,y = fit,color=Treatment))+geom_point(Pred1, mapping=aes(x = C,y = I,color=Treatment)) +geom_ribbon(data=Pred1,aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)
  
  pred1<-predict(tlresults[[1]]$M1[[1]], interval="confidence")
  pred2<-predict(tlresults[[1]]$M1[[2]], interval="confidence")
  pred3<-predict(tlresults[[1]]$M1[[3]], interval="confidence")
  Pred2<-NA
  if (any(!is.na(pred3)) & any(!is.na(pred2))){
    FIT<-list(pred1[,1],pred2[,1],pred3[,1])
    LOW<-list(pred1[,2],pred2[,2],pred3[,2])
    HI<-list(pred1[,3],pred2[,3],pred3[,3])
  }else if(any(is.na(pred2)) & any(!is.na(pred3))){
    FIT<-list(pred1[,1],pred3[,1])
    LOW<-list(pred1[,2],pred3[,2])
    HI<-list(pred1[,3],pred3[,3])
  }else{
    FIT<-list(pred1[,1],pred2[,1])
    LOW<-list(pred1[,2],pred2[,2])
    HI<-list(pred1[,3],pred2[,3])
  }
  Pred2<-data.frame(unlist(FIT),unlist(LOW),unlist(HI))
  Pred2<-Pred2[1:length(df_$C),]
  Pred2<-data.frame(Pred2,df_$C,df_$I)
  names(Pred2)<-c("fit","lower","upper","C","I")
  Pred2$Treatment<-tlresults[[1]]$dataset[1]
  #Area under the curve using trapezoid rule
  P1_AUC <- pracma::trapz(Pred1$C,Pred1$I)
  P2_AUC <- pracma::trapz(Pred2$C,Pred2$I)
 
  
  PLR_P2<-PLR_P1+geom_point(Pred2, mapping=aes(x = C,y = I,color=Treatment)) +geom_ribbon(data=Pred2,aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
    xlab("Temperature (\u00B0C)")+ylab("Relative Intensity")
  if(overlay=="TRUE"){
    AUCd<-round(P1_AUC-P2_AUC,2)
    Tm1_tr<-Pred1[which.min(abs(Pred1$fit - 0.5)),'C']
    Tm2_c<-Pred2[which.min(abs(Pred2$fit - 0.5)),'C']
    Tm_d<-Tm1_tr-Tm2_c
    p<-expression(paste(Delta, "AUCdiff"))
    if(AUCd>0){
      P1_AUC <- pracma::trapz(Pred1$C,Pred1$lower)
      P2_AUC <- pracma::trapz(Pred2$C,Pred2$upper)
      AUCd<-round(P1_AUC-P2_AUC,2)
    }else{
      P1_AUC <- pracma::trapz(Pred1$C,Pred1$upper)
      P2_AUC <- pracma::trapz(Pred2$C,Pred2$lower)
      AUCd<-round(P1_AUC-P2_AUC,2)
    }
    AUCd<-as.numeric(AUCd)
    PLR_P2<-PLR_P1+geom_point(Pred2, mapping=aes(x = C,y = I,color=Treatment)) +geom_ribbon(data=Pred2,aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
      xlab("Temperature (\u00B0C)")+ylab("Relative Intensity")+ annotate("text", x=62, y=1, label= paste("\u03A3","RSS= ",round(sum(sum(tlresults[[2]]$rss),sum(tlresults[[1]]$rss)),3)))+
      annotate("text", x=62, y=0.95, label=  paste("\u0394", "AUC = ",AUCd))+ annotate("text", x=62, y=0.9, label= paste("\u0394","Tm = ",Tm_d,"\u00B0C"))
    #bquote(Value~is~sigma~R^{2}==.(r2.value)))
  PLR_P2<-grid.arrange(PLN,PLR_P2, ncol=2)
  print(PLR_P2)
  }else if(overlay=="FALSE"){
   PLR<-PLR_P2+geom_point(data=Pred,mapping=aes(x=C,y=I))+geom_ribbon(data=Pred,aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+facet_wrap("Treatment")  
   print(PLR)
  }
}

plotTL<-tlCI(tlresults,DF,df_,df_1,overlay=TRUE)

plotTL<-tlCI(tlresults,DF,df_,df_1,overlay=FALSE)


