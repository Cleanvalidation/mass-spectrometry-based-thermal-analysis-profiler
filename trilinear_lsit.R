library(tidyverse)
library(minpack.lm)
library(rlist)
library(data.table)
library(knitr)
library(ggthemes)
library(gridExtra)
library(grid)
library(readxl)
library(nls2)
library(minpack.lm)
library(nlstools)
library(pkgcond)
library(rlist)
library(pracma)

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
    rlist::list.group(Accession,sample) %>%
    dplyr::mutate(value = value / value[temperature == min(temperature)]) %>%dplyr::ungroup()
  
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
    rlist::list.group(Accession, sample) %>%
    dplyr::filter(n() == 10) %>%  ## all 10 temperatures are present
    dplyr::summarise(T7 = value[temperature == temperatures[[7]]] / value[temperature == temperatures[[1]]],
                     T9 = value[temperature == temperatures[[9]]] / value[temperature == temperatures[[1]]],
                     T10 = value[temperature == temperatures[[10]]] / value[temperature == temperatures[[1]]]) %>% filter(T7 >= 0.4 & T7 <= 0.6 & T9 < 0.3 & T10 < 0.2)%>%dplyr::ungroup()
  
  ## split[[i]] by sample group and filter
  l.bytype <- split[[i]].data.frame(df.jointP, df.jointP$sample)
  
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
  d<-length(df.fit$fit)
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
    rlist::list.group(sample, Accession) %>%
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
  
  df <- df %>%dplyr::ungroup()
  
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
        rlist::list.group(t) %>%
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

d$LineRegion<-NA


DF<-d %>% base::split.data.frame(.$uniqueID,.$dataset) 
d_<-d %>% dplyr::filter(CC == 0) %>%  base::split.data.frame(.$uniqueID,.$dataset) 
d_1<-d %>% dplyr::filter(CC > 0) %>%  base::split.data.frame(.$uniqueID,.$dataset) 

# #prepare a list of datasets
# datalist<-unique(d$dataset)
# #prepare a list of concentration values
# conclist<-unique(d$CC)
# #prepare a list of proteins

###selection of sample compound
d<-d %>%
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
# DF<-d
# DF<-DF[order(DF$C),]
#Alternative hypothesis

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
    BS<-dplyr::sample_n(Data0,n,replace=TRUE)
    #  BS1<-sample_n(data.frame(Data1),n,replace=TRUE)
    #generate mean intensities per $C value
    BA<-BS %>%dplyr::group_by(C) %>% dplyr::summarise(n=n(),mean=mean(I,na.rm=TRUE),var=var(I,na.rm=TRUE))
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

DLR<-function(d){
  
  #preallocate final result as a list
  df_n<-vector(mode = "list", length(d))
  df1<-df_n
  df2<-df1
  df3<-df1
  df_1<-df_n
  df0<-df_n
  df_0<-df_n
  
  df_1<-lapply(d, function(x) {x %>% dplyr::group_by(C) %>% dplyr::summarise(.,
                                                                            I=mean(I,na.rm=TRUE))
    
  })
  #Add ranks to the list
  df_1<-lapply(df_1, function(x) {
    dplyr::mutate(x,LineRegion=dplyr::ntile(desc(I),3))
  })

  d<-lapply(d,function(x) x %>% arrange(C))
  for(i in 1:length(d)){
    
    #assign line regions on original dataset
    #here we would like to split[[i]] the data by regions to calculate linear models
    df1[[i]]<-d[[i]] %>% dplyr::filter(I>=tail(df_1[[i]]$I[df_1[[i]]$LineRegion==1],1)) %>% dplyr::mutate(LineRegion=1) 
    df2[[i]]<-d[[i]] %>% dplyr::filter(I<=tail(df_1[[i]]$I[df_1[[i]]$LineRegion==1],1),I>=tail(df_1[[i]]$I[df_1[[i]]$LineRegion==2],1)) %>% dplyr::mutate(LineRegion=2)
    df3[[i]]<-d[[i]]%>%dplyr::filter(I<=head(df_1[[i]]$I[df_1[[i]]$LineRegion==3],1)) %>% dplyr::mutate(LineRegion=3)
 
  
  
    #Inital guess:define the changepoint as the last points of regions 1 and 2
    change1<-tail(df1[[i]]$C,1)
    # 
    change2<-tail(df2[[i]]$C,1)# 
    # 
    
    #Set mean and variance of the top plateau
    Eq1<-mean(df1[[i]]$I)
    vEq1<-var(df1[[i]]$I)
    
    #define the number of samples in the blank
    nblank<-nrow(df1[[i]])
    
    #calculate t statistic
    alpha <- 0.05
    tstat<-numeric(1)
    
    tstat<-qt(1-alpha,nblank-1)*sqrt(vEq1*nblank+(vEq1*nblank)/(nblank-1))
    
    #define confidence intervals for the blank
    CI_1H<-Eq1+tstat*sqrt(vEq1)
    CI_1L<-Eq1-tstat*sqrt(vEq1)
    
    #define the start of region 2 to overlay CI control
    df2[[i]]$LineRegion<-ifelse(df2[[i]]$I<CI_1L,2,1)
    
    
    #if the intensity value for LR 2 is below CI, keep the LineRegion
    
    #Now repeat for change 2
    #Set mean and variance of the bottom plateau control
    Eq1<-mean(df3[[i]]$I)
    vEq1<-var(df3[[i]]$I)
    #define the number of samples in the bottom plateau
    nblank<-nrow(df3[[i]])
    
    #calculate t statistic
    alpha <- 0.05
    tstat<- qt(1-alpha,nblank-1)*sqrt(vEq1*nblank+(vEq1*nblank)/(nblank-1))
    
    
    #define confidence intervals for the blank
    CI_1H<-Eq1+tstat*sqrt(vEq1)
    
    CI_1L<-Eq1-tstat*sqrt(vEq1)
    
    #If Eq11 is nan set CI_1H1 as the max Intensity in region 3
    #define the start of region 2 to overlay CI
    ph1<-df2[[i]]$I[df2[[i]]$I %in% tail(df2[[i]]$I)]
    ph1<-ifelse(ph1>=CI_1H,2,3)
    df2[[i]]$LineRegion[df2[[i]]$I %in% tail(df2[[i]]$I)]<-ph1
    
    #if the intensity values of df2 fall above the 95% intervals
    #LR 2 is conserved, otherwise, set as LR 3
    
    df_0[[i]]<-rbind(df1[[i]],df2[[i]],df3[[i]])
    
    df_0[[i]]<-df_0[[i]][order(df_0[[i]]$C),]
    #df1<-df1 %>% filter(CC>0)
    df_0[[i]]$LineRegion<-as.factor(df_0[[i]]$LineRegion)
    
    df_n[[i]]<-df_0[[i]]
  }
  return(df_n)
}

CP<-function(df_0){
  
  df_n<-vector(mode = "list", length(df_0))
  ctest<-df_n
  dap<-df_n
  Split<-df_n
  #changepoint check control:
  for(i in 1:length(df_0)){
    df_0[[i]]<-df_0[[i]] %>% dplyr::arrange(C)
    ctest<-NA
    ctest<-df_0[[i]] %>%dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup()
    
    ###############################################
    
    #get the data points that belong to two regions
    Split<-ctest %>%subset(n<(max(n))) %>% data.frame()
    Split<-Split %>% subset(Split$C[duplicated(Split$C)] %in% Split$C)
    Split$LineRegion<-as.numeric(Split$LineRegion)
    
    mode<-function(x){
      ux<-unique(x)
      ux[which.max(tabulate(match(x,ux)))]
    }
    #find duplicates and subtract 1 from LineRegion
    #find duplicates and get the first one
    dap<-list()
    dap[[i]]<-df_0[[i]] %>% subset(df_0[[i]]$C %in% Split$C) 
    dap[[i]]$LineRegion<-as.numeric(dap[[i]]$LineRegion)
    ##dap[[i]]<-dap[[i]] %>% dplyr::group_by(C) %>% dplyr::mutate(LineRegion=ifelse(LineRegion<mode(LineRegion),mode(LineRegion),LineRegion))
    
    
    dap[[i]]<-dap[[i]] %>%dplyr::group_by(C) %>%  dplyr::mutate(LineRegion=ifelse(as.numeric(LineRegion)>mode(as.numeric(LineRegion)),ifelse(LineRegion==2 & I<0.5 & max(as.numeric(LineRegion==3)),mode(as.numeric(LineRegion))+1,mode(as.numeric(LineRegion))),as.numeric(LineRegion)))
    df_0[[i]]$LineRegion[df_0[[i]]$C %in% dap[[i]]$C]<-dap[[i]]$LineRegion
    df_0[[i]]$LineRegion<-as.factor(df_0[[i]]$LineRegion)
   
    #check again
    # df_0[[i]]<-df_0[[i]] %>% dplyr::arrange(C)
    # ctest<-NA
    # ctest<-df_0[[i]] %>%dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup()
    # 
    # ###############################################
    # 
    # #get the data points that belong to two regions
    # Split<-ctest %>%subset(n<(max(n))) %>% data.frame()
    # Split<-Split %>% subset(Split$C[duplicated(Split$C)] %in% Split$C)
    # Split$LineRegion<-as.numeric(Split$LineRegion)
    # 
    # #find duplicates and subtract 1 from LineRegion
    # #find duplicates and get the first one
    # dap<-list()
    # dap[[i]]<-df_0[[i]] %>% subset(df_0[[i]]$C %in% Split$C) 
    # dap[[i]]$LineRegion<-as.numeric(dap[[i]]$LineRegion)
    # #dap[[i]]<-dap[[i]] %>% dplyr::group_by(C) %>% dplyr::mutate(LineRegion=ifelse(min(LineRegion)<max(LineRegion),
    # dap[[i]]<-dap[[i]] %>% dplyr::group_by(C) %>% dplyr::mutate(LineRegion=ifelse(min(LineRegion)<max(LineRegion),
    #                                                                               ifelse(LineRegion<max(LineRegion) & I<0.5,LineRegion+1,LineRegion-1),LineRegion))
    # 
    # df_0[[i]]$LineRegion[df_0[[i]]$C %in% dap[[i]]$C]<-dap[[i]]$LineRegion
    # df_0[[i]]$LineRegion<-as.factor(df_0[[i]]$LineRegion)
    # 
    # return(df_0)
  }
 
  return(df_0)
}


  
  #preallocate list
results<-vector(mode = "list", length(d_))
results_t<-vector(mode = "list",length(d_1))
results_n<-vector(mode = "list",length(DF))

results<-DLR(d_)#First guess at line regions
results_t<-DLR(d_1)
results_n<-DLR(DF)

RS<-vector(mode = "list", length(d_))
RST<-vector(mode = "list", length(d_1))
RSN<-vector(mode = "list", length(DF))

RS<-CP(results)#reassign changepoints shared between line regions
RST<-CP(results_t)#reassign changepoints shared between line regions
RSN<-CP(results_n)#reassign changepoints shared between line regions

saveRDS(RS,'RS.Rds')
saveRDS(RST,'RST.Rds')
saveRDS(RSN,'RSN.Rds')

df_<-RS
df_1<-RST
DF<-RSN

# df_<-readRDS("~/test_01/Canonical/CS7290-/internal data/RS.Rds")
# df_1<-readRDS("~/test_01/Canonical/CS7290-/internal data/RST.Rds")
# DF<-readRDS("~/test_01/Canonical/CS7290-/internal data/RSN.Rds")


tlstat<-function(DF,df,df1){
  
  
  #preallocate list results
  #control/vehicle
  i=1
  
  mean1<-list()
  mean1[[i]]<-data.frame(LineRegion=factor(1),formula=list(1),slope=numeric(1),intercept=numeric(1),rss=numeric(1),Rsq=numeric(1),dataset=factor(1),uniqueID=character(1))

  #treated
 
  mean1_1<-mean1
  #null
  mean3<-mean1
  #Calculate linear models per line region
  for(i in 1:length(df)){
   
   mean1[[i]]<-df[[i]] %>% 
      dplyr::group_by(LineRegion) %>% 
      dplyr::do(.,M1 = (lm(.$I ~ .$C, . ) ) ) %>% 
      dplyr::mutate(.,slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),dataset="vehicle") 
   mean1[[i]]$LineRegion<-seq_along(mean1[[i]]$M1) 
    
    
  }
  
  #define linear models with outputs
  for (i in 1:length(mean1)){
    
   mean1[[i]]<- mean1[[i]] %>% dplyr::rowwise(.)%>% dplyr::mutate(slope=as.numeric(coef(M1)[2]),
                                                                             intercept=as.numeric(coef(M1)[1]),
                                                       rss=deviance(M1),Rsq=summary(M1)$r.squared, 
      dataset=stringr::str_c(c(df[[i]]$dataset[1],"vehicle"),collapse=" "),uniqueID=df[[i]]$uniqueID[1])
  
  }
  
  #Treated
 
  for(i in 1:length(df1)){
   
    mean1_1[[i]]<-df1[[i]] %>% 
      dplyr::group_by(LineRegion) %>% 
      dplyr::do(.,M1 = (lm(.$I ~ .$C, . ) ) ) %>% 
      dplyr::mutate(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),dataset="treated") 
    mean1_1[[i]]$LineRegion<-seq_along(mean1_1[[i]]$M1)
    
   
  }
  
  #define linear models with outputs
  for (i in 1:length(mean1_1)){
 mean1_1[[i]]<-mean1_1[[i]]%>% dplyr::rowwise(.)%>% dplyr::mutate(slope=as.numeric(coef(M1)[2]),
                                                             intercept=as.numeric(coef(M1)[1]),
                                                             rss=deviance(M1),Rsq=summary(M1)$r.squared, 
                                                             dataset=stringr::str_c(c(df1[[i]]$dataset[1],"treated"),collapse=" "),uniqueID=df1[[i]]$uniqueID[1])
 
  }
  
  
  # null hypothesis
  for(i in 1:length(DF)){

 mean3[[i]]<- DF[[i]] %>% 
        dplyr::group_by(LineRegion) %>% 
        dplyr::do(.,M1 = (lm(.$I ~ .$C, . ) ) ) %>% 
        dplyr::mutate(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),dataset="null") 
 mean3[[i]]$LineRegion<-seq_along(mean3[[i]]$M1)
  }
    
    
    #define nls linear models with outputs
    for (i in 1:length(mean3)){
mean3[[i]]<-mean3[[i]]%>% dplyr::rowwise(.)%>% dplyr::mutate(slope=as.numeric(coef(M1)[2]),
                                                         intercept=as.numeric(coef(M1)[1]),
                                                         rss=deviance(M1),Rsq=summary(M1)$r.squared, 
                                                         dataset=stringr::str_c(c(DF[[i]]$dataset[1],"null"),collapse=" "),uniqueID=DF[[i]]$uniqueID[1])
      
    }
   for(i in 1:length(df)){
   results[[i]]<-data.frame()
   results[[i]]<-rbind(mean1[[i]],mean1_1[[i]],mean3[[i]])
                      
  
   }
   results
  }


tlresults<-list()
tlresults<-tlstat(DF,df_,df_1)
tlresults1<-tlresults#save unfiltered data
#apply filters prior to hypothesis testing
#tlresults<-tlresults %>% keep(function(x) min(as.numeric(x$Rsq),na.rm=TRUE) > 0.8)#the linear region have the largest slope < 0.03
tlresults<-tlresults %>% keep(function(x) sum(x$rss[stringr::str_detect(tolower(x$dataset), pattern = "null")],na.rm=TRUE) < 3)#remove data with extremely large RSS values 

tlresults<-tlresults %>% keep(function(x) sum(x$rss[!stringr::str_detect(tolower(x$dataset), pattern = "null")],na.rm=TRUE) < 1)
tlresults<-tlresults %>% keep(function(x) sum(x$rss[stringr::str_detect(tolower(x$dataset), pattern = "null")]) > sum(x$rss[!stringr::str_detect(tolower(x$dataset), pattern = "null")]))#remove data with extremely large RSS values 
tlresults<-tlresults %>% keep(function(x) max(x$slope[x$LineRegion==2],na.rm=TRUE) < -0.04)#the linear region have the largest slope < 0.03
#tlresults<-tlresults %>% keep(function(x) length(x$slope)>5)#remove list values with less than 5 rows
#tlresults<-tlresults %>% keep(function(x) abs(max(x$slope[!x$LineRegion==2] ,na.rm=TRUE)) < 0.1)#keeps plateau values where the min abs(slope) < 0.06
#steepest slope in vehicle and treatment has to be less than 0.06C


Nsum<-lapply(tlresults, function(x) x %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "null")) %>% 
               dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS) %>% data.table::first(.$RSS)) 

#get the summed rss values for vehicle
Rssv<-lapply(tlresults, function(x) x %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% 
               dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS) %>% data.table::first(.$RSS)) 
#get the summed rss values for treated
Rsst<-lapply(tlresults, function(x) x %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% 
              dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS) %>% data.table::first(.$RSS)) 
#find the rss difference between treated and vehicle 

Dsum<- map2(Rssv,Rsst,~ .y-.x) %>% rbindlist(.) %>% data.frame() %>% dplyr::mutate(rank=dplyr::ntile(.,7)) 
Dsum<-Dsum %>% dplyr::mutate(id=rownames(Dsum))
names(Dsum)[1]<-"RSSd"

#keep data where the difference in RSS is less than the null
#nsum converted to data frame

Nsum<-rbindlist(Nsum) 
Nsum<-Nsum%>% data.frame(.) %>% dplyr::mutate(id=rownames(Nsum),rank=dplyr::ntile(.,7))
names(Nsum)[1]<-"RSSn"

#mutate data frame
Dsum<-Dsum %>% dplyr::inner_join(Nsum,by = c("id" = "id"))
names(Dsum)<-c("RSSd","rank","RSSn","n_rank")
test<-data.frame()
test<-Dsum[which(Dsum$rank==7 & Dsum$RSSn>Dsum$RSSd),] %>% data.frame()#get the stable proteins (+ = Rsstreated-Rssvehicle)

rssdec<-data.frame()
rssdec<-data.frame(data.table::fsort(test$RSSd,decreasing=TRUE))#decreasing Rss differences
names(rssdec)<-"Rssd"
ids<-test %>% dplyr::mutate(id=rownames(test))

orows<-data.frame()
orows <- rssdec %>%
  dplyr::select(Rssd) %>%
  dplyr::inner_join(ids, by = c('Rssd'='RSSd'))

orows$id<-sapply(orows$id, function(x) as.numeric(as.character(x)))
  #order by RSS differences while keeping original rownames for index
#create an external data frame for stabilized proteins

#get the rows that follow the conditional

#find values for stabilized proteins and rank them


Df1<-tlresults[orows$id] #divide 1=highly destabilized,4=noeffect,7=highly stabilized
df1<-list()
#get uniqueID and dataset for stable proteins with decreasing RSS differences
df1<-lapply(Df1,function(x) x %>% dplyr::select(uniqueID,dataset) %>% head(.,1))
df1<-data.frame(rbindlist(df1))
df1$dataset<-tidyr::separate(df1,dataset," ")
df1$dataset<-df1$dataset[,2]
#unlist to data.frame
#order the original data by RSS differences
#
df2<-df_ %>% keep(function(x) (head(x$uniqueID,1) %in% df1$uniqueID & head(x$dataset,1) %in% df1$dataset)) 
df2<-data.frame(rbindlist(df2))

#inner join to get all replicates in order of stability
df_final <- df1 %>%
  dplyr::select(dataset,uniqueID) %>%
  dplyr::inner_join(df2, by = c('uniqueID' = 'uniqueID', 'dataset'='dataset'))


df_final$uniqueID<-as.factor(df_final$uniqueID)
df_final$dataset<-as.factor(df_final$dataset)


#plot confidence intervals


tlCI<-function(i,Df1,DF,df_,df_1,overlay=TRUE){
  #null data: third list subset ==3
  DF1<-rbindlist(DF)
  null<-data.frame()
  
  i<-i
  df_<-df_
  df_1<-df_1
  DF1<-DF1 %>% subset(uniqueID %in% df_final$uniqueID[i])
  
  null<-Df1[[i]] %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "null"))
  #first subset list is the element chosen from the list of proteins 
  pred1<-predict(null$M1[[1]], interval="confidence")
  pred2<-predict(null$M1[[2]], interval="confidence")
  pred3<-predict(null$M1[[3]], interval="confidence")
  
  
  
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  
 
  Pred<-NA
  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred<-rbind(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred<-rbind(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
  Pred<-rbind(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred<-rbind(pred1,pred3)
  }
  rownames(Pred)<-1:nrow(Pred)

  #Pred<-Pred[1:length(DF1$C),]##############
  Pred<-data.frame(Pred,DF1$C[1:nrow(Pred)],DF1$I[1:nrow(Pred)])################
  names(Pred)<-c("fit","lower","upper","C","I")
  
  Pred$Treatment<-null$dataset[1]##################
  Pred<-na.omit(Pred)
  rownames(Pred)<-1:nrow(Pred)
  PLN<-ggplot2::ggplot(Pred, ggplot2::aes(x = C,y = fit,color=Treatment)) +ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1],"null"))+ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::annotate("text", x=62, y=1, label= paste("RSS= ",round(sum(null$rss),3)))
  
  
  DF_f<-rbindlist(df_)%>% subset(uniqueID %in% df_final$uniqueID[i] & CC==0) %>% as.data.frame(.)
  
  vehicle<-data.frame()
  vehicle<-Df1[[i]] %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "vehicle"))
  
  pred1<-predict(vehicle$M1[[1]], interval="confidence")
  pred2<-predict(vehicle$M1[[2]], interval="confidence")
  pred3<-predict(vehicle$M1[[3]], interval="confidence")
  Pred1<-NA
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)

  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred1<-rbind(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred1<-rbind(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
    Pred1<-rbind(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred1<-rbind(pred1,pred3)
}
  
  #Pred<-Pred[1:length(DF1$C),]##############
  Pred1<-data.frame(Pred1,DF_f$C[1:nrow(Pred1)],DF_f$I[1:nrow(Pred1)])################
  names(Pred1)<-c("fit","lower","upper","C","I")
  
  Pred1$Treatment<-vehicle$dataset[1]##################
  Pred1<-na.omit(Pred1)
  rownames(Pred1)<-1:nrow(Pred1)

  
  PLR_P1<-ggplot2::ggplot(Pred1, ggplot2::aes(x = C,y = fit,color=Treatment))+ggplot2::geom_point(Pred1, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +ggplot2::geom_ribbon(data=Pred1,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)
  
  DF_f1<-data.frame()
  DF_f1<-rbindlist(df_1) 
  DF_f1<-DF_f1%>% subset(uniqueID %in% df_final$uniqueID[i] & CC>0) 
  
  treated<-data.frame()
  treated<-Df1[[i]] %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "treated"))
  
  pred1<-predict(treated$M1[[1]], interval="confidence")
  pred2<-predict(treated$M1[[2]], interval="confidence")
  pred3<-predict(treated$M1[[3]], interval="confidence")

  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  Pred2<-NA
  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred2<-rbind(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred2<-rbind(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
    Pred2<-rbind(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred2<-rbind(pred1,pred3)
  }
  rownames(Pred2)<-1:nrow(Pred2)
  
  #Pred<-Pred[1:length(DF1$C),]##############
  Pred2<-data.frame(Pred2,DF_f1$C[1:nrow(Pred2)],DF_f1$I[1:nrow(Pred2)])################
  names(Pred2)<-c("fit","lower","upper","C","I")
  
  Pred2$Treatment<-treated$dataset[1]##################
  Pred2<-na.omit(Pred2)
  rownames(Pred2)<-1:nrow(Pred2)
  #Area under the curve using trapezoid rule
  P1_AUC <- pracma::trapz(Pred1$C,Pred1$I)
  P2_AUC <- pracma::trapz(Pred2$C,Pred2$I)
  
  
  PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")
  if(overlay=="TRUE"){
    AUCd<-round(P2_AUC-P1_AUC,2)
    Tm1<-data.frame()
    Tm2<-data.frame()
    
    
    Tm1<-Pred1[which.min(abs(Pred1$fit - 0.5)),'C']#pred1 is vehicle
    Tm2<-Pred2[which.min(abs(Pred2$fit - 0.5)),'C']#pred2 is treated
    Tm_d<-Tm2 -Tm1
    p<-expression(paste(Delta, "AUCdiff"))
    if(AUCd>0){
      P1_AUC <- pracma::trapz(Pred1$C,Pred1$lower)
      P2_AUC <- pracma::trapz(Pred2$C,Pred2$upper)
      AUCd<-round(P2_AUC-P1_AUC,2)
    }else{
      P1_AUC <- pracma::trapz(Pred1$C,Pred1$upper)
      P2_AUC <- pracma::trapz(Pred2$C,Pred2$lower)
      AUCd<-round(P2_AUC-P1_AUC,2)
    }
    AUCd<-as.numeric(AUCd)
    PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1],"alternative"))+
      ggplot2::annotate("text", x=62, y=1, label= paste("\u03A3","RSS= ",round(sum(Df1[[i]] %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% dplyr::select(rss) %>% sum(.),Df1[[i]] %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% dplyr::select(rss) %>% sum(.)),3)))+
      ggplot2::annotate("text", x=62, y=0.9, label=  paste("\u0394", "AUC = ",AUCd))+ ggplot2::annotate("text", x=62, y=0.8, label= paste("\u0394","Tm = ",Tm_d,"\u00B0C"))
    #bquote(Value~is~sigma~R^{2}==.(r2.value)))
    PLR_P2<-grid.arrange(PLN,PLR_P2, ncol=2)
    print(PLR_P2)
  }else if(overlay=="FALSE"){
    PLR<-PLR_P2+ggplot2::geom_point(data=Pred,mapping=ggplot2::aes(x=C,y=I))+ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1],"alternative"))+facet_wrap("Treatment") 
    print(PLR)
  }
}
i=1
plotTL<-tlCI(i,Df1,DF,df_,df_1,overlay=TRUE)

plotTL<-tlCI(i,Df1,DF,df_,df_1,overlay=FALSE)


