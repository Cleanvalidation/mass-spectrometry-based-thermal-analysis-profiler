
library(minpack.lm)
library(rlist)
library(data.table)
library(knitr)
library(ggthemes)
library(gridExtra)
library(grid)
library(readxl)
library(nls2)
library(stats)
library(ggplot2)
library(pkgcond)
library(rlist)
library(pracma)
library(fs)
library(tidyverse)
library(splines)
library(mgcv)
library(purrr)
library(nlstools)
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
#' @importFrom stringr stringr::str_extract
#' @export
#' 
read_cetsa <- function(f){
  df.raw <- readxl::read_excel(f)
  df <- df.raw%>%
    dplyr::select(Accession, tidyselect::starts_with('Abundance')) %>%
    tidyr::gather('id', 'value', -Accession) %>%
    dplyr::mutate(sample = stringr::str_extract(id, '(?<=Abundance: ).*(?=:)')) %>%
    dplyr::mutate(temp_ref = stringr::str_extract(id, '(?<=: )([0-9]+[N|C]?)'))
  df
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
    dplyr::mutate(value = value / value[temperature == min(temperature)]) %>% unique(.)
  
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
  df$Accession<-as.factor(df$Accession)
  df$temperature<-as.factor(df$temperature)
  df.jointP <- df %>%
    dplyr::group_by(Accession,sample) %>% dplyr::mutate(n=dplyr::n()) %>% #copies the number of temperatures per group
    dplyr::filter(n>=10) %>% #removes groups with less than 10 temperature channels
    dplyr::group_split(.) #split into groups
  
  df.jointP<-lapply(df.jointP,function(x) x %>% dplyr::mutate(.,T7 = value[temperature == temperatures[7]]/value[temperature == temperatures[1]],
                                                              T9 = value[temperature == temperatures[9]]/value[temperature == temperatures[1]],
                                                              T10 = value[temperature == temperatures[10]]/ value[temperature == temperatures[1]]) %>% 
                      dplyr::filter(T7 >= 0.4 & T7 <= 0.6 & T9 < 0.3 & T10 < 0.2))#normalization from TPP
  
  #convert to df # dplyr::bind_rows #fix
  df.jointP<-data.table::rbindlist(df.jointP)
  ## split[[i]] by sample group and filter
  l.bytype <- split.data.frame(df.jointP, df.jointP$sample)
  
  ## determine which group contains the greatest number of curves and use this for normalization
  n.filter <- lapply(l.bytype, nrow)
  df.normP <- l.bytype[[which.max(n.filter)]]
  norm.accessions <- df.normP$Accession
  
  ## calculate median for each sample group
  
  df.mynormset <- df %>% base::subset(Accession %in% norm.accessions)
  
  df.median <- df %>%
    dplyr::group_by(sample,temperature) %>%
    dplyr::summarise(value = median(value))
  
  df.median$temperature<-as.numeric(levels(df.median$temperature))[df.median$temperature]
  ## fit curves to the median data
  df.fit <- df.median %>%
    dplyr::group_by(sample) %>% 
    dplyr::do(fit = cetsa_fit(d = ., norm = FALSE))# do is to carry out an operation on a df
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
  
  df.fittedVals <- df.fittedVals %>% tidyr::gather()
  colnames(df.fittedVals)<-c('sample','fitted_values')
  
  df.fittedValst<-df.fittedVals %>% dplyr::group_by(sample)
  
  ## calculate ratios between the fitted curves and the median values
  df.out <- df.median %>%
    data.frame(dplyr::full_join(df.fittedValst,df.median,'sample')) %>%
    dplyr::mutate(correction = df.fittedValst$fitted_values / value) %>% dplyr::select('sample','temperature','value','fitted_values','correction')
  
  ## apply normalization factor to data
  df$temperature<-as.numeric(levels(df$temperature))[df$temperature]
  df <- df %>% dplyr::left_join(df.out %>% dplyr::select(sample,temperature,correction), by = c('sample', 'temperature')) %>%
    dplyr::mutate(norm_value = value * correction) %>% dplyr::ungroup(.)
  
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
#'
#' @importFrom tibble rowid_to_column
#'
#' @export



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
    dplyr::select(Accession, sample, params) %>%
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
                      algorithm = "brute-force",#note: check other ones
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
        dplyr::summarise(av = mean(y, na.rm = T), sd = ifelse(dplyr::n() == 1, 0, sd(y, na.rm = T)))
      
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
# d <- readRDS("tppData_Cliff_preprocessed.Rds")
# #standardize names
# d<-dplyr::rename(d,C=temperature,I=relAbundance,CC=compoundConcentration)
# 
# d$LineRegion<-NA
# 


# #prepare a list of datasets
# datalist<-unique(d$dataset)
# #prepare a list of concentration values
# conclist<-unique(d$CC)
# #prepare a list of proteins
setwd("~/CS7290-/internal data")
# # df<- read_excel("eFT_30K_all5samples_PROTEINS.xlsx")
df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4','F5'), sample_name = c('MEK_1','MEK_2', 'MEK_3','DMSO_1','DMSO_3'), stringsAsFactors = FALSE)

f<-"~/CS7290-/internal data/eFT_30K_all5samples_PROTEINS.xlsx"
df_raw <- read_cetsa(f)
df_clean <- clean_cetsa(df_raw, temperatures = df.temps, samples = df.samples)#assigns temperature and replicate values
df_norm <- normalize_cetsa(df_clean , df.temps$temperature) %>% unique(.)#normalizes according to Franken et. al. without R-squared filter
df_norm <- df_norm %>% dplyr::select(-value,-correction)
#rename columns
df_norm <- df_norm %>% dplyr::rename(dataset = "sample")
colnames(df_norm)<-c("uniqueID","dataset","C","I")
<<<<<<< HEAD
d<-df_norm %>% unique(.) %>% dplyr::ungroup(.)#get unique data
d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)#concentration values are defined in uM
d$dataset<-ifelse(d$dataset=="F4" | d$dataset=="F5","vehicle","treated")#dataset is defined by PD sample definition
DF<-d %>% unique(.) %>% dplyr::group_split(uniqueID) #split null dataset only by protein ID
d_<-d %>% dplyr::filter(CC == 0) %>%  unique(.) %>% dplyr::group_split(uniqueID,dataset) #split vehicle dataset
d_1<-d %>% dplyr::filter(CC > 0) %>% unique(.) %>% dplyr::group_split(uniqueID,dataset) #split treated dataset
=======
d<-df_norm %>% unique(.) %>% dplyr::ungroup(.)#save unique data 
d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)#concentration
d$dataset<-ifelse(d$dataset=="F4" | d$dataset=="F5","vehicle","treated")#dataset
DF<-d %>% unique(.) %>% dplyr::group_split(uniqueID) 
d_<-d %>% dplyr::filter(CC == 0) %>%  unique(.) %>% dplyr::group_split(uniqueID,dataset) 
d_1<-d %>% dplyr::filter(CC > 0) %>% unique(.) %>%   dplyr::group_split(uniqueID,dataset) 
>>>>>>> 48b0e0212b914fc1957c76becaf225697947f6bc

#keep data with less than 2 missing values
DF<-DF %>% purrr::keep(function(x) length(unique(x$C))>=9)#keep values with at least 9 temperature channels (TMT10)
d_<-d_ %>% purrr::keep(function(x) length(unique(x$C))>=9)
d_1<-d_1 %>% purrr::keep(function(x) length(unique(x$C))>=9)
#convert to data frame for uniqueID presence
DF<-data.table::rbindlist(DF)
d_<-data.table::rbindlist(d_)
d_1<-data.table::rbindlist(d_1)
#make sure uniqueIDs are present for treated and vehicle
CID<-intersect(DF$uniqueID,d_$uniqueID)
CID<-intersect(CID,d_1$uniqueID)
DF<-DF %>% subset(uniqueID %in% CID)
DF$LineRegion<-1
d_<-d_%>% subset(uniqueID %in% CID)
d_$LineRegion<-1
d_1<-d_1%>% subset(uniqueID %in% CID)
d_1$LineRegion<-1
#split dataset into equal-sized lists
DF<-DF %>%  dplyr::group_split(uniqueID) 
d_<-d_ %>% dplyr::group_split(uniqueID,dataset) 
d_1<-d_1 %>% dplyr::group_split(uniqueID,dataset) 



DLR<-function(d){
  #preallocate final result as a list
  df_n<-vector(mode = "list", length(d))
  df_n[[1]]<-data.frame()
  df1<-df_n
  df2<-df1
  df3<-df1
  df_1<-df_n
  df0<-df_n
  df_0<-df_n
 
  df_1<-lapply(d, function(x) {x %>% dplyr::group_by(C) %>% dplyr::mutate(I=mean(I,na.rm=TRUE))
    
  })
  df_1 <-lapply(df_1,function(x){x %>% dplyr::mutate(C = as.factor(C),I=as.factor(I),CC=as.factor(CC))
  }) 
  LR<-lapply(df_1, function(x) {dplyr::ntile(dplyr::desc(x$I),3)%>%
      as.data.frame(.) %>% dplyr::rename("LineRegion"=".")})
  df_1<-purrr::map(df_1,function(x){x %>% dplyr::select(-LineRegion)})
  #Add ranks to the list
  df_1<-purrr::map2(df_1,LR, function(x,y) {c(x,y) %>% as.data.frame(.)})
  df_1 <-lapply(df_1,function(x){x %>% dplyr::mutate(C = as.factor(C),I=as.data.frame.numeric(I),CC=as.factor(CC),
                                                     LineRegion=as.data.frame.numeric(LineRegion))})
  
  d<-lapply(d,function(x) x %>% dplyr::arrange(C)) 
  df1<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==1) %>% as.data.frame()})
  df2<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==2) %>% as.data.frame()})
  df3<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==3) %>% as.data.frame()})
  
  change1<-lapply(df1,function(x){x %>% tail(x$C,1)})
  change2<-lapply(df2,function(x){x %>% tail(x$C,1)})
  for(i in 1:length(d)){                                                                                                                                                                                                                                                        
    
    #Inital guess:define the changepoint as the last points of regions 1 and 2
    change1<-tail(df1[[i]]$C,1)
    # 
    change2<-tail(df2[[i]]$C,1)# 
    # 
    
    #Set mean and variance of the top plateau
    Eq1<-lapply(df1[[i]]$I, mean, na.rm = TRUE)
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

  for(i in seq_along(df_0)){
    df_0[[i]]<-df_0[[i]] %>%dplyr::arrange(C)
    ctest<-NA
    ctest<-df_0[[i]] %>%dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup()
    #get the data points that belong to two regions
    Split<-ctest %>%subset(n<(max(n))) %>% data.frame()
    Split<-Split %>% subset(Split$C[duplicated(Split$C)] %in% Split$C)
    Split$LineRegion<-as.numeric(Split$LineRegion)
    
    dap<-list()
    
    dap[[i]]<-df_0[[i]] %>% subset(df_0[[i]]$C %in% Split$C) 
    dap[[i]]$LineRegion<-as.numeric(dap[[i]]$LineRegion)
    
    
    Dap<-dap[[i]] %>%dplyr::group_split(C)
    dap<-lapply(Dap,function(x) x %>% dplyr::mutate(LineRegion=max(LineRegion)))
    dap<-rbindlist(dap) %>% as.data.frame(.)
    df_0[[i]]$LineRegion[df_0[[i]]$C %in% dap$C]<-dap$LineRegion
    df_0[[i]]$LineRegion<-as.factor(df_0[[i]]$LineRegion)
  }
  return(df_0)
}

#preallocate list
results<-vector(mode = "list", length(d_))
results_t<-vector(mode = "list",length(d_1))
results_n<-vector(mode = "list",length(DF))
gdc <- memoise::cache_filesystem("~/Google Drive/.rcache")
results<-suppressWarnings(DLR(d_))#First guess at line regions
results_t<-suppressWarnings(DLR(d_1))
results_n<-suppressWarnings(DLR(DF))
#keep only IDs found in both datasets

results<-results %>% lapply(function(x) na.omit(x))
results_t<-results_t %>% lapply(function(x) na.omit(x))
results_n<-results_n %>% lapply(function(x) na.omit(x))

 #reassign changepoints shared between line regions
df_<-suppressWarnings(CP(results))
df_1<-suppressWarnings(CP(results_t))# 
DFN<-suppressWarnings(CP(results_n))# 
#remove results to save space 
rm(results,results_t,results_n,d_,d_1,DF)

#boostrap function
#returns bootstrapped intensities
#n is the sampling window
#N is the number of iterations
BStrap<-function(Data0,n,N){
  n<-n
  N<-N
  BS<-NA
  BA<-NA
  BSvar<-list(NA)
  BSvar[[1]]<-data.frame()
  Data1<-lapply(Data0,function(x)x %>%dplyr::group_by(C) %>% 
                  dplyr::summarise(I=mean(I),LineRegion=LineRegion)%>% unique(.))

  L<-purrr::map2(Data0,Data1,function(x,y)
  
  BS<-lapply(Data0, function(x){x
    boot::boot.ci(x$I,conf=0.95,type="norm",R=n)})
  #Bootstrap"sample with replacement":
  BS<-lapply(Data0, function(x)x %>% dplyr::select(uniqueID,C,I,LineRegion) %>%  dplyr::sample_n(n,replace=TRUE))
  #generate mean intensities per $C value
  BS<-lapply(BS,function(x) x[order(x$C),])
  
  return(BS)
}


BSvarN<- suppressWarnings(BStrap(DFN,1000,1000))
BSvar<-suppressWarnings(BStrap(df_,1000,1000))
BSvar1<-suppressWarnings(BStrap(df_1,1000,1000))

BSvar<-lapply(BSvar,function(x)x[order(x$C),])
BSvar1<-lapply(BSvar1,function(x)x[order(x$C),])
BSvarN<-lapply(BSvarN,function(x)x[order(x$C),])
#append names to list 

df_<-lapply(df_,function(x)x[order(x$C),])
df_1<-lapply(df_1,function(x)x[order(x$C),])
DFN<-lapply(DFN,function(x)x[order(x$C),])



tlstat<-function(DF,df,df1,PI=FALSE){
  #preallocate list results
  #control/vehicle
  i<-1
  DF<-DF
  df<-df
  df1<-df1
  
  mean1<-list()
  mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="vehicle",uniqueID=df[[i]]$uniqueID[1],Tm=rep(0,1))
  if(!isTRUE(PI)){
    for(i in 1:length(df)){
      df[[i]]<-unique(df[[i]])
      mean1[[i]]<-df[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>%
        dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                      Tm=with(df[[i]], stats::approx(df[[i]]$I,df[[i]]$C, xout=max(df[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=purrr::map(M1,function(x){deviance(x)}),
                      Rsq=purrr::map(M1,function(x){summary(x)$r.squared}), 
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="vehicle",
                      uniqueID=df[[i]]$uniqueID[1])
      
      
    }
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    for(i in 1:length(df1)){
      df1[[i]]<-unique(df1[[i]])
      mean1_1[[i]]<-df1[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>% 
        dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                      Tm=with(df1[[i]], stats::approx(df1[[i]]$I,df1[[i]]$C, xout=max(df1[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=map(M1,function(x){deviance(x)}),
                      Rsq=map(M1,function(x){summary(x)$r.squared}), 
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="treated",
                      uniqueID=df1[[i]]$uniqueID[1])
      
      
      
    }
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    DF<-lapply(DF,function(x) x %>% na.omit())
    for(i in 1:length(DF)){
      DF[[i]]<-unique(DF[[i]])
      mean3[[i]]<-DF[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>% 
        dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                      Tm=with(DF[[i]], stats::approx( DF[[i]]$I,DF[[i]]$C, xout=max(DF[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=map(M1,function(x){deviance(x)}),
                      Rsq=map(M1,function(x){summary(x)$r.squared}),
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="null",
                      uniqueID=DF[[i]]$uniqueID[1])
      
      
      
      
    }
    mean1<-lapply(mean1,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    mean1_1<-lapply(mean1_1,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    mean3<-lapply(mean3,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean3<-mean3 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    
    #split by uniqueIDs
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
    results<-rlist::list.rbind(c(mean1,mean1_1,mean3)) %>% dplyr::group_split(uniqueID)
  } else if (isTRUE(PI)){
    for(i in 1:length(df)){
      df[[i]]<-unique(df[[i]])
      mean1[[i]]<-df[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>%
        dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                      Tm=with(df[[i]], stats::approx(df[[i]]$I,df[[i]]$C, xout=max(df[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=map(M1,function(x){deviance(x)}),
                      Rsq=map(M1,function(x){summary(x)$r.squared}), 
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="vehicle",
                      uniqueID=df[[i]]$uniqueID[1])
      
      
    }
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    for(i in 1:length(df1)){
      df1[[i]]<-unique(df1[[i]])
      mean1_1[[i]]<-df1[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>% 
        dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                      Tm=with(df1[[i]], stats::approx(df1[[i]]$I,df1[[i]]$C, xout=max(df1[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=map(M1,function(x){deviance(x)}),
                      Rsq=map(M1,function(x){summary(x)$r.squared}), 
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="treated",
                      uniqueID=df1[[i]]$uniqueID[1])
      
      
      
    }
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1), dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    
    for(i in 1:length(DF)){
      DF[[i]]<-unique(DF[[i]])
      mean3[[i]]<-DF[[i]] %>% data.frame(.) %>% 
        dplyr::group_nest(LineRegion,uniqueID) %>% 
        dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                      CI=map(M1,function(x){predict(x,interval="confidence")[,1]} %>% data.frame(.)),
                      Tm=with(DF[[i]], approx( DF[[i]]$I,DF[[i]]$C, xout=max(DF[[i]]$I, na.rm=TRUE)-0.5))$y,
                      slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                      intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                      rss=map(M1,function(x){deviance(x)}),
                      Rsq=map(M1,function(x){summary(x)$r.squared}), 
                      AUC = purrr::map2(data,CI,function(x,y) pracma::trapz(x$C,y[,'fit'])),
                      dataset="treated",
                      uniqueID=DF[[i]]$uniqueID[1])
      
      
      
    }
    mean1<-lapply(mean1,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest_legacy(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    mean1_1<-lapply(mean1_1,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest_legacy(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    mean3<-lapply(mean3,function(x) x %>% dplyr::select(-data,-CI) %>% tidyr::unnest_legacy(Tm=Tm,slope=slope,intercept=intercept,rss=rss,Rsq=Rsq))
    
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean3<-mean3 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    
    #split by uniqueIDs
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
    results<-rlist::list.rbind(c(mean1,mean1_1,mean3)) %>% dplyr::group_split(uniqueID)
  }
  results
}


spstat<-function(DF,df,df1,PI=FALSE){
  #plot spline results
  if(!isTRUE(PI)){
    i<-1
    mean1<-list()
    mean1[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC = rep(0,1),dataset="vehicle",uniqueID=df[[i]]$uniqueID[1])
    df<-df %>% purrr::keep(function(x)length(unique(x$C))>9)
    mean1<- lapply(df,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                         spar=as.numeric(M1[[1]]$spar),
                                                         Tm=with(df[[i]], stats::approx(df[[i]]$I,df[[i]]$C, xout=max(df[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                         lambda=as.numeric(M1[[1]]$lambda),
                                                         df =round(as.numeric(M1[[1]]$df),3),
                                                         rss = as.numeric(M1[[1]]$pen.crit),
                                                         knots = M1[[1]]$fit$nk,
                                                         AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                         dataset="vehicle",
                                                         uniqueID=x$uniqueID[1],
                                                         n = nrow(x)))
    
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1])
    
    df1<-df1 %>% purrr::keep(function(x)length(unique(x$C))>9)
    mean1_1<- lapply(df1,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                            spar=as.numeric(M1[[1]]$spar),
                                                            Tm=with(df1[[i]], stats::approx(df1[[i]]$I,df1[[i]]$C, xout=max(df1[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                            lambda=as.numeric(M1[[1]]$lambda),
                                                            df =round(as.numeric(M1[[1]]$df),3),
                                                            rss = as.numeric(M1[[1]]$pen.crit),
                                                            knots = M1[[1]]$fit$nk,
                                                            AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                            
                                                            dataset="treated",
                                                            uniqueID=x$uniqueID[1],
                                                            n = nrow(x)))
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1])
    
    DF<-DF %>% purrr::keep(function(x)length(unique(x$C))>9)
    mean3<- lapply(DF,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                         spar=as.numeric(M1[[1]]$spar),
                                                         Tm=with(DF[[i]], stats::approx(DF[[i]]$I,DF[[i]]$C, xout=max(DF[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                         lambda=as.numeric(M1[[1]]$lambda),
                                                         df =round(as.numeric(M1[[1]]$df),3),
                                                         rss = as.numeric(M1[[1]]$pen.crit),
                                                         knots = M1[[1]]$fit$nk,
                                                         AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                         dataset="null",
                                                         uniqueID=x$uniqueID[1],
                                                         n = nrow(x)))
    
    
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean3<-mean3 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    
    #split by uniqueIDs
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
    results<-rlist::list.rbind(c(mean1,mean1_1,mean3)) %>% dplyr::group_split(uniqueID)
  } else if (isTRUE(PI)){
    i<-1
    mean1<-list()
    mean1[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC = rep(0,1),dataset="vehicle",uniqueID=df[[i]]$uniqueID[1])
    
    mean1<- lapply(df,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                         spar=as.numeric(M1[[1]]$spar),
                                                         Tm=with(df[[i]], stats::approx(df[[i]]$I,df[[i]]$C, xout=max(df[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                         lambda=as.numeric(M1[[1]]$lambda),
                                                         df =round(as.numeric(M1[[1]]$df),3),
                                                         rss = as.numeric(M1[[1]]$pen.crit),
                                                         knots = M1[[1]]$fit$nk,
                                                         AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                         dataset="vehicle",
                                                         uniqueID=x$uniqueID[1],
                                                         n = nrow(x)))
    
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC= rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1])
    
    
    mean1_1<- lapply(df1,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                            spar=as.numeric(M1[[1]]$spar),
                                                            Tm=with(df1[[i]], stats::approx(df1[[i]]$I,df1[[i]]$C, xout=max(df1[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                            lambda=as.numeric(M1[[1]]$lambda),
                                                            df =round(as.numeric(M1[[1]]$df),3),
                                                            rss = as.numeric(M1[[1]]$pen.crit),
                                                            knots = M1[[1]]$fit$nk,
                                                            AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                            dataset="treated",
                                                            uniqueID=x$uniqueID[1],
                                                            n = nrow(x)))
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(spar=rep(0,1),Tm=rep(0,1),lambda=rep(0,1),df=rep(0,1),rss=rep(0,1),knots=rep(0,1),AUC=rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1])
    
    
    mean3<- lapply(DF,function(x) x %>% dplyr::summarise(M1 = list(stats::smooth.spline(x=x$C,y=x$I)),
                                                         spar=as.numeric(M1[[1]]$spar),
                                                         Tm=with(DF[[i]], stats::approx(DF[[i]]$I,DF[[i]]$C,ties = "ordered", xout=max(DF[[i]]$I, na.rm=TRUE)-0.5))$y,
                                                         lambda=as.numeric(M1[[1]]$lambda),
                                                         df =round(as.numeric(M1[[1]]$df),3),
                                                         rss = as.numeric(M1[[1]]$pen.crit),
                                                         knots = M1[[1]]$fit$nk,
                                                         AUC = pracma::trapz(predict(M1[[1]]$fit)$x,predict(M1[[1]]$fit)$y),
                                                         dataset="null",
                                                         uniqueID=x$uniqueID[1],
                                                         n = nrow(x)))
    
   
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    mean3<-mean3 %>% as.data.frame(.) %>% subset(uniqueID %in% CID)
    
    #split by uniqueIDs
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
   
  }
  #mutate(d1 = paramsAlternative - paramsNull,
  #d2 = n1 - paramsAlternative
  
  ma<-purrr::map2(mean1,mean1_1,function(x,y) x %>% dplyr::mutate(n1=sum(x$n,y$n),DTm=y$Tm-x$Tm,RSS=y$rss+x$rss,DF=x$df+y$df,DAUC=y$AUC-x$AUC))
  testResults <-purrr::map2(ma,mean3,function(x,y) x %>%
                              #(rssDiff/d1)/(rss1/d2)
                              dplyr::mutate(ftest = ((x$RSS-y$rss)/(x$DF-y$df))/(y$rss/(x$n1-x$DF)),
                                            d1 = (x$DF-y$df),
                                            d2 = (x$n1-x$DF),
                                            pV = ifelse(is.na(stats::pf(ftest, df1 = d1,
                                                                        df2 = d2)),NA,1 - stats::pf(ftest, df1 = (x$DF-y$df), df2 = (x$n1-x$DF))),
                                            pA = stats::p.adjust(pV,"BH"))
                            )
  
  check<-dplyr::bind_rows(testResults) %>% data.frame(.)
  ggplot(check)+
    geom_density(aes(x=ftest),fill = "steelblue",alpha = 0.5) + 
    geom_line(aes(x=ftest,y= df(ftest,df1=d1,df2=d2)),color="darkred",size = 1.5) +
    theme_bw() +
    coord_cartesian(xlim=c(0,10))
  
 

  #t1<-check[t$.>1]#remove negative or F test values <1
  
  
  
  mean1<-mean1 %>% purrr::keep(function(x) x$uniqueID %in% check)
  mean1_1<-mean1_1 %>% purrr::keep(function(x) x$uniqueID %in% check)
  mean3<-mean3 %>% purrr::keep(function(x) x$uniqueID %in% check)
  results<-rlist::list.rbind(c(mean1,mean1_1,mean3)) %>% dplyr::group_split(uniqueID)
  results
}
#gettrilinear results
tlresults<-list()
tlresults_PI<-list()
#confidence intervals
tlresults<-tlstat(DFN,df_,df_1,PI=FALSE)#place null, vehicle and treated lists with no prediction intervals

tlf<-function(tlresults){
##Apply Filters
#####################

tlresults1<-tlresults#save unfiltered data
#apply filters prior to hypothesis testing
tlresults<-tlresults %>% keep(function(x) min(as.numeric(x$Rsq),na.rm=TRUE) >= 0.40)#the linear region have the largest slope < 0.03
tlresults<-tlresults %>% keep(function(x) mean(as.numeric(x$slope),na.rm=TRUE) <= -0.02)
#tlresults<-tlresults %>% keep(function(x)  sum(data.frame(x)[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <10)#move data with extremely large RSS values 
# tlresults<-tlresults %>% keep(function(x) sum(data.frame(x)[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <1.3)
tlresults<-tlresults %>% keep(function(x) sum(data.frame(x)[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) > sum(data.frame(x)[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE))#remove data with extremely large RSS values 
tlresults<-tlresults %>% keep(function(x) mean(data.frame(x)$Tm[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle")],na.rm=TRUE) < mean(data.frame(x)$Tm[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated")],na.rm=TRUE))
#tlresults<-tlresults %>% keep(function(x) max(data.frame(x)$slope[x$LineRegion==2],na.rm=TRUE) < -0.03)#the linear region have the largest slope < 0.03
#tlresults<-tlresults %>% keep(function(x) length(x$slope)>8)#remove list values with less than 5 rows
#tlresults<-tlresults %>% keep(function(x) abs(max(x$slope[!x$LineRegion==2] ,na.rm=TRUE)) < 0.1)#eeps plateau values where the min abs(slope) < 0.06
#steepest slope in vehicle and treatment has to be less than 0.06C

Nsum<-list()
Nsum[[1]]<-data.frame(RSS=0,Tm=0)

Nsum<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "null")) %>% 
               dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))

#get the summed rss values for vehicle
Rssv<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% 
               dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>%head(.,1))
#get the summed rss values for treated
Rsst<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% 
               dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(.$rss))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
#find the rss difference between treated and vehicle 

Rssv<-lapply(Rssv,function(x)na.omit(x))
Rsst<-lapply(Rsst,function(x)na.omit(x))
#find common IDs
CID<-intersect(rbindlist(Rsst)$uniqueID,rbindlist(Rssv)$uniqueID)
#keep common IDs
Rssv<-Rssv %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID)) 
Rsst<-Rsst %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))
Nsum<-Nsum %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))                           
K1<-data.frame(rbindlist(purrr::map2(Rsst,Rssv,function(x,y) data.frame(RSSd = x$RSS-y$RSS, Tma = x$Tm - y$Tm)))) 
K2<-rbindlist(Rssv)$uniqueID %>% as.data.frame(.)
Dsum<-data.frame(K1,K2)
names(Dsum)<-c("RSSd","Tma","uniqueID")
Dsum<-Dsum %>% dplyr::mutate(rank = dplyr::ntile(Dsum$Tma,7))
#keep data where the difference in RSS is less than the null
#nsum converted to data frame
Nsum<-data.frame(RSSn=data.table::rbindlist(Nsum))
names(Nsum)<-c("RSSn","Tmn","dataset","uniqueID")
Nsum<-Nsum %>% dplyr::filter(uniqueID %in% CID)
Nsum<-Nsum %>% dplyr::mutate(id=rownames(Nsum))

Nsum$dataset<-as.factor(Nsum$dataset)
#mutate data frame
#Cliff 1D
Dsum1<-Dsum %>% dplyr::left_join(Nsum,by = c("uniqueID"="uniqueID"))
#Childs
Dsum2<-Dsum %>% dplyr::right_join(Nsum,by = c("uniqueID"="uniqueID"))

Dsum<-Dsum2
Dsum$RSSd<-Dsum1$RSSd
Dsum$Tma<-Dsum1$Tma
<<<<<<< HEAD
Dsum<-Dsum %>% dplyr::mutate(rank = dplyr::ntile(Dsum$Tma,7))
=======
Dsum<-Dsum %>% dplyr::mutate(rank = ntile(Dsum$Tma,7))
<<<<<<< HEAD
>>>>>>> 48b0e0212b914fc1957c76becaf225697947f6bc
=======
>>>>>>> 48b0e0212b914fc1957c76becaf225697947f6bc
Dsum<-dplyr::arrange(Dsum, dplyr::desc(Tma), dplyr::desc(RSSd))  %>% dplyr::filter(RSSd>0) 

test<-data.frame()
test<-Dsum[which(Dsum$RSSn>Dsum$RSSd),] %>% data.frame()#get the stable proteins (+ = Rsstreated-Rssvehicle)

# rssdec<-data.frame()
# rssdec<-data.frame(data.table::fsort(test$RSSd,decreasing=TRUE))#decreasing Rss differences
# names(rssdec)<-"Rssd"
# 
tmdec<-data.frame()
tmdec<-data.table::fsort(test$Tma,decreasing=TRUE) %>% data.frame()
names(tmdec)<-"Tm"

test<-tmdec %>% inner_join(test,by=c("Tm"="Tma"))#orders data by decreasing Tm


orows<-data.frame()
orows <- test

orows$id<-sapply(orows$id, function(x) as.numeric(as.character(x)))
#order by RSS differences while keeping original rownames for index
#create an external data frame for stabilized proteins

Df1<-tlresults[orows$id] #divide 1=highly destabilized,4=noeffect,7=highly stabilized
df1<-list()
#get uniqueID and dataset for stable proteins with decreasing RSS differences
df1<-lapply(Df1,function(x) x %>% dplyr::select(uniqueID,dataset) %>% head(.,1))
df1<-data.frame(rbindlist(df1))

#unlist to data.frame
#order the original data by RSS differences
#
df2<-rbindlist(DFN) %>% as.data.frame(.) %>% right_join(df1, by = c('uniqueID' = 'uniqueID') ) %>% dplyr::select(-dataset.y)
names(df2)<-c("uniqueID","Dataset","C","I","CC","LineRegion")
return(list(df1,df2,Df1))
}
#return filtered lists
res<-tlf(tlresults)
tlCI<-function(i,df1,df2,Df1,overlay=TRUE){
  null<-data.frame()
  i<-i
  df1<-df1
  Df1<-Df1[[i]]
  
  
  DF1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated'))
  
  null<-Df1 %>% subset(dataset == "null")
  
  pred1<-predict(null$M1[[1]], interval="confidence")
  pred2<-predict(null$M1[[2]], interval="confidence")
  pred3<-predict(null$M1[[3]], interval="confidence")
  Pred1<-NA
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  
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
  
  PLN<-ggplot2::ggplot(Pred, ggplot2::aes(x = C,y = fit,color=Treatment)) +
    ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1$uniqueID[1],"null"))+
    ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ 
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
    ggplot2::annotate("text", x=62, y=1, label= paste("RSS= ",round(sum(null$rss),3)))
  

  
  DF_f<-df2 %>%subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated')) %>% subset(dataset=="vehicle")
  
  vehicle<-Df1 %>% subset(dataset == "vehicle")
  
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
  DF_f1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated')) %>% subset(dataset=="treated")
  
  
  treated<-data.frame()
  treated<-Df1 %>% subset(dataset == "treated")
  
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
  #Residuals
  rn<-data.frame(c(residuals(null$M1[[1]]),residuals(null$M1[[2]]),residuals(null$M1[[3]])))
  Pred<-cbind(Pred,rn)
  names(Pred)<- names(Pred)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]]),residuals(vehicle$M1[[3]])))
  Pred1<-cbind(Pred1,rn)
  names(Pred1)<- names(Pred1)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]]),residuals(treated$M1[[3]])))
  Pred2<-cbind(Pred2,rn)
  names(Pred2)<- names(Pred2)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  Preds<-rbind(Pred1,Pred2)
  PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = residuals,color=Treatment)) +
    ggplot2::geom_point()+ ggplot2::ggtitle(paste(Df1$uniqueID[1]))+
    ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
  print(PLrs)
  PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
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
    PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
      ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Df1$uniqueID[1],"alternative"))+
      ggplot2::annotate("text", x=62, y=1, label= paste("\u03A3","RSS= ",round(sum(Df1 %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% dplyr::select(rss) %>% sum(.),Df1 %>% dplyr::filter(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% dplyr::select(rss) %>% sum(.)),3)))+
      ggplot2::annotate("text", x=62, y=0.9, label=  paste("\u0394", "AUC = ",AUCd))+ ggplot2::annotate("text", x=62, y=0.8, label= paste("\u0394","Tm = ",Tm_d,"\u00B0C"))
    #bquote(Value~is~sigma~R^{2}==.(r2.value)))
    par(mfrow=c(2,2))
    print(PLR_P2)
  }else if(overlay=="FALSE"){
    PLR<-PLR_P2+ggplot2::geom_point(data=Pred,mapping=ggplot2::aes(x=C,y=I))+ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ggplot2::ggtitle(paste(Df1$uniqueID[1],"alternative"))+facet_wrap("Treatment") 
    print(PLR)
  }
}
i=2

plotTL<-tlCI(i,res[[1]],res[[2]],res[[3]],overlay=TRUE)
#df1 <- only IDs in order desc(stability)
#df2<-original data in order  
#Df1 <- ordered spline results 
###############################
#get spline results
spresults<-list()
spresults_PI<-list()

spresults<-spstat(DFN,df_,df_1,PI=FALSE)

spf<-function(spresults){
spresults1<-spresults
##############################
#Apply filters 
#RSS null > RSS treated
sp<-spresults %>% keep(function(x) x$rss[x$dataset=="null"]>sum(x$rss[!x$dataset=="null"],na.rm=TRUE))
#keep the stabilized shifts
spresults<-sp %>% keep(function(x) x$Tm[x$dataset=="treated"]<x$Tm[x$dataset=="vehicle"])
#keep the positive AUC differences
spresults<-spresults %>% keep(function(x) mean(x$AUC[x$dataset=="treated"],na.rm=TRUE)>mean(x$AUC[!x$dataset=="vehicle"],na.rm=TRUE))
spresults<-spresults %>% keep(function(x) max(x$lambda)<1)

#get Tm and RSS differences
sp<-lapply(spresults, function(x) x %>% dplyr::mutate(Tmd= x$Tm[x$dataset == "vehicle"] - x$Tm[x$dataset=="treated"],
                                                      RSSd = x$rss[x$dataset == "null"] - sum(x$rss[!x$dataset=="null"]),
                                                      AUCd = x$AUC[x$dataset == "treated"]- x$AUC[x$dataset == "vehicle"]))
#conserve list indexes
sl<-lapply(seq_along(sp),function(x) as.numeric({paste(x)}))
#insert list index column
sp<-map2(sp,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))
sp<-data.table::rbindlist(sp) %>% as.data.frame(.)
sp<-dplyr::dplyr::arrange(sp,dplyr::desc(AUCd),dplyr::desc(RSSd),dplyr::desc(Tmd)) %>% dplyr::select(uniqueID,id) %>% unique(.) 
#arrange results by decreasing AUCd, RSSd and Tmd and standardize the order in spresults
#Df1 holds the model results and stats for splines 

df1<-sp$uniqueID
df2<-rbindlist(DFN) %>% data.frame(.) 
df2<-sp%>% inner_join(df2, by = c("uniqueID"="uniqueID"))
Df1<-spresults[sp$id]
return(list(df1,df2,Df1))
}
res_sp<-spf(spresults)
# #filter bootstrapped data
# BSvar <- BSvar %>% keep(function(x) x$uniqueID[1] %in% df1)
# BSvar1<- BSvar1%>% keep(function(x) x$uniqueID[1] %in% df1)
# BSvarN<- BSvarN%>% keep(function(x) x$uniqueID[1] %in% df1)

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}
###############################
spCI<-function(i,df1,df2,Df1,overlay=TRUE,alpha){
  null<-data.frame()
  i<-i
  df2$uniqueID<-as.character(df2$uniqueID)
  #get original data
  ###########################################
  df1<-df1[i]
  DF1<-df2%>% subset(uniqueID %in% df1)
  null<-Df1[[i]] %>% subset(uniqueID %in% df1& dataset %in% "null")
  ###########################################
  DF_f<-df2%>% subset(uniqueID %in% df1 & dataset %in% "vehicle")
  vehicle<-Df1[[i]] %>% subset(uniqueID == df1 & dataset == "vehicle")
  ###########################################
  DF_f1<-df2%>% subset(uniqueID %in% df1 & dataset %in% "treated")
  treated<-Df1[[i]] %>% subset(uniqueID == df1 & dataset == "treated")
  
  ###########################################
  #get confidence intervals for all conditions
  ###########################################
  
  #return fit and confidence intervals
    
    BSVarN<-NA
    BSVar<-NA
    BSVar1<-NA
    BSvarN<-NA
    BSvar<-NA
    BSvar1<-NA
    BSvarN<-df2 %>% subset(uniqueID == df1 ) 
    BSvar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")
    BSvar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")
    BSVarN<-df2 %>% subset(uniqueID == df1 ) %>% dplyr::group_by(C) %>% dplyr::summarise(I=mean(I))
    BSVar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")%>% dplyr::group_by(C) %>% dplyr::summarise(I=mean(I))
    BSVar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")%>% dplyr::group_by(C) %>% dplyr::summarise(I=mean(I))
    BSVarN<-BSVarN %>% dplyr::mutate(dataset="null") 
    BSVar<-BSVar %>% dplyr::mutate(dataset="vehicle")
    BSVar1<-BSVar1 %>% dplyr::mutate(dataset="treated")
    BSVarN$dataset<-as.factor(BSVarN$dataset)
    BSVar$dataset<-as.factor(BSVar$dataset)
    BSVar1$dataset<-as.factor(BSVar1$dataset)
    fit <-  stats::smooth.spline(x = BSVar$C, y=BSVar$I,cv=TRUE)
    fit1<-  stats::smooth.spline(x = BSVar1$C, y=BSVar1$I,cv=TRUE)
    fitN<-  stats::smooth.spline(x = BSVarN$C, y=BSVarN$I,cv=TRUE)
    
    #####try GAM
    #fit penalized splines
    m <- mgcv::gam(I ~ s(C), data = BSVar , method = "ML")
    m1<-  mgcv::gam(I ~ s(C), data =BSVar1, method = "ML")
    mn<-  mgcv::gam(I ~ s(C), data = BSVarN, method = "ML")
    
    #####try GAM
    
    #Plot boostrapped  residuals with 95%CI
    #PLP<-plot(m, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
    #generate random values from a multivariate normal distribution
   
    #get some parmeters
    Vb <- vcov(m)
    newd <- with(BSVar, data.frame(C = seq(min(C), max(C), length = 20)))
    pred <- predict(m, newd, se.fit = TRUE)
    se.fit <- pred$se.fit
    #get some parmeters
    Vb1<- vcov(m1)
    newd1<- with(BSVar1,data.frame(C = seq(min(C), max(C), length = 30)))
    pred1<- predict(m1,newd1,se.fit = TRUE)
    se.fit1<- pred1$se.fit
    #generate std
    set.seed(42)
    N <- 10000
    #sample n from mvn dist
    BUdiff <- mgcv::rmvn(N, mu = rep(0, nrow(Vb)), Vb )
    #sample n from mvn dist
    BUdiff1<-  mgcv::rmvn(N, mu = rep(0, nrow(Vb1)),Vb1)
    #calculate deviation
    Cg <- predict(m, newd, type = "lpmatrix")
    simDev <- Cg %*% t(BUdiff)
    #calculate deviation
    Cg1<- predict(m1,newd1,type = "lpmatrix")
    simDev1<- Cg1%*% t(BUdiff1)
    #calculate abs deviation
    absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
    #calculate abs deviation
    absDev1<- abs(sweep(simDev1,1, se.fit, FUN = "/"))
    #max abs dev
    masd <- apply(absDev, 2L, max)
    #max abs dev
    masd1<- apply(absDev1,2L, max)
    #95% crit values
    crit <- quantile(masd, prob = alpha/2)
    #95% crit values
    crit1<- quantile(masd1,prob = alpha/2)
    #plot CI
    pred <- transform(cbind(data.frame(pred), newd),
                      uprP = fit + (2 * se.fit),
                      lwrP = fit - (2 * se.fit),
                      uprS = fit + (crit * se.fit),
                      lwrS = fit - (crit * se.fit))
    pred$dataset<-"vehicle"
    pred$dataset<-as.factor(pred$dataset)
    pred$CI<-"vehicle"
    pred$CI<-as.factor(pred$CI)

    plot<-ggplot(pred,mapping= ggplot2::aes(x = C,y=fit ))+
      geom_point(BSVar, mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      geom_ribbon(aes(ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2) +
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle("")
    
    pred1<- transform(cbind(data.frame(pred1),newd1),
                      uprP = fit + (2 * se.fit),
                      lwrP = fit - (2 * se.fit),
                      uprS = fit + (crit * se.fit),
                      lwrS = fit - (crit * se.fit))
    pred1$dataset<-"treated"
    pred1$dataset<-as.factor(pred1$dataset)
    pred1$CI<-"treated"
    pred1$CI<-as.factor(pred1$CI)
    pred1$AUC<-pracma::trapz(pred1$fit)-pracma::trapz(pred1$fit)
    pred1$AUC<-round(pred1$AUC[1],3)
    if( pred1$AUC[1] > 5){
      pred1$AUC<-pracma::trapz(pred1$lwrP)-pracma::trapz(pred$uprP)#AUC diff in stabilized CI
      
      pred1$AUC<-pracma::trapz(pred1$lwrP)-pracma::trapz(pred$uprP)#AUC diff in stabilized CI
      
    }else if ( pred1$AUC[1]< -5){
      pred1$AUC<-pracma::trapz(pred$lwrP)-pracma::trapz(pred1$uprP)#AUC diff in destabilized CI
    }else{
      pred1$AUC<-pracma::trapz(pred1$fit)-pracma::trapz(pred$fit) #AUC diff in fit
    }
    pred1$Tm<-stats::approx(pred1$fit,pred1$C,0.5)$y-stats::approx(pred$fit,pred$C,0.5)$y# Tm difference (+ stabilized)
    pred1$RSS<- deviance(m1)-deviance(m)#RSS diff(+ stabilized)
    pred1$Tm<-round(pred1$Tm[1],1)
    pred1$RSS<- round(pred1$RSS,3)
    #Residuals
    
    Pred<-data.frame(m$fitted.values,m$residuals)
    names(Pred)<- c("fit","rn")
    Pred$dataset<-as.factor("vehicle")
    BSVar$dataset<-as.factor("vehicle")
    Pred1<-data.frame(m1$fitted.values,m1$residuals)
    names(Pred1)<-c("fit","rn")
    Pred1$dataset<-as.factor("treated")
    Preds<-rbind(Pred,Pred1)
    BSVar1$dataset<-as.factor("treated")
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=dataset)) +ggplot2::geom_point()+ ggplot2::ggtitle(paste(Df1$uniqueID[1]))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    print(PLrs)
    plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
      ggplot2::geom_point(BSvar,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(pred,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(c(as.character(df1),"alternative"))+
      ggplot2::annotate("text", x=62, y=1, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])))+
      ggplot2::annotate("text", x=62, y=0.9, label=  paste("\u0394", "AUC = ",round(pred1$AUC[1],3)))+ggplot2::annotate("text", x=62, y=0.8, label= paste("\u0394","Tm = ",round(pred1$Tm[1],3),"\u00B0C"))
    
    plot<-plot1+
      ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
      ggplot2::labs(y = "Relative Solubility",
                    x = "Temperature")
    print(plot)  
  }

#generate 95%CI for splines
Pred<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],overlay=TRUE,0.05)
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
                       seed = NULL, alwaysPermute = FALSE, maxAttempts = 50){
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
#GoF by RSS comparison for each protein
#RSS difference

computeRSSdiff <- function(x,y,treatment,maxAttempts = 100, repeatsIfNeg = 50){
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
      do({
        fit = computeRSS(x=.$x, y = .$y,start = start1, seed=repeats,
                         maxAttempts = maxAttempts,
                         alwaysPermute = alwaysPermute)
        
      }) %>% ungroup
    rss0 <- nullResults$rss
    rss1 <-sum(altResults$rss)
    rssDiff <- rss0-rss1
    
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

  sigfit<-function(df_,df_1,i){
    nlm1<-list(rep(NA,length(df_)))
    nlm2<-nlm1
    dfc<-nlm1
    dfy<-nlm1
    result<-nlm1
    result1<-nlm1
    Pred<-nlm1
    Pred1<-list(rep(NA,length(df_1)))
    df_<-lapply(df_,function(x)x %>% dplyr::mutate(dataset="vehicle"))
    
    df_1<-lapply(df_1,function(x)x %>% dplyr::mutate(dataset="treated"))
    #order by C
    df_<-lapply(df_,function(x)x[order(x$C),])
    df_1<-lapply(df_1,function(x)x[order(x$C),])
    
    df<-lapply(df_,function(x) x %>% dplyr::group_by(C) %>% dplyr::mutate(I = mean(I)))
    df1<-lapply(df_1,function(x) x %>% dplyr::group_by(C) %>% dplyr::mutate(I = mean(I)))
    #sigmoidal fit
    nlm1<-df %>% purrr::keep(function(x)!inherits(fitSingleSigmoid(x$C,x$I),'try-error')) #flag sigmoidal data
    nlm2<-lapply(nlm1,function(x)x %>% group_by(uniqueID) %>% dplyr::mutate(fit = try(fitSingleSigmoid(x$C,x$I)))) #fit sigmoids
    CT<-nlm2 %>% lapply(function(x) !coef(x)['Pl'][[1]]==0)#find values where Pl = 0 
    nlm2<-nlm2[CT==TRUE]#filter nls data with Pl = 0
    nlm1<-nlm1[CT==TRUE]
 
    Ct<-purrr::map(nlm1,function(x) seq(from=min(x$C),to=max(x$C),length.out=length(unique(x$C))))
    Ct<-purrr::map2(nlm1,Ct,function(x,y) rep(y,length(x$C)/length(y))) 
    Ct<-lapply(Ct,function(x)x[order(x)])
    
    dfc <- purrr::map(nlm2,function(x) try(nlstools::confint2(x,level=0.95)))#collect predicted values
    nlm2<-purrr::map2(Ct,dfc,function(x,y)x %>% dplyr::mutate(Pl=y[1],a=y[2],b=y[3],Pl1=y[4],a1=y[5],b1=y[6]))
    result<- purrr::map(nlm2,function(x)x %>%  dplyr::mutate(LOW = (1-as.numeric(x$Pl[1]))/(1+exp((as.numeric(x$b[1])-as.numeric(x$a[1])/x$CP)))+as.numeric(x$Pl[1]), HI = (1-as.numeric(x$Pl1[1]))/(1+exp((as.numeric(x$b1[1])-as.numeric(x$a1[1])/x$CP)))+as.numeric(x$Pl1[1])))#get lower CI
    #result<-purrr::map2(nlm2,result,function(x,y)cbind(x,y))
    
    AUC<-data.table::rbindlist(purrr::map(result,function(x)pracma::trapz(x$C,x$I) %>% as.data.frame(.))) %>% as.data.frame(.)
    Tm<- lapply(result,function(x) try(stats::approx(x$I,x$C,0.5)$y[1]))
    T<- Tm %>% purrr::map(function(x) inherits(x,'try-error'))
    Tm<-Tm[T==FALSE]
    
    result<-purrr::map2(result[1:length(Tm)],Tm,function(x,y)cbind(x,data.frame(y)))
    result<-purrr::map(result[1:length(Tm)],function(x) x %>% dplyr::rename("Tm"="y"))
    result<-purrr::map2(result[1:length(Tm)],AUC[1:length(Tm),],function(x,y)cbind(x,data.frame(y)))
    Pred<-purrr::map(result[1:length(Tm)],function(x) x %>% dplyr::rename("AUC"="y"))

    #sigmoidal fit
    nlm1<-df_1 %>% purrr::keep(function(x)!inherits(fitSingleSigmoid(x$C,x$I),'try-error')) #remove non-sigmoidal data
    nlm2<-lapply(nlm1,function(x)fitSingleSigmoid(x$C,x$I)) #fit sigmoids
    CT<-nlm2 %>% lapply(function(x) !coef(x)['Pl'][[1]]==0)#find values where Pl = 0 
    nlm2<-nlm2[CT==TRUE]#filter nls data with Pl = 0
    nlm1<-nlm1[CT==TRUE]
    
    Ct<-purrr::map(nlm1,function(x) seq(from=min(x$C),to=max(x$C),length.out=length(unique(x$C))))
    Ct<-purrr::map2(nlm1,Ct,function(x,y) rep(y,length(x$C)/length(y))) 
    Ct<-lapply(Ct,function(x)x[order(x)])
    nlm1<-purrr::map2(nlm1,Ct,function(x,y) try(cbind(x,y)))
    
    nlm<-lapply(nlm1,function(x) inherits(x,'try-error'))
    nlm1<-nlm1[nlm==FALSE]
    nlm1<-purrr::map(nlm1,function(x) x %>% dplyr::rename("CP"="y"))
    nlm2<-nlm2[nlm==FALSE]
    dfc <- purrr::map(nlm2,function(x) try(nlstools::confint2(x,level=0.95)))#collect predicted values
    nlm2<-purrr::map2(nlm1,dfc,function(x,y)x %>% dplyr::mutate(Pl=y[1],a=y[2],b=y[3],Pl1=y[4],a1=y[5],b1=y[6]))
    result1<- purrr::map(nlm2,function(x)x %>%  dplyr::mutate(LOW = (1-as.numeric(x$Pl[1]))/(1+exp((as.numeric(x$b[1])-as.numeric(x$a[1])/x$CP)))+as.numeric(x$Pl[1]), HI = (1-as.numeric(x$Pl1[1]))/(1+exp((as.numeric(x$b1[1])-as.numeric(x$a1[1])/x$CP)))+as.numeric(x$Pl1[1])))#get lower CI
    #result<-purrr::map2(nlm2,result,function(x,y)cbind(x,y))
    
    
    AUC<-data.table::rbindlist(purrr::map(result1,function(x)pracma::trapz(x$C,x$I) %>% as.data.frame(.))) %>% as.data.frame(.)
    Tm<- lapply(result1,function(x) try(stats::approx(x$I,x$C,0.5)$y[1]))
    T<- Tm %>% purrr::map(function(x) inherits(x,'try-error'))
    Tm<-Tm[T==FALSE]
    
    result1<-purrr::map2(result1[1:length(Tm)],Tm,function(x,y)cbind(x,data.frame(y)))
    result1<-purrr::map(result1[1:length(Tm)],function(x) x %>% dplyr::rename("Tm"="y"))
    result1<-purrr::map2(result1[1:length(Tm)],AUC[1:length(Tm),],function(x,y)cbind(x,data.frame(y)))
    Pred1<-purrr::map(result1[1:length(Tm)],function(x) x %>% dplyr::rename("AUC"="y"))
    

    #find common IDs between vehicle and treated with converging models
    ID<-data.table::rbindlist(result)$uniqueID %>% as.data.frame(.) %>% unique(.)
    names(ID)<-"uniqueID"
    ID1<-data.table::rbindlist(result1)$uniqueID %>% as.data.frame(.) %>% unique(.)
    names(ID1)<-"uniqueID"
    IID<-intersect(ID,ID1)#get intersecting IDs which have convergence info
  
    #Filter data based on common IDs
   
    Pred<- Pred %>% purrr::keep(.,function(x) x$uniqueID[1] %in% IID$uniqueID)
    Pred1<- Pred1 %>% purrr::keep(.,function(x) x$uniqueID[1] %in% IID$uniqueID)
    #Get differences in Tm, AUC and sum of RSS values for vehicle and treated
    P<-data.frame(data.table::rbindlist(Pred))
    P<- P %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(id = dplyr::n())
    
    P1<-data.frame(data.table::rbindlist(Pred1))
    P1<- P1 %>% dplyr::group_by(uniqueID) %>% dplyr::mutate(id = dplyr::n())
    
    
    Pred<-purrr::map2(Pred,Pred1,function(x,y)x %>% dplyr::mutate(Tm = round(y$Tm-x$Tm,1),AUC = round(y$AUC-x$AUC,2),RSS=round(sum(x$RSS+y$RSS),1)))
    
    #Check sigmoidal fit
    P<- ggplot2::ggplot(Pred[[i]] %>% as.data.frame(), ggplot2::aes(x =CP,y = I,color=dataset)) +
      ggplot2::geom_point(Pred[[i]] %>% as.data.frame(),mapping=ggplot2::aes(x=CP,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred[[i]],mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) 
     
    P1<-P+
      ggplot2::geom_point(Pred1[[i]] %>% as.data.frame(),mapping=ggplot2::aes(x=CP,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred1[[i]],mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
      ggplot2::labs(y = "Relative Solubility",
                    x = "Temperature")+
      ggplot2::annotate("text", x=62, y=0.9, label=  paste("\u0394", "AUC = ",Pred[[i]]$AUC))+ggplot2::annotate("text", x=62, y=0.8, label= paste("\u0394","Tm = ",Pred[[i]]$Tm,"\u00B0C"))
    
    print(P1)
    #Venn Creation
    #find common IDs in trilinear vs sigmoidal
    TRvsSig<-intersect(res[[1]]$uniqueID,IID$uniqueID)
    #find common IDs trilinear vs spline
    TRvsSpl<-intersect(res[[1]]$uniqueID,res_sp[[1]])
    TRvsSp_n<-sum(res[[1]]$uniqueID %in% res_sp[[1]])
    #individual proteins found in trilinear
    TRvsSp_no<-sum(!res[[1]]$uniqueID %in% res_sp[[1]])
    #Individual proteins found in spline
    SpvsTR_no<-sum(!res_sp[[1]] %in% res[[1]]$uniqueID)
    
    #Get protein IDs obtained by each model
    #trilinear
    TRP<-df1$uniqueID[!df1$uniqueID==TRvsSig]
    #sigmoidal
    SIGP<-IID$uniqueID[!IID$uniqueID==TRvsSig]
    ## lty - outline of cirlces
    
    ## fill - colour
    
    ## alpha - colour transparency
    
    grid.newpage()
    draw.single.venn(22, category = "Dog People", lty = "blank", fill = "cornflower blue", 
                     alpha = 0.5)
  }
  
  
        