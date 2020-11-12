# # install after setting up renv
# install.packages("minpack.lm")
# install.packages("rlist")
# install.packages("data.table")
# install.packages("knitr")
# install.packages("ggthemes")
# install.packages("gridExtra")
# install.packages("grid")
# install.packages("readxl")
# install.packages("nls2")
# install.packages("stats")
# install.packages("ggplot2")
# install.packages("pkgcond")
# install.packages("rlist")
# install.packages("pracma")
# install.packages("fs")
# install.packages("tidyverse")
# install.packages("splines")
# install.packages("mgcv")
# install.packages("purrr")
# install.packages("nlstools")
# install.packages("readxl")
# install.packages("nls2")
# install.packages("stringr")
# install.packages("stringi")
# install.packages("mice")

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
library(stringr)
library(stringi)
library(mice)
library(dplyr)
library(ComplexUpset)

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
read_cetsa <- function(f,PSM=FALSE,Batch=TRUE){
  #if the input is peptides, f is a list and PSM = TRUE otherwise f is a data frame and PSM = false
  if (isTRUE(PSM)){
    file.list<-f
    i=1
    
    for ( i in seq(file.list)){
      df<-list()
      df.raw <- readxl::read_excel(file.list[i])
      df[[i]] <- df.raw%>%
        dplyr::select(names(df.raw)[stringr::str_detect(names(df.raw),'Master')],
                      tidyselect::starts_with('Abundance')|tidyselect::starts_with('1'),
                      tidyselect::starts_with('Annotated'),
                      tidyselect::starts_with('Isolation'),
                      tidyselect::starts_with('Ion'),
                      tidyselect::starts_with('Charge'),
                      tidyselect::starts_with('XCorr')) %>%
        dplyr::rename("Accession"="Master Protein Accessions") 
      
      names(df[[i]])<-names(df[[i]]) %>% 
        stringr::str_replace_all(" ","_")
      
      df[[i]]<-df[[i]]<-df[[i]]%>% 
        gather('Accession','value',names(df[[i]])[str_detect(names(df[[i]]),'1')]) %>% 
        dplyr::mutate(sample_id = as.factor(paste("F",as.character(i),sep="")),
                      temp_ref = stringr::str_extract(Accession, '([:digit:][:digit:][:digit:]+[N|C]?)'))
      
      
    }
    return(dplyr::bind_rows(df))
  }else{
    df.raw <- readxl::read_excel(f)
    df <- df.raw%>%
      dplyr::select(Accession, tidyselect::starts_with('Abundance')) %>%
      tidyr::gather('id', 'value', -Accession,description) %>%
      dplyr::mutate(sample = stringr::str_extract(id, '(?<=Abundance: ).*(?=:)')) %>%
      dplyr::mutate(temp_ref = stringr::str_extract(id, '([:digit:][:digit:][:digit:]+[N|C]?)'))
    df
  }
  if (!isTRUE(Batch)){
    df.raw <- readxl::read_excel(file.list[i])
    df[[i]] <- df.raw %>% 
      dplyr::select(Accession, tidyselect::starts_with('Abundance')) %>%
      tidyr::gather('id', 'value', -Accession,description) %>%
      dplyr::mutate(sample_id = as.factor(paste("F",as.character(i),sep="")),
                    temp_ref = stringr::str_extract(Accession, '([:digit:][:digit:][:digit:]+[N|C]?)'))
    df
  }
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
clean_cetsa <- function(df, temperatures = NULL, samples = NULL,carrier=FALSE,PSM=FALSE) {
  samples$sample_id<-as.factor(samples$sample_id)
  samples$sample_name<-as.factor(samples$sample_name)
  if (isTRUE(carrier)){
    r<-df %>% dplyr::filter(!temp_ref=="131C")%>% dplyr::mutate(rank=dplyr::ntile(value,3)) 
    if (isTRUE(PSM)){
      df<-df %>% dplyr::full_join(r, by= c("Accession","Isolation_Interference_in_Percent","value","sample_id","temp_ref","Charge","XCorr","Ions_Matched","Ion_Inject_Time_in_ms","Annotated_Sequence"))
      
    }else{
      df<-df %>% dplyr::full_join(r, by= c("Accession","id","value","sample_id","temp_ref","missing"))
      
    }
    
  }else{
    df<-df %>% dplyr::mutate(rank=dplyr::ntile(value,3))
  }
  if (is.null(temperatures)) warning('No temperature data')
  if (!is.null(samples)) {
    df <- df %>%
      dplyr::left_join(samples, by = c('sample_id' = 'sample_id'))
  }
  if (!is.null(temperatures)) {
    df <- df %>%
      dplyr::left_join(temperatures, by = 'temp_ref')
    
  } else {
    df <- df %>%
      dplyr::rename(temperature = temp_ref)
  }
  if (any(is.na(df$value))){
    df<-df %>% dplyr::rename("sample"="sample_id")
    df<-df %>% dplyr::group_split(Accession, sample) 
    df<-lapply(df,function(x) x%>% dplyr::mutate(missing=sum(is.na(.$value))))
    df<-dplyr::bind_rows(df)
  }
  
  df <- df %>%
    dplyr::select(Accession, sample, temperature,value,missing,rank,sample_name) %>%
    dplyr::filter(!is.na(temperature),!is.na(value)) %>%
    dplyr::group_split(Accession,sample) 
  df<-lapply(df,function(x) x%>%
               dplyr::mutate(value = .$value /.$ value[temperature == min(.$temperature)]) %>% unique(.))
  df<-dplyr::bind_rows(df)
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
 
  df.jointP <- suppressWarnings(df %>%
    dplyr::group_by(Accession,sample) %>% dplyr::mutate(n=dplyr::n()) %>% #copies the number of temperatures per group
    dplyr::filter(n>=10) %>% dplyr::mutate(.,T7 = mean(value[temperature == temperatures[7]]/value[temperature == temperatures[1]]),
                                           T9 = mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]]),
                                           T10 = mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]])) %>%
    dplyr::filter(T7 >= 0.4 & T7 <= 0.6 & T9 < 0.3 & T10 < 0.2))#normalization from TPP


  ## split[[i]] by sample group and filter
  l.bytype <- split.data.frame(df.jointP, df.jointP$sample)
  
  ## determine which sample (F1 through FN) contains the greatest number of curves and use this for normalization
  n.filter <- lapply(l.bytype, nrow)
  df.normP <- l.bytype[[which.max(n.filter)]]
  norm.accessions <- df.normP$Accession
  
  ## calculate median for each sample group
  
  df.mynormset <- df %>% base::subset(Accession %in% norm.accessions)
  
  df.median <- df %>%
    dplyr::group_by(sample,temperature) %>%
    dplyr::mutate(value = median(value,na.rm=TRUE))
  
  
  ## fit curves to the median data for each sample (F1 through FN)
  df.fit <- df.median %>%
    dplyr::group_by(sample) %>% 
    dplyr::filter(temperature<68) %>% 
    dplyr::do(fit = cetsa_fit(d = ., norm = FALSE)) %>% 
    dplyr::mutate(fitted_values = ifelse(!is.logical(fit),list(data.frame(fitted_values=predict(fit))),NA),
                  temperature=ifelse(!is.logical(fit),list(data.frame(temperature=fit$data$t)),NA)) %>% 
    dplyr::select(sample,fitted_values,temperature) 
  
  
  
  ## calculate the fitted values
  d<-length(df.fit$fit)
  #unnest fitted values from list and name value column
  check<-df.fit %>% unnest(c(fitted_values,temperature))
  
  check<-check %>% unique(.) 
  check$sample<-as.factor(check$sample)

  test<-df.median %>% dplyr::group_by(sample,temperature) %>% dplyr::full_join(check,c('sample','temperature'))
  ## calculate ratios between the fitted curves and the median values
  df.out <- test %>%
    dplyr::mutate(correction = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
    dplyr::select('sample','temperature','value','fitted_values','correction')
  df.out<-df.out %>% dplyr::select(-fitted_values,-value)
  ## apply normalization factor to data
 
  df$correction<-df.out$correction
  df <- df %>% 
    dplyr::mutate(norm_value = ifelse(is.na(correction),value,value * correction)) %>% dplyr::ungroup(.)
  
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
#READ PSM FILES
f<-"~/Cliff_MD/Old data"
#read new data
setwd(f)

file.list <- list.files(pattern='*.xlsx')
s1 <- str_replace(file.list, ".xlsx","")#remove Excel extension
s1 <- str_replace_all(s1, "-","")#remove punct
s1<-str_replace_all(s1,"[:upper:][:digit:][:digit:][:digit:]","")
s1<-str_split(s1,"[:digit:][:digit:][:digit:][:digit:]",simplify=TRUE)
s1<-s1[,ncol(s1)]
s1<-str_replace_all(s1,"[:digit:][:digit:]","")
samples <- do.call(paste0, expand.grid(factor(c('F'),levels = c('F')),1:length(s1)))
df.samples <- data.frame(sample_id = samples, sample_name = as.factor(s1))
#check that sample_names match to file list
chk<-mapply(grepl, s1, file.list)
df_raw <- read_cetsa(file.list,PSM=TRUE,Batch=TRUE)



##################################

# #READ PROTEIN DATA not processed in batch
setwd("~/CS7290-/internal data")

samples <- do.call(paste0, expand.grid(factor(c('F'),levels = c('F')),1:24))
df.t <- function(n){
  TMT<-data.frame(NA)
  if(n==10){
    TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67), stringsAsFactors = FALSE)
  }else{
    TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,68), stringsAsFactors = FALSE)
  }
  return(TMT)
}

df.temps<-df.t(11)
#df.samples <- data.frame(sample_id = samples, sample_name = c('MEK_1','MEK_2', 'MEK_3','DMSO_1','DMSO_3'), stringsAsFactors = FALSE)
df.samples <- data.frame(sample_id = samples, sample_name = c(do.call(paste0, expand.grid(factor(c('CFAIMS'),levels = c('CFAIMS')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('CFAIMSPHI'),levels = c('CFAIMSPHI')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('CNOFAIMS'),levels = c('CNOFAIMS')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('CNOFAIMSPHI'),levels = c('CNOFAIMSPHI')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('NCFAIMS'),levels = c('NCFAIMS')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('NCFAIMSPHI'),levels = c('NCFAIMSPHI')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('NCNOFAIMS'),levels = c('NCNOFAIMS')),1:3)),
                                                              do.call(paste0, expand.grid(factor(c('NCNOFAIMSPHI'),levels = c('NCNOFAIMSPHI')),1:3))))                                                              



f<-"~/CS7290-/internal data/ALLP.xlsx"
#read new data
setwd("~/Cliff_MD")
f<-"~/Cliff_MD"
file.list <- list.files(pattern='*.xlsx')
s1 <- str_replace(file.list, ".xlsx","")#remove Excel extension
s1 <- str_replace_all(s1, "-","")#remove punct
s1<-str_replace_all(s1,"[:upper:][:digit:][:digit:][:digit:]","")
s1<-str_split(s1,"[:digit:][:digit:][:digit:][:digit:]",simplify=TRUE)
s1<-s1[,ncol(s1)]
s1<-str_replace_all(s1,"[:digit:][:digit:]","")
samples <- do.call(paste0, expand.grid(factor(c('F'),levels = c('F')),1:length(s1)))
df.samples <- data.frame(sample_id = samples, sample_name = as.factor(s1))
#check that sample_names match to file list
chk<-mapply(grepl, s1, file.list)
df_raw <- read_cetsa(file.list,PSM=TRUE,Batch = FALSE)


#record missing values
missing<-sum(is.na(df_raw$value))
#annotate protein data with missing values
MID<-df_raw[is.na(df_raw$value),]
#flag missing values in data
df_raw<-df_raw %>%  dplyr::mutate(missing = ifelse(is.na(value),1,0))
#assign TMT channel and temperature data
df_clean <- clean_cetsa(df_raw, temperatures = df.temps, samples = df.samples,carrier=TRUE,PSM=TRUE)#assigns temperature and replicate values



df_norm <- normalize_cetsa(df_clean, df.temps$temperature) #normalizes according to Franken et. al. without R-squared filter
df_norm <- df_norm %>% dplyr::select(-value,-correction)

#plot the targets to check values
PLN<-ggplot2::ggplot(df_norm %>% dplyr::filter(Accession=="Q02750"), ggplot2::aes(x = temperature,y = norm_value,color=sample_name)) +
ggplot2::geom_point(ggplot2::aes(x=temperature,y=norm_value))+ ggplot2::ggtitle(paste("Q02750"))+
ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ylim(0,2)

print(PLN)

PLN<-ggplot2::ggplot(df_norm %>% dplyr::filter(Accession=="P36507"), ggplot2::aes(x = temperature,y = norm_value,color=sample_name)) +
ggplot2::geom_point(ggplot2::aes(x=temperature,y=norm_value))+ ggplot2::ggtitle(paste("P36507"))+
ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ylim(0,2)

print(PLN)
################
df_norm_l<- df_norm %>% dplyr::group_split(Accession)
#High, med and low abundance
df_norm <- normalize_cetsa(df_clean, df.temps$temperature) #normalizes according to Franken et. al. without R-squared filter
df_norm_l<- df_norm %>% dplyr::group_split(Accession)

HighP<-df_norm_l %>% purrr::keep(function(x) getmode(x$rank)==3)
MedP<-df_norm_l %>% purrr::keep(function(x) mode(x$rank)==2)
LowP<-df_norm_l %>% purrr::keep(function(x) mode(x$rank)==1)


df_raw<-df_raw %>% dplyr::left_join(df.temps,by="temp_ref")
df_raw<-df_raw %>% dplyr::left_join(df.samples, by ="sample_id")
df_<-df_raw %>% dplyr::filter(temperature<68)
 
upMV <- function(df_) {
  d<-df_ #%>% dplyr::filter(sample_name==unique(.$sample_name)[1])
  d$Accession<-as.factor(d$Accession)
  d$sample_name<-as.factor(d$sample_name)
  d$missing<-as.integer(d$missing)
  d<-d %>% dplyr::filter(missing==1)
   #%>% dplyr::group_split(sample_name)
  d1<-d %>% dplyr::select(-value,-id,-temp_ref)
  d<-spread(d1,temperature,missing,fill=FALSE)
  
  
  d1<-d %>% dplyr::select(-Accession,-sample_id)
  
  p<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[1]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
    base_annotations=list(
      'Missing Values'=intersection_size(
        counts=TRUE,
        aes=aes(fill=names(d1$sample_name))
      )
    ),
    set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[1])
  p2<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[2]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
           base_annotations=list(
             'Missing Values'=intersection_size(
               counts=TRUE,
               aes=aes(fill=names(d1$sample_name))
             )
           ),
           set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[2])
  p3<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[3]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
           base_annotations=list(
             'Missing Values'=intersection_size(
               counts=TRUE,
               aes=aes(fill=names(d1$sample_name))
             )
           ),
           set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[3])
  p4<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[4]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
            base_annotations=list(
              'Missing Values'=intersection_size(
                counts=TRUE,
                aes=aes(fill=names(d1$sample_name))
              )
            ),
            set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[4])
  p5<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[5]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
            base_annotations=list(
              'Missing Values'=intersection_size(
                counts=TRUE,
                aes=aes(fill=names(d1$sample_name))
              )
            ),
            set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[5])
  p6<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[6]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
            base_annotations=list(
              'Missing Values'=intersection_size(
                counts=TRUE,
                aes=aes(fill=names(d1$sample_name))
              )
            ),
            set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[6])
  p7<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[7]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
            base_annotations=list(
              'Missing Values'=intersection_size(
                counts=TRUE,
                aes=aes(fill=names(d1$sample_name))
              )
            ),
            set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[7])
  p8<-upset(d1 %>% dplyr::filter(sample_name==unique(d1)$sample_name[8]),names(d1)[which(!names(d1) %in% "sample_name")],name="Shared Missing Values",
            base_annotations=list(
              'Missing Values'=intersection_size(
                counts=TRUE,
                aes=aes(fill=names(d1$sample_name))
              )
            ),
            set_sizes=FALSE
  )+ggtitle(unique(d1$sample_name)[8])
  #combine plots
  pdf("UpsetMV.pdf", height=10, width=20)
  p
  p2
  p3
  p4
  p5
  p6
  p7
  p8
  dev.off()
}
upMV(df_)


summaryMV <- function(df_) {
  d<-df_ #%>% dplyr::filter(sample_name==unique(.$sample_name)[1])
  d$Accession<-as.factor(d$Accession)
  d$sample_name<-as.factor(d$sample_name)
  d$missing<-as.integer(d$missing)
  d<-d %>% dplyr::filter(missing==1)
  #%>% dplyr::group_split(sample_name)
  d2<-d
  d<-spread(d,sample_name,missing,fill=FALSE)
  d1<-d %>% dplyr::select(-value,-id,-temp_ref,-temperature,-Accession,-sample_id)
  
  
  p<-upset(d1,names(d1),name="Shared Missing Values",
           base_annotations=list(
             'Missing Values'=intersection_size(
               counts=TRUE)),
           set_sizes=FALSE
  )+ggtitle("Summary of Missing Values per Condition")
 
  pdf("UpsetSummary.pdf", height=10, width=20)
  p
  
  dev.off()
}
summaryMV(df_)


plotVenn <- function(a, ...) {
  grid.newpage()
  if (nrow(a) == 2) {
    out <- draw.single.venn(a$missing[1],a$sample_name[1],cex=3, lty = rep("blank",2), fill = c("pink","light blue"), alpha = rep(0.5, 2))          
    grid.arrange(gTree(children=out), top="Missing Values")
  }
  if (nrow(a) == 3) {
    out <- draw.pairwise.venn(a$missing[1], a$missing[2], a$missing[3],a$sample_name[1:2],cex=c(3,3,3), lty = rep("blank",2), fill = c("pink","light blue"), alpha = rep(0.5, 2))          
    grid.arrange(gTree(children=out), top="Missing Values")
    }
  if (nrow(a) == 5) {
    out <- draw.triple.venn(a$missing[1], a$missing[2],a$missing[3],a$missing[4],a$missing[5], a$sample_name[1:3], 
                            lty = rep("blank",3), fill = c("pink","light blue","mediumorchid"),cex=rep(2,5), alpha = rep(0.5, 3))
    grid.arrange(gTree(children=out), top="Missing Values")
  }
  if (length(a) == 15) {
    out <- draw.quad.venn(a$missing[1], a$missing[2],a$missing[3],a$missing[4],a$missing[5],a$missing[6],a$missing[7],a$missing[8],a$missing[9],a$missing[10],a$missing[11],a$missing[12],a$missing[13],a$missing[14],a$missing[15], a$sample_name[1:4], 
                            lty = rep("blank",4), fill = c("skyblue", "pink1", 
                                                           "mediumorchid", "orange"),cex=rep(1,15), alpha = rep(0.5, 4))
    grid.arrange(gTree(children=out), top="Missing Values")
  }
  if (!exists("out")) 
    out <- "Please check your missing data"
  return(out)
}
draw.pairwise.venn(sum(df_raw$missing==TRUE), 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
                                                                                   2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                        0), cat.dist = rep(0.025, 2))
ggplot(df_raw,aes(y=missing,x=sample_name,fill=sample_name,alpha=0.2))+
  facet_grid(~temperature)+
  geom_hist(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
  geom_boxplot(width=0.1) +
  ggplot2::ylab("missing")+
  ggplot2::xlab("sample")+
  theme(axis.text.x = element_text(angle = 90))




#top 5 for each class
HighP<-df_norm_l %>% purrr::keep(function(x) any(x$rank==3,na.rm=TRUE)) %>% head(.,5)
MedP<-df_norm_l %>% purrr::keep(function(x) all(x$rank==2 ,na.rm=TRUE))%>% head(.,5)
LowP<-df_norm_l %>% purrr::keep(function(x) all(x$rank==1,na.rm=TRUE))%>% head(.,5)
#
PLN<-ggplot2::ggplot(HighP[[3]], ggplot2::aes(x = temperature,y = norm_value,color=sample_name)) +
ggplot2::geom_point(ggplot2::aes(x=temperature,y=norm_value))+ ggplot2::ggtitle(paste(HighP[[3]]$Accession[1]))+
ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")

print(PLN)

PLN<-ggplot2::ggplot(MedP[[5]], ggplot2::aes(x = temperature,y = norm_value,color=sample_name)) +
  ggplot2::geom_point(ggplot2::aes(x=temperature,y=norm_value))+ ggplot2::ggtitle(paste(MedP[[5]]$Accession[1]))+
  ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")

print(PLN)

PLN<-ggplot2::ggplot(LowP[[3]], ggplot2::aes(x = temperature,y = norm_value,color=sample_name)) +
  ggplot2::geom_point(ggplot2::aes(x=temperature,y=norm_value))+ ggplot2::ggtitle(paste(LowP[[3]]$Accession[1]))+
  ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")

print(PLN)

#######
#rename columns
df_norm <- df_norm %>% dplyr::rename(dataset = "sample")
colnames(df_norm)<-c("uniqueID","dataset","C","missing","rank","sample_name","I")# uniqueID is protein accession, dataset is sample from PD output,C is temperature, I is intensity

d<-df_norm %>% dplyr::ungroup(.)#get unique data
#d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)#concentration values are defined in uM
#Cliff append names
df.samples<-df.samples %>% dplyr::rename("dataset"="sample_id")
d<-d %>% dplyr::left_join(df.samples,by="dataset")
#remove 67 temp
#d<-d %>% dplyr::filter(C<67)

#create dataset column
CFAIMSPHI <- d %>%
  dplyr::mutate(Dataset = .$sample_name)


CFAIMSPHI<-str_replace(d$sample_name,"[:digit:]","")
#Set unique names as dataset for plots
d$Dataset<-CFAIMSPHI
#
d<-d %>% dplyr::filter(Dataset == "CNOFAIMSPHI"| Dataset=="CNOFAIMS")
#d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)
#study the effect of FAIMS with
d$dataset<-ifelse(stringr::str_detect(tolower(d$Dataset),'cnofaimsphi')=="TRUE" ,"treated","vehicle")

#d$dataset<-ifelse(stringr::str_detect(tolower(d$sample_name),'cfaimsphi')=="TRUE" | stringr::str_detect(tolower(d$sample_name),'ncnofaims')=="TRUE","vehicle","treated")
#d$dataset<-ifelse(d$dataset=="F4" | d$dataset=="F5","vehicle","treated")#dataset is defined by PD sample definition
#study the effect of FAIMS as the dataset
d$CC<-ifelse(stringr::str_detect(tolower(d$sample_name),'cnofaimsphi')=="TRUE" ,1,0)

DF<-d %>% dplyr::group_split(uniqueID) #split null dataset only by protein ID
d_<-d %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,dataset) #split vehicle dataset
d_1<-d %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,dataset) #split treated dataset

#remove values where the lowest temperature has rank 3 (highest intensity)
DF<-DF %>% purrr::keep(function(x) any(x$rank==3))
# #keep data with less than 2 missing values
# DF<-DF %>% purrr::keep(function(x) length(unique(x$C))>=9)#keep values with at least 9 temperature channels (TMT10)
# d_<-d_ %>% purrr::keep(function(x) length(unique(x$C))>=9)
# d_1<-d_1 %>% purrr::keep(function(x) length(unique(x$C))>=9)
#convert to data frame for uniqueID presence
DF<-dplyr::bind_rows(DF)
d_<-dplyr::bind_rows(d_)
d_1<-dplyr::bind_rows(d_1)#1 
#make sure uniqueIDs are present for treated and vehicle
CID<-intersect(d_1$uniqueID,d_$uniqueID)
CID<-intersect(CID,DF$uniqueID)
DF<-DF %>% subset(.$uniqueID %in% CID)
DF$LineRegion<-1
d_<-d_%>% subset(.$uniqueID %in% CID)
d_$LineRegion<-1
d_1<-d_1%>% subset(.$uniqueID %in% CID)
d_1$LineRegion<-1
#split dataset into equal-sized lists
DF<-DF %>%  dplyr::group_split(uniqueID) 
d_<-d_ %>% dplyr::group_split(uniqueID,dataset) 
d_1<-d_1 %>% dplyr::group_split(uniqueID,dataset) 

simNA<-function(df_norm,df.samples,pct = 0.1){
  #####################Start simulations
  #introduce mv
  #generate pattern to impute I only
  #turn uniqueID's into nrows
  df_norm1<-df_norm %>% dplyr::group_split(uniqueID,dataset)
  df_norm1<-df_norm1 %>% purrr::keep(function(x) nrow(x) > 5)
  df_norm1<-purrr::map2(df_norm1,seq(length(df_norm1)),function(x,y) x %>% dplyr::mutate(n=y))
  df_nor<-lapply(df_norm1,function(x) x %>% dplyr::select(-uniqueID,-dataset,-missing,-rank))
  pattern<-lapply(df_norm1,function(x) data.frame(n=rep(0,length(unique(x$dataset))),C=rep(0,length(unique(x$dataset))),I=rep(1,length(unique(x$dataset)))))
  
  Chk<-purrr::map2(df_nor,pattern,function(x,y)mice::ampute(x, prop = 0.1, patterns = y, freq = NULL, mech = "MAR",
              weights = NULL, cont = TRUE, type = NULL, odds = NULL,
              bycases = TRUE, run = TRUE))
  
  Chk1<-df_norm1
  
  Chk<-purrr::map2(Chk1,Chk,function(x,y) x %>% dplyr::full_join(y$amp,by=c("C","n")) %>% dplyr::mutate(I=I.y) %>% dplyr::select(-n,-I.x,-I.y))
  Chk<-dplyr::bind_rows(Chk)
  d<-Chk %>% dplyr::ungroup(.)#get unique data
  #d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)#concentration values are defined in uM
  #Cliff append names
  #df.samples<-df.samples %>% dplyr::rename("dataset"="sample_id")
  d<-d %>% dplyr::left_join(df.samples,by="dataset")
  #remove 67 temp
  #d<-d %>% dplyr::filter(C<67)
  #create dataset column
  CFAIMSPHI <- d %>%
    dplyr::mutate(Dataset = .$sample_name)
  
  
  CFAIMSPHI<-str_replace(d$sample_name,"[:digit:]","")
  #Set unique names as dataset for plots
  d$Dataset<-CFAIMSPHI
  #
  d<-d %>% dplyr::filter(Dataset == "NCNOFAIMSPHI"| Dataset=="NCNOFAIMS")
  #d$CC<-ifelse(d$dataset=="F4" | d$dataset=="F5",0,1)
  #study the effect of FAIMS with
  d$dataset<-ifelse(stringr::str_detect(tolower(d$Dataset),'ncnofaimsphi')=="TRUE" ,"treated","vehicle")
  
  #d$dataset<-ifelse(stringr::str_detect(tolower(d$sample_name),'cfaimsphi')=="TRUE" | stringr::str_detect(tolower(d$sample_name),'ncnofaims')=="TRUE","vehicle","treated")
  #d$dataset<-ifelse(d$dataset=="F4" | d$dataset=="F5","vehicle","treated")#dataset is defined by PD sample definition
  #study the effect of FAIMS as the dataset
  d$CC<-ifelse(stringr::str_detect(tolower(d$sample_name),'ncnofaimsphi')=="TRUE" ,1,0)
  
  DF<-d %>% dplyr::group_split(uniqueID) #split null dataset only by protein ID
  d_<-d %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,dataset) #split vehicle dataset
  d_1<-d %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,dataset) #split treated dataset
  
  
  # #keep data with less than 2 missing values
  # DF<-DF %>% purrr::keep(function(x) length(unique(x$C))>=9)#keep values with at least 9 temperature channels (TMT10)
  # d_<-d_ %>% purrr::keep(function(x) length(unique(x$C))>=9)
  # d_1<-d_1 %>% purrr::keep(function(x) length(unique(x$C))>=9)
  #convert to data frame for uniqueID presence
  DF<-dplyr::bind_rows(DF)
  d_<-dplyr::bind_rows(d_)
  d_1<-dplyr::bind_rows(d_1)#1 
  #make sure uniqueIDs are present for treated and vehicle
  CID<-intersect(DF$uniqueID,d_$uniqueID)
  CID<-intersect(CID,d_1$uniqueID)
  DF<-DF %>% subset(.$uniqueID %in% CID)
  DF$LineRegion<-1
  d_<-d_%>% subset(.$uniqueID %in% CID)
  d_$LineRegion<-1
  d_1<-d_1%>% subset(.$uniqueID %in% CID)
  d_1$LineRegion<-1
  #split dataset into equal-sized lists
  DF<-DF %>%  dplyr::group_split(uniqueID) 
  d_<-d_ %>% dplyr::group_split(uniqueID,dataset) 
  d_1<-d_1 %>% dplyr::group_split(uniqueID,dataset) 
  
  L<-list(DF,d_,d_1)
  return(L)
}

test<-simNA(df_norm,df.samples,0.1)
saveRDS(test,'NAdata10.R')
DF<-test[[1]]
d_<-test[[2]]
d_1<-test[[3]]
################################
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

  d<-lapply(d,function(x) x %>% dplyr::mutate(I=as.numeric(I)))
  df_1<-lapply(d, function(x) {x %>%
      dplyr::group_by(C) %>%
      dplyr::mutate(missing = sum(is.na(.$I)),
                    I=mean(I,na.rm=TRUE)) %>% 
      dplyr::ungroup(.)
    
  })
  #rank intensity values using 3 regions,  rename column as LineRegion
  LR<-lapply(df_1, function(x) {dplyr::ntile(dplyr::desc(x$I),3)%>%
      as.data.frame(.) %>% dplyr::rename("LineRegion"=".")})
  df_1<-purrr::map(df_1,function(x){x %>% dplyr::select(-LineRegion)})#remove Line Region column from one dataset before merging
  
  #Add LR to the list
  df_1<-purrr::map2(df_1,LR, function(x,y) {c(x,y) %>% as.data.frame(.)})
  df_1 <-lapply(df_1,function(x){x %>% dplyr::mutate(C = C,I=I,CC=as.factor(CC))})
  
  #separate by Line Regions
  df1<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==1) %>% as.data.frame(.)})
  df2<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==2) %>% as.data.frame(.)})
  df3<-lapply(df_1,function(x){x %>% dplyr::filter(LineRegion==3) %>% as.data.frame(.)})
  
  #preallocate model data per line region
  LM1<-list(NA)
  LM2<-list(NA)
  LM3<-list(NA)
  df1<-lapply(df1,function(x) x[order(x$C),])
  df2<-lapply(df2,function(x) x[order(x$C),])
  df3<-lapply(df3,function(x) x[order(x$C),])
  # #Flag NA values
  df1<-lapply(df1,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df2<-lapply(df2,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df3<-lapply(df3,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  # #remove NA values
  df1<-lapply(df1,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df2<-lapply(df2,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df3<-lapply(df3,function(x) x %>% dplyr::filter(!is.na(x$I)))
  # #remove empty rows for proteins
  df1<-df1 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  df2<-df2 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  df3<-df3 %>% purrr::keep(function(x) as.logical(nrow(x)>0))
  #get common uniqueIDs
  d1<-dplyr::intersect(dplyr::bind_rows(df1)$uniqueID,dplyr::bind_rows(df2)$uniqueID)
  CID<-dplyr::intersect(d1,dplyr::bind_rows(df3)$uniqueID)
  #keep common uniqueIDs
  
  df1<-df1 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  df2<-df2 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  df3<-df3 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
  #find fitted curves L3<-purrr::map(df3,function(x) tryCatch(lm(formula = I~C,data = x ,na.action=na.omit), error = function(e){NA}))
  L1<-purrr::map(df1,function(x) tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA}))
  LM1<-purrr::map2(df1,L1,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L2<-purrr::map(df2,function(x) tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA}))
  LM2<-purrr::map2(df2,L2,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L3<-purrr::map2(df3,seq(df3),function(x,y) tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'),error=function(e)print(y)))
  LM3<-purrr::map2(df3,L3,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))

  #linear fit per line region
  LM1<-purrr::map2(df1,L1,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM2<-purrr::map2(df2,L2,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM3<-purrr::map2(df3,L3,function(x,y)x %>% dplyr::mutate(M1 = list(y)))


  #fitted curves
  x1<-lapply(LM1, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x2<-lapply(LM2, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x3<-lapply(LM3, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  
  #fit per line region with confidence intervals
  fit1<-purrr::map(LM1,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  fit2<-purrr::map(LM2,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  fit3<-purrr::map(LM3,function(x)x %>% dplyr::mutate(LM1= list(try(predict(x$M1[[1]],se.fit = TRUE)))))
  
  #keep last value for CI 
  fit1 <- purrr::map(fit1,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))
  fit2 <- purrr::map(fit2,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))
  fit3 <- purrr::map(fit3,function(x) x %>% dplyr::mutate(CI=try(tail(x$LM1[[1]]$se.fit,1))))

  #append # of fitted curves to original data (columns must have the same rows for map2)
  df1<-purrr::map2(df1,x1,function(x,y)x %>% dplyr::mutate(fitn=y))
  df2<-purrr::map2(df2,x2,function(x,y)x %>% dplyr::mutate(fitn=y))
  df3<-purrr::map2(df3,x3,function(x,y)x %>% dplyr::mutate(fitn=y))
  
  #append # of fitted curves to original data (columns must have the same rows for map2)
  df1<-purrr::map2(df1,fit1,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  df2<-purrr::map2(df2,fit2,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  df3<-purrr::map2(df3,fit3,function(x,y)x %>% dplyr::mutate(CI=y$CI))
  
  #Reassign line Regions if intensity falls within previous Line Region's CI
  
  df2<-purrr::map2(df1,df2,function(x,y)y %>%
                     dplyr::mutate(LineRegion=ifelse(any(y$I<tail(x$I-x$CI,1)),2,1))) 
 
  df3<-purrr::map2(df2,df3,function(x,y)y %>% 
                     dplyr::mutate(LineRegion=ifelse(any(y$I<tail(x$I-x$CI,1)),3,2))) 
  
  df1<-df1 %>% dplyr::bind_rows(.)
  df2<-df2 %>% dplyr::bind_rows(.)
  df3<-df3 %>% dplyr::bind_rows(.)
  #merge all prepared lists to one data frame
  df_0<-rbind(df1,df2,df3) 
  
  #define line Region as a factor
  df_0$LineRegion<-as.factor(df_0$LineRegion)
  df_0<-df_0 %>% dplyr::group_split(uniqueID)
  return(df_0)
}


#this function takes original data with replicates as an input
CP<-function(df_0,d){ #df_0 is the result data frame and d is the orginal data with replicates
  df_0<-lapply(df_0, function(x) x %>% dplyr::select(-missing,-CI))

  #remove data points with missing data for replicates
  # d<-d %>% purrr::keep(function(x){
  #   nrow(x)>=20
  # })
  #make sure uniqueIDs are consistent among data frames
  d<-dplyr::bind_rows(d)
  df_0<-dplyr::bind_rows(df_0)
  #keep the IDs in df_0 which are present in d
  df_0<-df_0 %>% dplyr::filter(uniqueID %in% d$uniqueID)
  
  #Split into lists once again
  d<-d %>% dplyr::group_split(uniqueID)
  
  #remove IDs that are not common in both datasets
  d<-d %>% purrr::keep(function(x) x$uniqueID[1] %in% df_0$uniqueID)
  #split into list
  df_0<-df_0 %>% dplyr::group_split(uniqueID)
  #For the original data (unlabeled LR) define LR with intensities
  df_0<-suppressWarnings(purrr::map2(d,df_0,function(x,y) x %>% #if the intensity in DF is greater than the max(LR2) label 1 else if the intensity is less than min (LR=2)label 3
                   dplyr::mutate(LineRegion=as.numeric(ifelse(x$I>=min(y$I[y$LineRegion==1]),1,ifelse(x$I<min(y$I[y$LineRegion==2]),3,2))))))
  
  df_n<-vector(mode = "list", length(df_0))
  ctest<-df_n
  dap<-data.frame()
  Split<-df_n
  #This function is to verify consistent line Region assignments for C (temperature) across replicates
  df_0<-lapply(df_0,function(x)x %>% dplyr::arrange(C) %>% dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup())
  
  #subset of the data with shared line region values, using purrr::map to keep size constant
  Split<-purrr::map(df_0,function(x)x %>%subset(n<max(n)) %>% data.frame(.) %>% dplyr::mutate(LineRegion=as.numeric((.$LineRegion))))
  #remove NA values
  Split<-dplyr::bind_rows(Split)
  df_0<-dplyr::bind_rows(df_0)
  #Join df_0 with the subset of values
  dap<-df_0 %>% dplyr::left_join(Split,by=c("C"="C","uniqueID"="uniqueID","dataset"="dataset","I"="I","Dataset"="Dataset","CC"="CC","n"="n"))
  dap<-dap %>% dplyr::group_split(uniqueID)
  dap<-lapply(dap,function(x)x %>% dplyr::mutate(LineRegion=ifelse(.$C<=tail(.$C[.$LineRegion.x==1],1),1,ifelse(.$C<=tail(.$C[.$LineRegion.x==2],1),2,3))) %>% dplyr::select(-LineRegion.x,-LineRegion.y,-sample_name.y,-missing.y))

  return(dap)
}

#preallocate list
results<-vector(mode = "list", length(d_))
results_t<-vector(mode = "list",length(d_1))
results_n<-vector(mode = "list",length(DF))
#enable folder for memoise to save outputs
#gdc <- memoise::cache_filesystem("~/Google Drive/.rcache")#initialize folder location to store results 
#and remove from local memory
results<-suppressWarnings(DLR(d_))#First guess at line regions
results_t<-suppressWarnings(DLR(d_1))
results_n<-suppressWarnings(DLR(DF))

#reassign change points shared between line regions
d<-DF
df_<-lapply(d,function(x) x %>% dplyr::filter(dataset=="vehicle"))
df_1<-lapply(d,function(x) x %>% dplyr::filter(dataset=="treated"))

df_<-suppressWarnings(CP(results,df_))
df_1<-suppressWarnings(CP(results_t,df_1))
DFN<-suppressWarnings(CP(results_n,d))# 
#remove results to save space 
rm(results,results_t,results_n,d_,d_1,DF)#10

#boostrap function
#returns bootstrapped intensities
#n is the sampling window
#N is the number of iterations
# BStrap<-function(Data0,n,N){
#   n<-n
#   N<-N
#   BS<-NA
#   BA<-NA
#   BSvar<-list(NA)
#   BSvar[[1]]<-data.frame()
#   Data0$I<-as.numeric(as.vector(Data0$I))
#   Data0<-Data0 %>% group_split(uniqueID)
#   Data1<-lapply(Data0,function(x) x %>% dplyr::group_by(C) %>% 
#                   dplyr::mutate(I=mean(I),LineRegion=LineRegion)%>% unique(.))
#   
#   BS<-lapply(Data1, function(x){x %>% boot::boot(x$I,statistic = bs,R=n,formula = )})
#   #Bootstrap"sample with replacement"
#   BS1<-lapply(Data0, function(x)x %>% dplyr::select(uniqueID,C,I,LineRegion) %>%  dplyr::sample_n(n,replace=TRUE))
#   #generate mean intensities per $C value
#   BS<-lapply(BS,function(x) x[order(x$C),])
#   
#   return(BS)
# }
# 
# 
# BSvarN<- suppressWarnings(BStrap(DFN,1000,1000))
# BSvar<-suppressWarnings(BStrap(df_,1000,1000))
# BSvar1<-suppressWarnings(BStrap(df_1,1000,1000))
# 
# BSvar<-lapply(BSvar,function(x)x[order(x$C),])
# BSvar1<-lapply(BSvar1,function(x)x[order(x$C),])
# BSvarN<-lapply(BSvarN,function(x)x[order(x$C),])
#append names to list 

df_<-lapply(df_ ,function(x)x[order(x$C),])
df_1<-lapply(df_1,function(x)x[order(x$C),])
DFN<-lapply(DFN,function(x)x[order(x$C),])

df_<-purrr::map2(df_,seq(df_),function(x,y)x %>% dplyr::mutate(N=y))
df_1<-purrr::map2(df_1,seq(df_1),function(x,y)x %>% dplyr::mutate(N=y))
DFN<-purrr::map2(DFN,seq(DFN),function(x,y)x %>% dplyr::mutate(N=y))
tlstat<-function(DF,df,df1,PI=FALSE){
  #preallocate list results
  #control/vehicle
  i<-1
  
  
  DF<-DF
  df<-df
  
  #convert to df for numeric variables C and I 
  DF<-dplyr::bind_rows(DF)
  df1<-dplyr::bind_rows(df1)
  df<-dplyr::bind_rows(df)
  #convert factor to numeric columns
  df1$C<-as.numeric(as.vector(df1$C))
  df$C<-as.numeric(as.vector(df$C))
  DF$C<-as.numeric(as.vector(DF$C))
  
  df1$I<-as.numeric(as.vector(df1$I))
  df$I<-as.numeric(as.vector(df$I))
  DF$I<-as.numeric(as.vector(DF$I))
  #convert back to list
  df1<- df1 %>% dplyr::group_split(uniqueID)
  df<-df %>% dplyr::group_split(uniqueID)
  DF<-DF %>% dplyr::group_split(uniqueID)
  
  if(!isTRUE(PI)){
    mean1<-list()
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- lapply(df,function(x) x %>% data.frame(.) %>% 
                     dplyr::group_nest(LineRegion,uniqueID) %>%
                     dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                   CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                   Tm=with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                   slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                   intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                   rss=purrr::map(M1,function(x){deviance(x)}),
                                   Rsq=purrr::map(M1,function(x){summary(x)$r.squared}), 
                                   dataset="vehicle",
                                   uniqueID=x$uniqueID[1]))
    
    mean1<-lapply(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- lapply(df1,function(x) x %>% data.frame(.) %>%
                       dplyr::group_nest(LineRegion,uniqueID) %>% 
                       dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                     CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                     Tm=with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                     slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                     intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                     rss=map(M1,function(x){deviance(x)}),
                                     Rsq=map(M1,function(x){summary(x)$r.squared}), 
                                     dataset="treated",
                                     uniqueID=x$uniqueID[1]))
    
    
    
    mean1_1<-lapply(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))

    
    mean3<- lapply(DF,function(x) x %>% data.frame(.) %>%
                     dplyr::group_nest(LineRegion,uniqueID) %>% 
                     dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                   CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                   Tm=with(x, stats::approx( x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                   slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                   intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                   rss=map(M1,function(x){deviance(x)}),
                                   Rsq=map(M1,function(x){summary(x)$r.squared}),
                                   dataset="null",
                                   uniqueID=x$uniqueID[1]))
    
    
    mean3<-lapply(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))

     #Apply lax Rsq and negative slope filter to remove flat melt curves
    mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    
    mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
    mean3<-mean3  %>% subset(uniqueID %in% CID)
    #split into lists by uniqueID
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    #Calculate rss0 and rss1 null vs alt
    rss0<-lapply(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
    rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
                                                            Tm = y$Tm[[1]]-x$Tm[[1]]))
    #params for null and alternative models
    pN<-lapply(mean3,function(x)x %>% dplyr::summarise(pN = 4))
    pA<-lapply(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
    
    #sum residuals
    n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = as.numeric(nrow(dplyr::bind_rows(x$data))) + as.numeric(nrow(dplyr::bind_rows(y$data)))))
    #degrees of freedom before
    d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
    d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA)) 
    #delta RSS
    rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS-y$RSS %>% data.frame(.))
    #bind rows
    rssDiff<-dplyr::bind_rows(rssDiff)$.
    rss0<-dplyr::bind_rows(rss0)$RSS
    rss1<-dplyr::bind_rows(rss1)$RSS
    d2<-dplyr::bind_rows(d2)$d2
    d1<-dplyr::bind_rows(d1)$d1
    #F-test
    Fvals<-(rssDiff/rss1)*(d2/d1)
    #append results to data
    ResF<-purrr::map2(mean1,Fvals,function(x,y) x %>% dplyr::mutate(Fvals=y))
    ResF<-purrr::map2(ResF,rss0,function(x,y) x %>% dplyr::mutate(rss0=y))
    ResF<-purrr::map2(ResF,rss1,function(x,y) x %>% dplyr::mutate(rss1=y))
    ResF<-purrr::map2(ResF,rssDiff,function(x,y) x %>% dplyr::mutate(rssDiff=y))
    ResF<-purrr::map2(ResF,d1,function(x,y) x %>% dplyr::mutate(d1=y))
    ResF<-purrr::map2(ResF,d2,function(x,y) x %>% dplyr::mutate(d2=y))
   
    #convert to df
    mean1<-dplyr::bind_rows(mean1)
    ResF<-dplyr::bind_rows(ResF)
   
    #convert results to list
    testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-dataset)
    testResults<-testResults%>% dplyr::left_join(ResF,by="uniqueID")

    #p-val
    testResults<-testResults %>%
      dplyr::mutate(pV = 1-pf(testResults$Fvals,df1=testResults$d1,df2=testResults$d2))
    testResults<-testResults %>% dplyr::mutate(pAdj = p.adjust(.$pV,method="BH"))
    
    #V is zero, so it would not work as a scaling factor
    ggplot(testResults)+
      geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
      geom_line(aes(x=Fvals,y= df(Fvals,df1=4,df2=8)),color="darkred",size = 1.5) +
      theme_bw() +
      coord_cartesian(xlim=c(0,10))+
      ggplot2::xlab("F-values")
    #testResults<-testResults %>% dplyr::filter(pAdj<0.05)
    #scale variables
    M<-median(testResults$rssDiff,na.rm=TRUE)
    V<-mad(testResults$rssDiff,na.rm=TRUE)
    #alternative scaling factor sig0-sq
    altScale<-0.5*V/M
    #filter out negative delta rss
    testResults<-testResults %>% dplyr::filter(rssDiff>0)
    #effective degrees of freedom
    ed1<-MASS::fitdistr(x=testResults$rssDiff, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    ed2<-MASS::fitdistr(x=testResults$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    #scale data
    testScaled <-testResults %>% 
      dplyr::mutate(rssDiff = .$rssDiff/altScale,
                    rss1 =.$rss1/altScale,
                    d1=ed1,
                    d2=ed2)
    #
    #new F-test
    testScaled<-testScaled %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
    Fvals<-testScaled$Fvals
    d1<-testScaled$d1
    d2<-testScaled$d2
   
    #scaled values 
    ggplot(testScaled)+
      geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
      geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
      theme_bw() +
      coord_cartesian(xlim=c(0,10))+
      ggplot2::xlab("F-values")
      #Define checked as filtered protein IDs
    check<-testScaled$uniqueID
    test<-testScaled %>% dplyr::filter(.$pAdj<0.05)
    ggplot(test)+
      geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
      geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
      theme_bw() +
      coord_cartesian(xlim=c(0,10))+
      ggplot2::xlab("F-values")
    
    mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
    mean3<-dplyr::bind_rows(mean3)
    mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    
  } else if (isTRUE(PI)){
    mean1<-list()
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- lapply(df,function(x) x %>% data.frame(.) %>% 
                     dplyr::group_nest(LineRegion,uniqueID) %>%
                     dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                   CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                   Tm=with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                   slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                   intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                   rss=map(M1,function(x){deviance(x)}),
                                   Rsq=map(M1,function(x){summary(x)$r.squared}),
                                   dataset="vehicle",
                                   uniqueID=x$uniqueID[1]))
    
    
    mean1<-lapply(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- lapply(df1,function(x) x %>% data.frame(.) %>%
                       dplyr::group_nest(LineRegion,uniqueID) %>% 
                       dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                     CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                     Tm=with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                     slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                     intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                     rss=map(M1,function(x){deviance(x)}),
                                     Rsq=map(M1,function(x){summary(x)$r.squared}),
                                     dataset="treated",
                                     uniqueID=x$uniqueID[1]))
    
    
    mean1_1<-lapply(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))

    
    mean3<- lapply(DF,function(x) x %>% data.frame(.) %>%
                     dplyr::group_nest(LineRegion,uniqueID) %>% 
                     dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                   CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                   Tm=with(x, stats::approx( x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y,
                                   slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                   intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                   rss=map(M1,function(x){deviance(x)}),
                                   Rsq=map(M1,function(x){summary(x)$r.squared}),
                                   dataset="null",
                                   uniqueID=x$uniqueID[1]))
    
    
    mean3<-lapply(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    
    
    
    
    #Apply lax Rsq and negative slope filter to remove flat melt curves
    mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
    
    mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(mean1)
    mean1_1<-dplyr::bind_rows(mean1_1)
    mean3<-dplyr::bind_rows(mean3)
    
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
    mean3<-mean3  %>% subset(uniqueID %in% CID)
    #split into lists by uniqueID
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    # #Calculate rss0 and rss1 null vs alt
    # rss0<-lapply(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
    # rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
    #                                                         Tm = y$Tm[[1]]-x$Tm[[1]]))
    # #params for null and alternative models
    # pN<-lapply(mean3,function(x)x %>% dplyr::summarise(pN = 4))
    # pA<-lapply(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
    # #sum residuals
    # n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = as.numeric(nrow(dplyr::bind_rows(x$data))) + as.numeric(nrow(dplyr::bind_rows(y$data)))))
    # #degrees of freedom before
    # d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
    # d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA)) 
    # #delta RSS
    # rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS-y$RSS %>% data.frame(.))
    # #bind rows
    # rssDiff<-dplyr::bind_rows(rssDiff)$.
    # rss0<-dplyr::bind_rows(rss0)$RSS
    # rss1<-dplyr::bind_rows(rss1)$RSS
    # d2<-dplyr::bind_rows(d2)$d2
    # d1<-dplyr::bind_rows(d1)$d1
    # #F-test
    # Fvals<-(rssDiff/rss1)*(d2/d1)
    # #append results to data
    # ResF<-purrr::map2(mean1,Fvals,function(x,y) x %>% dplyr::mutate(Fvals=y))
    # ResF<-purrr::map2(ResF,rss0,function(x,y) x %>% dplyr::mutate(rss0=y))
    # ResF<-purrr::map2(ResF,rss1,function(x,y) x %>% dplyr::mutate(rss1=y))
    # ResF<-purrr::map2(ResF,rssDiff,function(x,y) x %>% dplyr::mutate(rssDiff=y))
    # ResF<-purrr::map2(ResF,d1,function(x,y) x %>% dplyr::mutate(d1=y))
    # ResF<-purrr::map2(ResF,d2,function(x,y) x %>% dplyr::mutate(d2=y))
    # 
    # #convert to df
    # mean1<-dplyr::bind_rows(mean1)
    # ResF<-dplyr::bind_rows(ResF)
    # #convert results to list
    # testResults<-mean1 %>% dplyr::select(-slope,-data,-intercept,-LineRegion,-M1,-CI,-Tm,-rss,-Rsq,-AUC,-dataset)
    # testResults<-testResults%>% dplyr::left_join(ResF,by="uniqueID")
    # 
    # #p-val
    # testResults<-testResults %>%
    #   dplyr::mutate(pV = 1-pf(testResults$Fvals,df1=testResults$d1,df2=testResults$d2))
    # testResults<-testResults %>% dplyr::mutate(pAdj = p.adjust(.$pV,method="BH"))
    # 
    # #V is zero, so it would not work as a scaling factor
    # ggplot(testResults)+
    #   geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
    #   geom_line(aes(x=Fvals,y= df(Fvals,df1=4,df2=8)),color="darkred",size = 1.5) +
    #   theme_bw() +
    #   coord_cartesian(xlim=c(0,100))+
    #   ggplot2::xlab("F-values")
    # #scale variables
    # M<-median(testResults$rssDiff,na.rm=TRUE)
    # V<-mad(testResults$rssDiff,na.rm=TRUE)
    # #alternative scaling factor sig0-sq
    # altScale<-0.5*V/M
    # #filter out negative delta rss
    # testResults<-testResults %>% dplyr::filter(rssDiff>0)
    # #effective degrees of freedom
    # ed1<-MASS::fitdistr(x=testResults$rssDiff, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # ed2<-MASS::fitdistr(x=testResults$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
    # #scale data
    # testScaled <-testResults %>% 
    #   dplyr::mutate(rssDiff = .$rssDiff/altScale,
    #                 rss1 =.$rss1/altScale,
    #                 d1=ed1,
    #                 d2=ed2)
    # #
    # #new F-test
    # testScaled<-testScaled %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
    # Fvals<-testScaled$Fvals
    # d1<-testScaled$d1
    # d2<-testScaled$d2
    # 
    # #scaled values 
    # ggplot(testScaled)+
    #   geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
    #   geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
    #   theme_bw() +
    #   coord_cartesian(xlim=c(0,10))+
    #   ggplot2::xlab("F-values")
    # #Define checked as filtered protein IDs
    # check<-testScaled$uniqueID
    # test<-testScaled %>% dplyr::filter(.$pAdj<0.01)
    # ggplot(test)+
    #   geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) + 
    #   geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
    #   theme_bw() +
    #   coord_cartesian(xlim=c(0,10))+
    #   ggplot2::xlab("F-values")
    # 
    # mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
    # mean1_1<-dplyr::bind_rows(mean1_1)
    # mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
    # mean3<-dplyr::bind_rows(mean3)
    # mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
     results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    # 
  }
  return(results)
}


spstat<-function(DF,df,df1,PI=FALSE,Ftest=TRUE){
  
  df<-df %>% purrr::keep(function(x) is.data.frame(x))
  df1<-df1 %>% purrr::keep(function(x) is.data.frame(x))
  DF<-DF %>% purrr::keep(function(x) is.data.frame(x))
  
  df<-dplyr::bind_rows(df)
  df1<-dplyr::bind_rows(df1)
  DF<-dplyr::bind_rows(DF)
  
  #plot spline results
  
  df1$I<-as.numeric(as.vector(df1$I))
  df$I<-as.numeric(as.vector(df$I))
  DF$I<-as.numeric(as.vector(DF$I))
  #mutate to get CV values
  DF<-DF %>% dplyr::group_split(C,uniqueID) 
  DF<- lapply(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df<-df %>% dplyr::group_split(C,uniqueID,dataset) 
  df<- lapply(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df1<-df1 %>% dplyr::group_split(C,uniqueID,dataset) 
  df1<- lapply(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  
  #convert to data frame
  
  df<-dplyr::bind_rows(df)
  df1<-dplyr::bind_rows(df1)
  DF<-dplyr::bind_rows(DF)
  
  #aggregate column 
  #switch from factor to numeric
  #convert factor to numeric columns
  df1$C<-as.numeric(as.vector(df1$C))
  df$C<-as.numeric(as.vector(df$C))
  DF$C<-as.numeric(as.vector(DF$C))
  
  df1$I<-as.numeric(as.vector(df1$I))
  df$I<-as.numeric(as.vector(df$I))
  DF$I<-as.numeric(as.vector(DF$I))
  #convert back to list
  DF<-DF %>% dplyr::group_split(uniqueID) 
  df<-df %>% dplyr::group_split(uniqueID) 
  df1<-df1 %>% dplyr::group_split(uniqueID) 
  #remove NA vaues
  DF<-lapply(DF,function(x)na.omit(x))
  df<-lapply(df,function(x)na.omit(x))
  df1<-lapply(df1,function(x)na.omit(x))
  if(!isTRUE(PI)){
    
    #alternative spline fit method : Generalized Additive Models
    #fit penalized splines
    m <- lapply(df,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(x$I ~ s(x$C,k=5), data = x , method = "REML")))))
    m<-m %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m<-lapply(m,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                               sum = list(summary(.$M1[[1]])),
                                               Tm=ifelse(any(class(dplyr::first(.$M1[[1]]))=="gam"),try(with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                               rss=deviance(.$M1[[1]]),
                                               CV_pct = .$CV_pct,
                                               AUC = pracma::trapz(.$M1[[1]]$fit),
                                               rsq=summary(x$M1[[1]])$r.sq))
    #m<-lapply(m,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    m1 <- lapply(df1,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data = x , method = "REML")))))
    m1<-m1 %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m1<-lapply(m1,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                  sum = list(summary(.$M1[[1]])),
                                                  Tm=ifelse(any(class(dplyr::first(.$M1[[1]]))=="gam"),try(with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                                  rss=deviance(.$M1[[1]]),
                                                  CV_pct = .$CV_pct,
                                                  AUC = pracma::trapz(.$M1[[1]]$fit),
                                                  rsq=summary(x$M1[[1]])$r.sq))
    #m1<-lapply(df1,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    
    
    mn<- lapply(DF,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data =x, method = "REML")))))
    mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    mn<-lapply(mn,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                 sum = list(summary(.$M1[[1]])),
                                                 Tm=ifelse(any(class(dplyr::first(.$M1[[1]]))=="gam"),try(with(x, stats::approx(x$I,x$C, xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                                 rss=deviance(.$M1[[1]]),
                                                 CV_pct=.$CV_pct,
                                                 AUC = pracma::trapz(.$M1[[1]]$fit),
                                                 rsq=summary(.$M1[[1]])$r.sq))
    #mn<-lapply(mn,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    #convert to df and split by uniqueID 
    mean1<-dplyr::bind_rows(m)
    mean1_1<-dplyr::bind_rows(m1)
    mean3<-dplyr::bind_rows(mn)
    
    #obtain common uniqueIDs
    CID<-intersect(mean1$uniqueID,mean1_1$uniqueID)
    CID<-intersect(CID,mean3$uniqueID)
    #subset common uniqueIDs
    mean1<-mean1 %>% subset(uniqueID %in% CID)
    mean1_1<-mean1_1  %>% subset(uniqueID %in% CID)
    mean3<-mean3  %>% subset(uniqueID %in% CID)
    #split into lists by uniqueID
    mean1<-mean1 %>% dplyr::group_split(uniqueID)
    mean1_1<-mean1_1 %>% dplyr::group_split(uniqueID)
    mean3<-mean3 %>% dplyr::group_split(uniqueID)
    
    #Cliff
    results<-dplyr::bind_rows(mean1,mean1_1,mean3)
    if(isTRUE(Ftest)){
      #Calculate rss0 and rss1 null vs alt
      rss0<-lapply(mean3,function(x)data.frame(RSS = deviance(x$M1[[1]]),
                                               RSS1 = deviance(x$sig[[1]])))
      rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = (deviance(x$M1[[1]])+deviance(y$M1[[1]])),
                                                              RSS1 =(deviance(x$sig[[1]])+deviance(y$sig[[1]])),
                                                              Tm = y$Tm[[1]]-x$Tm[[1]]))
      #params for null and alternative models
      pN<-lapply(mean3,function(x)x %>% dplyr::summarise(pN = x$M1[[1]]$rank,
                                                         n=x$M1[[1]]$df.null,
                                                         pN1 = x$sig[[1]]$rank,
                                                         n1 = x$sig[[1]]$df.null))
      pA<-lapply(mean1_1,function(x)x %>% dplyr::summarise(pA = 2*(x$M1[[1]]$rank),
                                                           pA1= 2*(x$sig[[1]]$rank)))
      #data points
      n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = x$M1[[1]]$df.null + y$M1[[1]]$df.null,
                                                             n1_1 = x$sig[[1]]$df.null + y$sig[[1]]$df.null))
      #degrees of freedom before
      d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
      d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA))
      #DoF 'new_k' < k
      d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA1-y$pN1))
      d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1_1-y$pA1))
      
      #delta RSS
      rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS1-y$RSS1 %>% data.frame(.))
      #bind rows
      rssDiff<-dplyr::bind_rows(rssDiff)$.
      
      rss0<-dplyr::bind_rows(rss0)$RSS1
      rss1<-dplyr::bind_rows(rss1)$RSS1
      d2<-dplyr::bind_rows(d2)$d2
      d1<-dplyr::bind_rows(d1)$d1
      #F-test
      Fvals<-(rssDiff/rss1)*(d2/d1)
      #append results to data
      ResF<-purrr::map2(mean1,Fvals,function(x,y) x %>% dplyr::mutate(Fvals=y))
      ResF<-purrr::map2(ResF,rss0,function(x,y) x %>% dplyr::mutate(rss0=y))
      ResF<-purrr::map2(ResF,rss1,function(x,y) x %>% dplyr::mutate(rss1=y))
      ResF<-purrr::map2(ResF,rssDiff,function(x,y) x %>% dplyr::mutate(rssDiff=y))
      ResF<-purrr::map2(ResF,d1,function(x,y) x %>% dplyr::mutate(d1=y))
      ResF<-purrr::map2(ResF,d2,function(x,y) x %>% dplyr::mutate(d2=y))
      
      #convert to df
      mean1<-dplyr::bind_rows(mean1)
      ResF<-dplyr::bind_rows(ResF)
      #convert results to list
      testResults<-mean1 %>% dplyr::select(-M1,-sig)
      testResults<-testResults%>% dplyr::left_join(ResF,by="uniqueID")
      
      #p-val
      
      testResults<-testResults %>%
        dplyr::mutate(pV = as.numeric(1-pf(Fvals,df1=d1[1],df2=max(d2))))
      testResults<-testResults %>% dplyr::mutate(pAdj = p.adjust(.$pV,method="BH"))
      
      ggplot(testResults)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
      #scale variables
      M<-median(testResults$rssDiff,na.rm=TRUE)
      V<-mad(testResults$rssDiff,na.rm=TRUE)
      #alternative scaling factor sig0-sq
      altScale<-0.5*V/M
      #filter out negative delta rss
      testResults<-testResults %>% dplyr::filter(rssDiff>0)
      #effective degrees of freedom
      ed1<-MASS::fitdistr(x=testResults$rssDiff, densfun = "chi-squared", start = list(df=2))[["estimate"]]
      ed2<-MASS::fitdistr(x=testResults$rss1, densfun = "chi-squared", start = list(df=2))[["estimate"]]
      #scale data
      testScaled <-testResults %>%
        dplyr::mutate(rssDiff = .$rssDiff/altScale,
                      rss1 =.$rss1/altScale,
                      d1=ed1,
                      d2=ed2)
      #
      #new F-test
      testScaled<-testScaled %>% dplyr::mutate(Fvals=(rssDiff/rss1)*(d2/d1))
      Fvals<-testScaled$Fvals
      d1<-testScaled$d1
      d2<-testScaled$d2
      
      #scaled values
      ggplot(testScaled)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,1))+
        ggplot2::xlab("F-values")
      #Define checked as filtered protein IDs
      check<-testScaled$uniqueID
      test<-testScaled %>% dplyr::filter(.$pV<0.10)
      test$d1<-MASS::fitdistr(x=test$rssDiff, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      test$d2<-MASS::fitdistr(x=test$rss1, densfun = "chi-squared", start = list(df=1))[["estimate"]]
      
      ggplot(test)+
        geom_density(aes(x=Fvals),fill = "steelblue",alpha = 0.5) +
        geom_line(aes(x=Fvals,y= df(Fvals,df1=d1,df2=d2)),color="darkred",size = 1.5) +
        theme_bw() +
        coord_cartesian(xlim=c(0,10))+
        ggplot2::xlab("F-values")
      
      mean1<-mean1 %>% dplyr::filter(mean1$uniqueID %in% test$uniqueID)
      mean1_1<-dplyr::bind_rows(mean1_1)
      mean1_1<-mean1_1 %>% dplyr::filter(mean1_1$uniqueID %in% test$uniqueID)
      mean3<-dplyr::bind_rows(mean3)
      mean3<-mean3 %>% dplyr::filter(mean3$uniqueID %in% test$uniqueID)
      results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
      
    }
    
  }
  return(results)
}
#gettrilinear results
tlresults<-list()
tlresults_PI<-list()
#confidence intervals
tlresults<-tlstat(DFN,df_,df_1,PI=TRUE)#place null, vehicle and treated lists with no prediction intervals

tlf<-function(tlresults,DFN,APfilt=FALSE,PF=TRUE){
  ##Apply Filters
  #####################
  if(isTRUE(APfilt)){
    tlresults1<-tlresults#save unfiltered data
    #apply filters prior to hypothesis testing
    tlresults<-tlresults %>% keep(function(x) min(as.numeric(x$Rsq),na.rm=TRUE) >= 0.40)
    tlresults<-tlresults %>% keep(function(x) mean(as.numeric(x$slope),na.rm=TRUE) <= -0.02)
    #tlresults<-tlresults %>% keep(function(x)  sum(data.frame(x)[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <10)#move data with extremely large RSS values 
    # tlresults<-tlresults %>% keep(function(x) sum(data.frame(x)[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss'],na.rm=TRUE) <1.3)
    tlresults<-tlresults %>% keep(function(x) sum(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),na.rm=TRUE) > sum(unlist(x[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),na.rm=TRUE))#remove data with extremely large RSS values 
    tlresults<-tlresults %>% keep(function(x) mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'Tm']),na.rm=TRUE) < mean(unlist(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'Tm']),na.rm=TRUE))
    #tlresults<-tlresults %>% keep(function(x) max(data.frame(x)$slope[x$LineRegion==2],na.rm=TRUE) < -0.03)#the linear region have the largest slope < 0.03
    #tlresults<-tlresults %>% keep(function(x) length(x$slope)>8)#remove list values with less than 5 rows
    #tlresults<-tlresults %>% keep(function(x) abs(max(x$slope[!x$LineRegion==2] ,na.rm=TRUE)) < 0.1)#eeps plateau values where the min abs(slope) < 0.06
    #steepest slope in vehicle and treatment has to be less than 0.06C
  }
  Nsum<-list()
  Nsum[[1]]<-data.frame(RSS=0,Tm=0)
  #get the summed rss values for null
  Nsum<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "null")) %>% 
                 dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
  
  #get the summed rss values for vehicle
  Rssv<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% 
                 dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>%head(.,1))
  #get the summed rss values for treated
  Rsst<-lapply(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% 
                 dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
  #find the rss difference between treated and vehicle 
  
  Rssv<-lapply(Rssv,function(x)na.omit(x))
  Rsst<-lapply(Rsst,function(x)na.omit(x))
  #find common IDs
  CID<-intersect(dplyr::bind_rows(Rsst)$uniqueID,dplyr::bind_rows(Rssv)$uniqueID)
  #keep common IDs
  Rssv<-Rssv %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID)) 
  Rsst<-Rsst %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))
  Nsum<-Nsum %>% purrr::keep(function(x) isTRUE(x$uniqueID %in% CID))                           
  K1<-data.frame(dplyr::bind_rows(purrr::map2(Rsst,Rssv,function(x,y) data.frame(RSSd = x$RSS-y$RSS, Tma = x$Tm - y$Tm)))) 
  K2<-data.frame(uniqueID = dplyr::bind_rows(Rssv)$uniqueID)
  Dsum<-data.frame(K1,K2)
  
  Dsum<-Dsum %>% dplyr::mutate(rank = dplyr::ntile(Dsum$Tma,7))
  #keep data where the difference in RSS is less than the null
  #nsum converted to data frame
  Nsum<-data.frame(RSSn=dplyr::bind_rows(Nsum))
  names(Nsum)<-c("RSSn","Tmn","dataset","uniqueID")
  Nsum<-Nsum %>% dplyr::filter(uniqueID %in% CID)
  Nsum<-Nsum %>% dplyr::mutate(id=rownames(Nsum))
  
  Nsum$dataset<-as.factor(Nsum$dataset)
  #mutate data frame
  #join two data frames by uniqueID
  Dsum1<-Dsum %>% dplyr::left_join(Nsum,by = c("uniqueID"="uniqueID"))
  #Childs
  Dsum2<-Dsum %>% dplyr::right_join(Nsum,by = c("uniqueID"="uniqueID"))
  
  Dsum<-Dsum2
  Dsum$RSSd<-Dsum1$RSSd
  Dsum$Tma<-Dsum1$Tma
  Dsum<-Dsum %>% dplyr::mutate(rank = dplyr::ntile(Dsum$Tma,7))
  if (isTRUE(PF)){
    #rank the data by Tm change
    
    #arrange data from greater Tm and RSS difference to lowest
    Dsum<-dplyr::arrange(Dsum, dplyr::desc(Tma), dplyr::desc(RSSd))  %>% dplyr::filter(RSSd>0) 
    
    test<-data.frame()
    test<-Dsum[which(Dsum$RSSn>Dsum$RSSd),] %>% data.frame()#get the stable proteins (+ = Rsstreated-Rssvehicle)
    rssdec<-data.frame()
    rssdec<-data.frame(data.table::fsort(test$RSSd,decreasing=TRUE))#decreasing Rss differences
    names(rssdec)<-"Rssd"
    
    tmdec<-data.frame()
    tmdec<-data.table::fsort(test$Tma,decreasing=TRUE) %>% data.frame()
    names(tmdec)<-"Tm"
    
    test<-tmdec %>% inner_join(test,by=c("Tm"="Tma"))#orders data by decreasing Tm
    orows<-data.frame()
    orows <- test
    orows$id<-sapply(orows$id, function(x) as.numeric(as.character(x)))
    Df1<-tlresults[orows$id] #divide 1=highly destabilized,4=noeffect,7=highly stabilized
  }else{
    tlresults<-dplyr::bind_rows(tlresults)
    
    #order by RSS differences while keeping original rownames for index
    #create an external data frame for stabilized proteins
    Df1<-Dsum %>% dplyr::left_join(tlresults,by=c("uniqueID")) %>% dplyr::rename("dataset"="dataset.y") %>% dplyr::select(-dataset.x)
    Df1<-Df1 %>% dplyr::group_split(uniqueID)

    
  }
  
  df1<-list()
  
  #get uniqueID and dataset for stable proteins with decreasing RSS differences
  df1<-lapply(Df1,function(x) x %>% dplyr::select(uniqueID,dataset) %>% head(.,1))
  df1<-data.frame(dplyr::bind_rows(df1))
  
  #unlist to data.frame
  #order the original data by RSS differences
  #
  DFN<- dplyr::bind_rows(DFN)
  DFN$uniqueID<-as.vector(DFN$uniqueID)
  df1$uniqueID<-as.vector(df1$uniqueID)
  
  
  df2<-df1 %>% dplyr::right_join(DFN,by=c("uniqueID")) %>% dplyr::select(-dataset.x) %>%
    dplyr::rename("dataset"="dataset.y") 
  
  
  return(list(df1,df2,Df1))
}
#return filtered lists
res<-tlf(tlresults,DFN,APfilt=FALSE,PF=TRUE)
tlCI<-function(i,df1,df2,Df1,overlay=TRUE){
  null<-data.frame()
  i<-i
  df1<-df1
  Df1<-Df1[[i]]
  
  DF1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) 
  
  null<-Df1 %>% subset(dataset == "null")
  
  pred1<-predict(null$M1[[1]], interval="confidence") %>% as.data.frame(.)
  pred2<-predict(null$M1[[2]], interval="confidence")%>% as.data.frame(.)
  pred3<-predict(null$M1[[3]], interval="confidence")%>% as.data.frame(.)
  Pred1<-NA
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  
  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred<-dplyr::bind_rows(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred<-dplyr::bind_rows(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
    Pred<-dplyr::bind_rows(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred<-dplyr::bind_rows(pred1,pred3)
  }
  rownames(Pred)<-as.vector(1:nrow(Pred))
  
  #Pred<-Pred[1:length(DF1$C),]##############
  Pred<-cbind(Pred,DF1$C[1:nrow(Pred)],DF1$I[1:nrow(Pred)])################
  names(Pred)<-c("fit","lower","upper","C","I")
  
  Pred$Treatment<-null$dataset[1]##################
  Pred<-na.omit(Pred)
  Pred$C<-as.numeric(as.vector(Pred$C))
  Pred$I<-as.numeric(as.vector(Pred$I))
  PLN<-ggplot2::ggplot(Pred, ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1$uniqueID[1],"null"))+
    ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ 
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
    annotate("text", x=60, y=1,label=paste("RSS= ",round(sum(unlist(null$rss)),3)))

  
  DF_f<-df2 %>%subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated')) %>% subset(dataset=="vehicle")
  
  vehicle<-Df1 %>% subset(dataset == "vehicle")
  
  pred1<-predict(vehicle$M1[[1]], interval="confidence")%>% as.data.frame(.)
  pred2<-predict(vehicle$M1[[2]], interval="confidence")%>% as.data.frame(.)
  pred3<-predict(vehicle$M1[[3]], interval="confidence")%>% as.data.frame(.)
  Pred1<-NA
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  
  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred1<-dplyr::bind_rows(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred1<-dplyr::bind_rows(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
    Pred1<-dplyr::bind_rows(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred1<-dplyr::bind_rows(pred1,pred3)
  }
  
  #Pred<-Pred[1:length(DF1$C),]##############
  Pred1<-data.frame(Pred1,DF_f$C[1:nrow(Pred1)],DF_f$I[1:nrow(Pred1)])################
  names(Pred1)<-c("fit","lower","upper","C","I")
  
  Pred1$Treatment<-vehicle$dataset[1]##################
  Pred1<-na.omit(Pred1)
  rownames(Pred1)<-1:nrow(Pred1)
  Pred1$C<-as.numeric(as.vector(Pred1$C))
  Pred1$I<-as.numeric(as.vector(Pred1$I))
  
  PLR_P1<-ggplot2::ggplot(Pred1, ggplot2::aes(x = C,y = fit,color=Treatment))+ggplot2::geom_point(Pred1, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_ribbon(data=Pred1,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)
  
  DF_f1<-data.frame()
  DF_f1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated'))
  DF_f1<-DF_f1 %>% subset(dataset =="treated")
  
  
  treated<-data.frame()
  treated<-Df1 %>% subset(dataset == "treated")
  
  pred1<-predict(treated$M1[[1]], interval="confidence")%>% as.data.frame(.)
  pred2<-predict(treated$M1[[2]], interval="confidence")%>% as.data.frame(.)
  pred3<-predict(treated$M1[[3]], interval="confidence")%>% as.data.frame(.)
  
  pred1<-na.omit(pred1)
  pred2<-na.omit(pred2)
  pred3<-na.omit(pred3)
  Pred2<-NA
  FIT<- NA
  LOW<-NA
  HI<-NA
  if (nrow(pred1)>0 & nrow(pred2)>0 & nrow(pred3)>0){
    Pred2<-dplyr::bind_rows(pred1,pred2,pred3)
  } else if (nrow(pred2)>0 & nrow(pred3)>0){
    Pred2<-dplyr::bind_rows(pred2,pred3)
  } else if (nrow(pred1)>0 & nrow(pred2)>0){
    Pred2<-dplyr::bind_rows(pred1,pred2)  
  }else if (nrow(pred1)>0 & nrow(pred3)>0){
    Pred2<-dplyr::bind_rows(pred1,pred3)
  }
  rownames(Pred2)<-as.vector(1:nrow(Pred2))
  
  #Pred<-Pred[1:length(DF1$C),]##############
  Pred2<-data.frame(Pred2,DF_f1$C[1:nrow(Pred2)],DF_f1$I[1:nrow(Pred2)])################
  names(Pred2)<-c("fit","lower","upper","C","I")
  
  Pred2$Treatment<-treated$dataset[1]##################
  Pred2<-na.omit(Pred2)
  rownames(Pred2)<-as.vector(1:nrow(Pred2))
  #Area under the curve using trapezoid rule
  P1_AUC <- pracma::trapz(as.numeric(as.vector(Pred1$C)),as.numeric(as.vector(Pred1$I)))
  P2_AUC <- pracma::trapz(as.numeric(as.vector(Pred2$C)),as.numeric(as.vector(Pred2$I)))
  #Residuals
  rn<-data.frame(residuals=c(residuals(null$M1[[1]]),residuals(null$M1[[2]]),residuals(null$M1[[3]])))
  Pred<-cbind(Pred,rn)
  names(Pred)<- names(Pred)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]]),residuals(vehicle$M1[[3]])))
  Pred1<-cbind(Pred1,rn)
  names(Pred1)<- names(Pred1)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]]),residuals(treated$M1[[3]])))
  Pred2<-cbind(Pred2,rn[1:nrow(Pred2),])
  names(Pred2)<- names(Pred2)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  Preds<-rbind(Pred1,Pred2)
  Preds$C<-as.numeric(as.vector(Preds$C))
  Preds$I<-as.numeric(as.vector(Preds$I))
  Pred2$C<-as.numeric(as.vector(Pred2$C))
  Pred2$I<-as.numeric(as.vector(Pred2$I))
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
    Tm_d<-as.numeric(as.vector(Tm2))-as.numeric(as.vector(Tm1))
    p<-expression(paste(Delta, "AUCdiff"))
    if(AUCd>0){
      P1_AUC <- pracma::trapz(as.numeric(as.vector(Pred1$C)),as.numeric(as.vector(Pred1$lower)))
      P2_AUC <- pracma::trapz(as.numeric(as.vector(Pred2$C)),as.numeric(as.vector(Pred2$upper)))
      AUCd<-round(P2_AUC-P1_AUC,2)
    }else{
      P1_AUC <- pracma::trapz(as.numeric(as.vector(Pred1$C)),as.numeric(as.vector(Pred1$upper)))
      P2_AUC <- pracma::trapz(as.numeric(as.vector(Pred2$C)),as.numeric(as.vector(Pred2$lower)))
      AUCd<-round(P2_AUC-P1_AUC,2)
    }
    AUCd<-as.numeric(AUCd)
    PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
      ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Df1$uniqueID[1],"alternative"))+
      ggplot2::annotate("text", x=60, y=1, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "vehicle"),'rss']))+
                                                                                 sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "treated"),'rss'])),3)))+
      ggplot2::annotate("text", x=60, y=0.9, label=  paste("\u0394", "AUC = ",AUCd))+ ggplot2::annotate("text", x=60, y=0.8, label= paste("\u0394","Tm = ",Tm_d,"\u00B0C"))
    #bquote(Value~is~sigma~R^{2}==.(r2.value)))
    par(mfrow=c(2,2))
    print(PLR_P2)
  }else if(overlay=="FALSE"){
    PLR<-PLR_P2+ggplot2::geom_point(data=Pred,mapping=ggplot2::aes(x=C,y=I))+ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ggplot2::ggtitle(paste(Df1$uniqueID[1],"alternative"))+facet_wrap("Treatment") 
    print(PLR)
  }
}
i=1

plotTL<-tlCI(i,res[[1]],res[[2]],res[[3]],overlay=TRUE)
#df1 <- only IDs in order desc(stability)
#df2<-original data in order  
#Df1 <- ordered spline results 
###############################
#get spline results
spresults<-list()
spresults_PI<-list()

spresults<-spstat(DFN,df_,df_1,PI=FALSE,Ftest=FALSE)

spf<-function(spresults,filters = TRUE){
  #spresults<-spresults %>% dplyr::select(-sample_name.y,-Dataset.y) %>% dplyr::rename("Dataset"="Dataset.x","sample_name"="sample_name.x")
  spresults1<-spresults
  spresults<-spresults %>% dplyr::group_split(uniqueID)
  spresults<-spresults %>% purrr::keep(function(x) any(!is.na(unique(x$uniqueID))))
  
  if(!isTRUE(filters))
     {
       sl<-lapply(seq_len(length(spresults)),function(x) as.numeric({paste(x)})) 
       sp<-purrr::map2(spresults,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))  
       sp<-dplyr::bind_rows(sp)  
       df1<-data.frame(uniqueID = unique(sp$uniqueID))  
       df2<-dplyr::bind_rows(DFN)  
       df2$C<-as.numeric(as.vector(df2$C)) 
       df2$I<-as.numeric(as.vector(df2$I))  
       df2<-sp %>% left_join(df2, by = c("uniqueID"="uniqueID","dataset"="dataset","C"="C","I"="I","CC"="CC","sample_name.x"="sample_name.x","LineRegion"="LineRegion","missing.x"="missing.x","Dataset"="Dataset")) 
       Df1<-spresults
  }else{
    #Apply filters 
    #keep the positive AUC differences
    spresults<-spresults %>% keep(function(x) mean(x$AUC[x$dataset=="treated"],na.rm=TRUE)>mean(x$AUC[!x$dataset=="vehicle"],na.rm=TRUE))
    spresults<-spresults %>% keep(function(x) max(x$lambda)<1)
    
    #get Tm and RSS differences
    sp<-lapply(spresults, function(x) x %>% dplyr::mutate(Tmd= x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'Tm'][[1]] - x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'Tm'][[1]],
                                                          RSSd = sum(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']) - sum(x[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),
                                                          AUCd = x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'AUC'])[[1]]- x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'AUC'][[1]])
    #conserve list indexes
    sl<-lapply(seq_along(length(sp)),function(x) as.numeric({paste(x)}))
    
    #insert list index column
    sp<-map2(sp,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))
    sp<-dplyr::bind_rows(sp) +
      sp<-dplyr::arrange(sp,dplyr::desc(AUCd),dplyr::desc(RSSd),dplyr::desc(Tmd)) %>% dplyr::select(uniqueID,id) %>% unique(.) 
    #arrange results by decreasing AUCd, RSSd and Tmd and standardize the order in spresults
    #Df1 holds the model results and stats for splines 
    
    df1<-data.frame(uniqueID = unique(sp$uniqueID))
    df2<-dplyr::bind_rows(DFN) 
    df2$C<-as.numeric(as.vector(df2$C))
    df2$I<-as.numeric(as.vector(df2$I))
    df2<-sp %>% left_join(df2, by = c("uniqueID"="uniqueID","dataset"="dataset","C"="C","I"="I","CC"="CC","sample_name.x"="sample_name.x","LineRegion"="LineRegion","missing.x"="missing.x","Dataset"="Dataset")) 
    Df1<-spresults[sp$id]
  }
  ret<-list()
  ret[[1]]<-df1
  ret[[2]]<-df2
  ret[[3]]<-Df1
  return(ret)
}
res_sp<-spf(spresults,filters=FALSE)

#plot global histograms for variability 
his_sp<-function(Df1,df.temps,MD=FALSE){
  test<-dplyr::bind_rows(Df1)%>% dplyr::select(Dataset,CV_pct,dataset,C) %>% unique(.)
  test<-test %>% dplyr::rename("temperature"="C")
  test<-test %>% dplyr::left_join(df.temps,by="temperature")
  test<-test %>% dplyr::rename("C"="temperature")
  if(isTRUE(MD)){
    test$Dataset<-as.factor(test$Dataset)
    levels(test$Dataset)<-rev(levels(test$Dataset))
  ggplot(test,aes(y=CV_pct,x=Dataset,fill=Dataset,alpha=0.2))+
    facet_grid(~C)+
    geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
    geom_boxplot(width=0.1) +
    ggplot2::ylab("RSD%")+
    ggplot2::xlab("sample")+
    theme(axis.text.x = element_text(angle = 90))
  }else{  
    test$dataset<-as.factor(test$dataset)
    levels(test$dataset)<-rev(levels(test$dataset))
  ggplot(test,aes(y=CV_pct,x=dataset,fill=dataset,alpha=0.2))+
    facet_grid(~temp_ref)+
    geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
    geom_boxplot(width=0.1) +
    ggplot2::ylab("RSD%")+
    ggplot2::xlab("sample")+
    theme(axis.text.x = element_text(angle = 90))
  
    
  }
}
his_sp(res_sp[[3]],df.temps,MD=TRUE)
#plot global violin plots for Rsquared splines
rsq<-function(Df1,MD=FALSE){#input df_raw
  test<-dplyr::bind_rows(Df1) %>% dplyr::select(uniqueID,Dataset,dataset,rsq)
  
  if(MD==TRUE){
  test$Dataset<-as.factor(test$Dataset)
  levels(test$Dataset)<-rev(levels(test$Dataset))
  ggplot(test,aes(y=rsq,x=Dataset,fill=Dataset,alpha=0.2))+
    geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
    geom_boxplot(width=0.1) +
    ggplot2::ylab("R^2 splines")+
    ggplot2::xlab("sample")+
    theme(axis.text.x = element_text(angle = 90))
  }else{
    test$dataset<-as.factor(test$dataset)
    levels(test$dataset)<-rev(levels(test$dataset))
    ggplot(test,aes(y=rsq,x=dataset,fill=dataset,alpha=0.2))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("R^2 splines")+
      ggplot2::xlab("sample")+
      theme(axis.text.x = element_text(angle = 90))
  }
  
}
rsq(res_sp[[3]],MD=TRUE)

#plot global violin plots for Rsquared splines
AUC_V<-function(Df1,MD=FALSE){#input df_raw

  test<-dplyr::bind_rows(Df1) %>% dplyr::group_split(uniqueID,Dataset)
  test<-lapply(test,function(x) x[1,])
  test<-dplyr::bind_rows(test) %>% dplyr::select(uniqueID,Dataset,dataset,AUC) %>% unique(.)
  
  if(MD==TRUE){
    test$Dataset<-as.factor(test$Dataset)
    levels(test$Dataset)<-rev(levels(test$Dataset))
    ggplot(test,aes(y=AUC,x=Dataset,fill=Dataset,alpha=0.2))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("AUC splines")+
      ggplot2::xlab("sample")+
      theme(axis.text.x = element_text(angle = 90))+
      ggplot2::ylim(0,100)
  }else{
    test$dataset<-as.factor(test$dataset)
    levels(test$dataset)<-rev(levels(test$dataset))
    ggplot(test,aes(y=AUC,x=dataset,fill=dataset,alpha=0.2))+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("AUC splines")+
      ggplot2::xlab("sample")+
      theme(axis.text.x = element_text(angle = 90))+
      ggplot2::ylim(0,100) 
  }
  
}
AUC_V(res_sp[[3]],MD=TRUE)
 #plot global histograms for variability 
viol_int<-function(fdata,df_raw,df.samples,df.temps,MD=FALSE){#input df_raw
  #get filtered data 
  
  fdata<-dplyr::bind_rows(fdata) 
  df.samples<-dplyr::rename(df.samples,"Dataset"="sample_name")
  #fdata<-fdata %>% dplyr::rename("Dataset"="sample_name")
  fdata<-fdata %>% #rename sample_names
    dplyr::left_join(df.samples,by="Dataset") %>% 
    dplyr::select(-I)  #dataset.x is vehicle or treated 
  
  #dataset is F1 column
  #put equal column headers before joining
  test<-df_raw %>% dplyr::rename("RunID"="sample") 
  df.samples<-df.samples %>% dplyr::rename("RunID"="dataset")
  
  test<-test %>% dplyr::left_join(df.samples,by="RunID") 
  #filter by sample_name
  test<-test %>% dplyr::mutate(Dataset=str_replace(test$Dataset,"[:digit:]",""))
  test<-test %>% dplyr::filter(Dataset %in% fdata$Dataset)
  #add temperature values
  test<-test %>% 
    dplyr::full_join(df.temps,by="temp_ref") %>% dplyr::rename("C"="temperature","I"="value")
  fdata<-fdata %>% dplyr::rename("dataset"="dataset.x")
  fdata$dataset<-as.factor(fdata$dataset)
  fdata<-test %>% dplyr::full_join(fdata,by=c("Accession"="uniqueID","C"="C","Dataset"="Dataset"))
  # d<-max(nchar(unique(fdata$Dataset)))
  # fdata<-fdata %>% dplyr::mutate(dataset = ifelse(nchar(.$Dataset)<d,"vehicle","treated"))
  # 
  fdata$logI<-log(fdata$I,base=exp(2))
  
  fdata<-fdata %>% dplyr::select(Dataset,C,temp_ref,I,dataset,logI)
  fdata$temp_ref<-as.factor(fdata$temp_ref)
  fdata$C<-as.factor(fdata$C)
  if(isTRUE(MD)){
  fdata$Dataset<-as.factor(fdata$Dataset)
  levels(fdata$Dataset)<-rev(levels(fdata$Dataset))
  ggplot(fdata,aes(y=logI,x=Dataset,fill=Dataset,alpha=0.2))+
    facet_grid(~C)+
    geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
    geom_boxplot(width=0.1) +
    ggplot2::ylab("log2(Intensity)")+
    ggplot2::xlab("sample")+
    theme(axis.text.x = element_text(angle = 90))+
    ggplot2::ylim(3,11) 
  }else{
    fdata$dataset<-as.factor(fdata$dataset)
    levels(fdata$dataset)<-rev(levels(fdata$dataset))
    ggplot(fdata,aes(y=logI,x=dataset,fill=dataset,alpha=0.2))+
      facet_grid(~C)+
      geom_violin(na.rm=TRUE,show.legend="FALSE",color=NA)+theme_bw()+
      geom_boxplot(width=0.1) +
      ggplot2::ylab("log2(Intensity)")+
      ggplot2::xlab("sample")+
      theme(axis.text.x = element_text(angle = 90))+
      ggplot2::ylim(3,11)
  }
}
viol_int(res_sp[[3]],df_raw,df.samples,df.temps,MD=TRUE)
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
  #rename columns 
  df2<-df2 %>% dplyr::mutate(missing=df2$missing.x,sample_name=df2$sample_name.x)
  df2<-df2 %>% dplyr::rename("n"="n.x") %>% dplyr::select(-n.y)
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-as.character(df2$uniqueID)
  vmissing<-length(which(df2 %>% subset(dataset=="vehicle") %>% dplyr::select(missing) == 1))
  tmissing<-length(which(df2 %>% subset(dataset=="treated") %>% dplyr::select(missing) == 1))
  
  #get original data
  ###########################################
  df1<-df1$uniqueID[i]
  DF1<-df2%>% subset(.$uniqueID %in% df1)
  null<-Df1[[i]] %>% subset(.$uniqueID %in% df1 & dataset %in% "null")
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
    BSVarN<-df2 %>% subset(uniqueID == df1 ) %>%dplyr::group_by(C)# %>%  dplyr::mutate(I=mean(I))
    BSVar <-df2 %>% subset(uniqueID == df1 & dataset== "vehicle")%>%dplyr::group_by(C)# %>% dplyr::mutate(I=mean(I))
    BSVar1 <-df2 %>% subset(uniqueID == df1 & dataset== "treated")%>%dplyr::group_by(C)# %>% dplyr::mutate(I=mean(I))
    BSVarN<-BSVarN %>% dplyr::mutate(dataset="null") 
    BSVar<-BSVar %>% dplyr::mutate(dataset="vehicle")
    BSVar1<-BSVar1 %>% dplyr::mutate(dataset="treated")
    BSVarN$dataset<-as.factor(BSVarN$dataset)
    BSVar$dataset<-as.factor(BSVar$dataset)
    BSVar1$dataset<-as.factor(BSVar1$dataset)
    fit <-  stats::smooth.spline(x = BSVar$C, y=BSVar$I,cv=F)
    fit1<-  stats::smooth.spline(x = BSVar1$C, y=BSVar1$I,cv=F)
    fitN<-  stats::smooth.spline(x = BSVarN$C, y=BSVarN$I,cv=F)
    
    #####try GAM
    #fit penalized splines
    m <- mgcv::gam(I ~ s(C,k=5), data = BSVar , method = "ML")
    m1<-  mgcv::gam(I ~ s(C,k=5), data =BSVar1, method = "ML")
    mn<-  mgcv::gam(I ~ s(C,k=5), data = BSVarN, method = "ML")
    
    #####try GAM
    
    #Plot boostrapped  residuals with 95%CI
    #PLP<-plot(m, shade = TRUE, seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
    #generate random values from a multivariate normal distribution
   
    #get some parmeters
    Vb <- vcov(m)
    newd <- with(BSVar, data.frame(C = seq(min(C), max(C), length = 20)))%>% as.data.frame(.)
    pred <- predict(m, newd, se.fit = TRUE)%>% as.data.frame(.)
    se.fit <- pred$se.fit
    #get some parmeters
    Vb1<- vcov(m1) 
    newd1<- with(BSVar1,data.frame(C = seq(min(C), max(C), length = 30)))%>% as.data.frame(.)
    pred1<- predict(m1,newd1,se.fit = TRUE) %>% as.data.frame(.)
    se.fit1<- pred1$se.fit
    #generate std
    set.seed(42)
    N <- 1000
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
      geom_point(BSVar, mapping=ggplot2::aes(x=C,y=I,color = dataset,shape=factor(CC)))+
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
    
    pred1$missing<-sum(df2$missing)
    Pred<-data.frame(m$fitted.values,m$residuals)
    names(Pred)<- c("fit","rn")
    Pred$dataset<-as.factor("vehicle")
    BSVar$dataset<-as.factor("vehicle")
    Pred1<-data.frame(m1$fitted.values,m1$residuals)
    names(Pred1)<-c("fit","rn")
    Pred1$dataset<-as.factor("treated")
    Preds<-rbind(Pred,Pred1)
    BSVar1$dataset<-as.factor("treated")
    PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=dataset)) +ggplot2::geom_point()+ 
      ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1]," ",Df1[[i]]$Dataset[1]))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
    print(PLrs)
    plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
      ggplot2::geom_point(BSvar,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(c(as.character(df1),"alternative"))+
      # ggplot2::annotate("text", x=50, y=1, label= paste("missing values: vehicle",vmissing[1]))+
      # ggplot2::annotate("text", x=50, y=0.9, label= paste("missing values treated",BSvar1$Dataset.x[1],":",tmissing[1]))                    
      ggplot2::annotate("text", x=60, y=1, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])))+
      ggplot2::annotate("text", x=60, y=0.9, label=  paste("\u0394", "AUC = ",round(pred1$AUC[1],3)))+
      ggplot2::annotate("text", x=60, y=0.8, label= paste("\u0394","Tm = ",round(pred1$Tm[1],3),"\u00B0C"))
    
    plot<-plot1+
      ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
      ggplot2::labs(y = "Relative Solubility",
                    x = "Temperature (\u00B0C)")+
      coord_cartesian(ylim = c(-0.1, 1.1),xlim = c(37,67)) 
      
    print(plot)  
  }

#generate 95%CI for splines
Pred<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],overlay=TRUE,alpha=0.05)
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
          control = nls.control(maxiter = 20)),
      silent = TRUE)
}
repeatFits <- function(x,y, start= c(Pl = 0, a = 550, b=10),
                       seed = NULL, alwaysPermute = FALSE, maxAttempts = 20){
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
                       alwaysPermute = FALSE,maxAttempts = 20){
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
computeRSSdiff <- function(x,y,treatment,maxAttempts = 20, repeatsIfNeg = 30){
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

sigfit<-function(DFN,i,MD=TRUE){
    
    df_<-DFN
    df_1<-DFN
    DFN<-DFN
    nlm1<-list(rep(NA,length(df_)))
    nlm2<-nlm1
    dfc<-nlm1
    dfy<-nlm1
    result<-nlm1
    result1<-nlm1
    Pred<-nlm1
    Pred1<-list(rep(NA,1))
    
    df_<-df_ %>% purrr::keep(function(x) !is.null(nrow(x)))
    df_1<-df_1 %>% purrr::keep(function(x) !is.null(nrow(x)))
    DFN<-DFN %>% purrr::keep(function(x) !is.null(nrow(x)))
    #keep data with at least 9 rows
    df_<-df_ %>% purrr::keep(function(x)nrow(x)>9)
    df_1<-df_1 %>% purrr::keep(function(x) nrow(x)>9)
    DFN<-DFN %>% purrr::keep(function(x) nrow(x)>9)
    #convert to data frame
    df_<-dplyr::bind_rows(df_)#vehicle
    df_1<-dplyr::bind_rows(df_1)#treated
    DFN<-dplyr::bind_rows(DFN)#Null
    
    #subset dataset for vehicle and treated
    df_<-df_ %>% dplyr::filter(dataset=="vehicle")
    df_1<-df_1 %>% dplyr::filter(dataset=="treated")
    
    #intersect IDs
    d1<-dplyr::intersect(df_$uniqueID,df_1$uniqueID)
    CID<-dplyr::intersect(d1,DFN$uniqueID)
    #convert back to list
    df_<-df_ %>% dplyr::group_split(uniqueID)
    df_1<-df_1 %>% dplyr::group_split(uniqueID)
    DFN<-DFN %>% dplyr::group_split(uniqueID)
    #keep data with common IDs

    df_<-df_ %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
    df_1<-df_1 %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
    DFN<-DFN %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
    
    #convert to data frame
    df_<-dplyr::bind_rows(df_)#vehicle
    df_1<-dplyr::bind_rows(df_1)#treated
    DFN<-dplyr::bind_rows(DFN)#Null
    
   
    #convert to factor
    df_$C<-as.numeric(as.vector(df_$C))
    df_1$C<-as.numeric(as.vector(df_1$C))
    DFN$C<-as.numeric(as.vector(DFN$C))
    
    df_$I<-as.numeric(as.vector(df_$I))
    df_1$I<-as.numeric(as.vector(df_1$I))
    DFN$I<-as.numeric(as.vector(DFN$I))
    
    #convert back to list
    df_<-df_ %>% dplyr::group_split(uniqueID)
    df_1<-df_1 %>% dplyr::group_split(uniqueID)
    DFN<-DFN %>% dplyr::group_split(uniqueID)
    
    #order by C
    df_<-lapply(df_,function(x)x[order(x$C),])
    df_1<-lapply(df_1,function(x)x[order(x$C),])
    DFN<-lapply(DFN,function(x)x[order(x$C),])
    
    #sigmoidal fit for vehicle
    #remove non-sigmoidal behaving proteins
    nlm1<-df_ %>% purrr::keep(function(x)!inherits(fitSingleSigmoid(x$C,x$I),'try-error')) #flag and omit non-sigmoidal data
   
    #calculate fit for the subset of proteins with sigmoidal behavior
    nlm2<-lapply(nlm1,function(x)x %>% dplyr::mutate(fit = list(try(fitSingleSigmoid(x$C,x$I))))) #fit sigmoids
    
    #remove proteins with a plateau value not =  zero 
    CT<-nlm2 %>% purrr::keep(function(x) !coef(x$fit[[1]])[[1]]==0)#find values where Pl = 0 
    #get confidence intervals for sigmoidal function
    dfc <- purrr::map(CT,function(x) try(nlstools::confint2(x$fit[[1]],level=0.95)))#collect predicted values
    #obtain sigmoidal equation parameters
    nlm2<-purrr::map2(CT,dfc,function(x,y)x %>% dplyr::mutate(Pl=y[1],a=y[2],b=y[3],Pl1=y[4],a1=y[5],b1=y[6]))
    #ready confidence intervals for plot
    result<- purrr::map(nlm2,function(x)x %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                                                 HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                                                 nV=length(predict(x$fit[[1]]))))#get lower CI
    
    AUC<-purrr::map(result,function(x)pracma::trapz(x$C,x$I))
    Tm<- lapply(result,function(x) try(stats::approx(x$I,x$C,0.5)$y[1]))
    TM<- Tm %>% purrr::keep(function(x) !inherits(x,'try-error'))
    #append Tm 
    result<-purrr::map2(result,Tm,function(x,y)x %>% dplyr::mutate(Tm=y))
    #append AUC
    Pred<-purrr::map2(result,AUC,function(x,y)x %>% dplyr::mutate(AUC=y))
    
    #append RSS
    Pred<-purrr::map(Pred,function(x) x %>% dplyr::mutate(RSS=deviance(x$fit[[1]])))
    Pred<-purrr::map2(Pred,seq(Pred),function(x,y)x %>% dplyr::mutate(n=y))
    nlm1<-list()
    CT<-list()
    dfc<-list()
    
    #sigmoidal fit for treated
    nlm1<-df_1 %>% purrr::keep(function(x)!try(inherits(fitSingleSigmoid(x$C,x$I),'try-error'))) #flag and omit non-sigmoidal data
    
    nlm2<-lapply(nlm1,function(x)x %>% dplyr::mutate(fit = list(try(fitSingleSigmoid(x$C,x$I))))) #fit sigmoids
    CT<-nlm2 %>% purrr::keep(function(x) !coef(x$fit[[1]])[[1]]==0)#find values where Pl = 0 
   
    dfc <- purrr::map(CT,function(x) try(nlstools::confint2(x$fit[[1]],level=0.95)))#collect predicted values
    nlm2<-purrr::map2(CT,dfc,function(x,y)x %>% dplyr::mutate(Pl=y[1],a=y[2],b=y[3],Pl1=y[4],a1=y[5],b1=y[6]))
    #calculate CI
    result1<- purrr::map(nlm2,function(x)x %>% dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                                                   HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                                                   nV=length(predict(x$fit[[1]]))))#get lower CI
    AUC<-purrr::map(result1,function(x)pracma::trapz(x$C,x$I))
    Tm<- lapply(result1,function(x) try(stats::approx(x$I,x$C,0.5)$y[1]))
    TM<- Tm %>% purrr::keep(function(x) !inherits(x,'try-error'))
    #append Tm 
    result1<-purrr::map2(result1,Tm,function(x,y)x %>% dplyr::mutate(Tm=y))
    #append AUC
    Pred1<-purrr::map2(result1,AUC,function(x,y)x %>% dplyr::mutate(AUC=y))
    
    #append RSS
    Pred1<-purrr::map(Pred1,function(x) x %>% dplyr::mutate(RSS=deviance(x$fit[[1]])))
    Pred1<-purrr::map2(Pred1,seq(Pred1),function(x,y)x %>% dplyr::mutate(n=y))
    #find common IDs between vehicle and treated with converging models
    ID<-dplyr::bind_rows(result)$uniqueID 
    names(ID)<-"uniqueID"
    ID1<-dplyr::bind_rows(result1)$uniqueID %>% unique(.)
    names(ID1)<-"uniqueID"
    IID<-intersect(ID,ID1)#get intersecting IDs which have convergence info
    O_ID<-setdiff(ID,ID1)#get different proteins where one fit failed
    #Filter data based on common IDs
    Pred<- Pred %>% purrr::keep(.,function(x) x$uniqueID[1] %in% IID)
    Pred1<- Pred1 %>% purrr::keep(.,function(x) x$uniqueID[1] %in% IID)
     #remove fit column
    Pred<-lapply(Pred,function(x) x %>% dplyr::select(-fit))
    Pred1<-lapply(Pred1,function(x) x %>% dplyr::select(-fit))
    
    Pred<-purrr::map2(Pred,seq(Pred),function(x,y)x %>% dplyr::mutate(N=y))
    Pred1<-purrr::map2(Pred1,seq(Pred1),function(x,y)x %>% dplyr::mutate(N=y))
    
    chk<-dplyr::bind_rows(Pred) %>% dplyr::filter(uniqueID %in% "P36507")
    print(chk$N)
    chk<-dplyr::bind_rows(Pred) %>% dplyr::filter(uniqueID %in% "Q02750")
    print(chk$N)
    
    PRED<-Pred
    PRED1<-Pred1
    Pred<-Pred[[i]] %>% dplyr::select(uniqueID,dataset,Dataset,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,
                                      AUC,RSS,LOW,HI,N,Rsq)
    Pred1<-Pred1[[i]]%>% dplyr::select(uniqueID,dataset,Dataset,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,
                                       AUC,RSS,LOW,HI,N,Rsq)
    Pred$LOW<-Pred$LOW[[1]]
    Pred$HI<-Pred$HI[[1]]
    
    Pred1$LOW<-Pred1$LOW[[1]]
    Pred1$HI<-Pred1$HI[[1]]
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    Pred1$dAUC<-as.double(round(Pred1$AUC[1]-Pred$AUC[1],2))
    Pred1$RSS<-round(sum(Pred1$RSS[1]+Pred$RSS[1]),3)
    
    #
    if(isTRUE(MD)){
      #Check sigmoidal fit
      P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=Dataset))+
        ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = Dataset))+
        ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=Dataset), alpha = 0.2 ) +
        ggplot2::annotate("text", x=60, y=0.8, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
        ggplot2::annotate("text", x=60, y=0.7, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1]))+
        ggplot2::annotate("text", x=60, y=0.9, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))
      
      P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = Dataset))+
        ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=Dataset), alpha = 0.2 ) 
      
      
      print(P1)
    }else{
      #Check sigmoidal fit
      P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=dataset))+
        ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
        ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
        ggplot2::annotate("text", x=60, y=0.8, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
        ggplot2::annotate("text", x=60, y=0.7, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
        ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1]))+
        ggplot2::annotate("text", x=60, y=0.9, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))
      
      P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
        ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) 
      
      
      print(P1)
      
      
    }
    return(dplyr::bind_rows(Pred,Pred1))
}
PlSig<-sigfit(DFN,i,MD=TRUE)
VennMet<-function(res,res_sp,MID){
    #Venn Creation
    #find common IDs in trilinear vs sigmoidal
    TRvsSig<-intersect(res$uniqueID,IID$uniqueID)
    #find common IDs trilinear vs spline
    TRvsSpl<-intersect(res$uniqueID,res_sp)
    TRvsSp_n<-sum(res$uniqueID %in% res_sp)
    #individual proteins found in trilinear
    TRvsSp_no<-sum(!res$uniqueID %in% res_sp)
    #Individual proteins found in spline
    SpvsTR_no<-sum(!res_sp %in% res$uniqueID)
    
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
  
  
        
        