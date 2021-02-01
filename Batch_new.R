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
library(DBI)
library(furrr)
library(tibble)
library(ComplexUpset)
library(patchwork)
library(caret)
library(ggpubr)

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
read_cetsa <- function(protein_path,peptide_path,Prot_Pattern,PSM=FALSE,Batch=TRUE){
  file.list<-protein_path
  i=1
  peptide_path<-as.character(peptide_path)
  protein_path<-as.character(protein_path)
  find<-c('[:digit:][:digit:][:digit:][N|C]|[:digit:][:digit:][:digit:]')
  #if the input is peptides, f is a list and PSM = TRUE otherwise f is a data frame and PSM = false
  if (isTRUE(PSM)){
    df<-list()
    df1<-list()
    df.raw<-data.frame()
    setwd(peptide_path)
    f<-list.files(peptide_path,pattern="PSMs")
    for ( i in seq(f)){
      #first get PSM file read
      df.raw <- readxl::read_excel(f[i])
      df[[i]] <- df.raw%>%
        dplyr::select(names(df.raw)[stringr::str_detect(names(df.raw),'Master')],
                      tidyselect::starts_with('Abundance')|tidyselect::starts_with('1'),
                      tidyselect::starts_with('Annotated'),
                      tidyselect::contains('Isolation'),
                      tidyselect::contains('Ion'),
                      tidyselect::contains('Charge'),
                      tidyselect::contains('PEP'),
                      tidyselect::contains('Modifications'),
                      tidyselect::contains('Cleavages'),
                      tidyselect::starts_with('XCorr'),
                      tidyselect::contains('Delta'),
                      tidyselect::contains('File ID'),
                      tidyselect::contains('S/N'),
                      tidyselect::contains('Spectrum')) %>%
        dplyr::rename("Accession"="Master Protein Accessions") %>% 
        dplyr::select(Accession,
                      tidyselect::starts_with('Abundance'),
                      tidyselect::starts_with('Annotated'),
                      tidyselect::contains('Isolation'),
                      tidyselect::contains('Ion'),
                      tidyselect::contains('Charge'),
                      tidyselect::contains('PEP'),
                      tidyselect::contains('Modifications'),
                      tidyselect::contains('Cleavages'),
                      tidyselect::starts_with('XCorr'),
                      tidyselect::contains('Delta'),
                      tidyselect::contains('File ID'),
                      tidyselect::contains('S/N'),
                      tidyselect::contains('Spectrum'),
                      -tidyselect::contains('Grouped')) %>% 
        tidyr::gather('id', 'value', -Accession) %>% 
        dplyr::mutate(temp_ref = stringr::str_extract_all(id,find))
      
    }
    df2<-dplyr::bind_rows(df) %>% tidyr::unnest(temp_ref)
    names(df2 )<-names(df2 ) %>% 
      stringr::str_replace_all(" ","_")
    #then link proteins to peptide sample_id
    setwd(protein_path)
    g<-list.files(protein_path,pattern=as.character(Prot_Pattern))
    df<-list()
    
      for (i in seq(g)){
        df[[i]]<- readxl::read_excel(g[i])
        
        df[[i]]<- df[[i]] %>% 
          dplyr::select(Accession,tidyselect::starts_with('Abundance'),-tidyselect::contains('Grouped')) %>% 
          tidyr::gather('id', 'value', -Accession) %>%
          dplyr::mutate(sample_id = ifelse(!is.na(str_extract(id, str_c('F',"[:digit:][:digit:]"))),str_extract(id, str_c('F',"[:digit:][:digit:]")),stringr::str_extract_all(id,str_c('F','[:digit:]'))),
                        temp_ref = stringr::str_extract_all(id,find),
                        missing=ifelse(is.na(value),1,0),
                        Spectrum_File=g[i])
        class(df[[i]]$sample_id)<-"character"
      }
      df<-dplyr::bind_rows(df) %>% 
        tidyr::unnest(cols = temp_ref) %>%
        dplyr::select(Accession,sample_id,Spectrum_File,value) %>%
        dplyr::rename(Protein_value=value)
     
      df2<-df %>% dplyr::left_join(df2,by=c("Accession","Spectrum_File"))
      
      
      return(df2)
    }else{
    i<-1
    file.list<-g
    df2<-list()
    df<-list()
    df2[[i]]<-data.frame()
    df[[i]]<-data.frame()
    
    for ( i in seq(file.list)){
      df<- readxl::read_excel(file.list[i])
      df2[[i]] <- df %>% 
        dplyr::select(Accession,tidyselect::starts_with('Abundance'),-tidyselect::contains('Grouped')) %>% 
        tidyr::gather('id', 'value', -Accession) %>%
        dplyr::mutate(sample_id = ifelse(!is.na(str_extract(id, str_c('F',"[:digit:][:digit:]"))),str_extract(id, str_c('F',"[:digit:][:digit:]")),stringr::str_extract_all(id,str_c('F','[:digit:]'))),
                      temp_ref = stringr::str_extract_all(id,find),
                      missing=ifelse(is.na(value),1,0))%>% tidyr::unnest(cols = c(sample_id, temp_ref))
    }
    df2<-dplyr::bind_rows(df2)
    
    
    df_raw_D_R<-df2 %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(value,3))%>% dplyr::select(sample_id,Accession,rank)
    df2<-df2 %>% dplyr::left_join(df_raw_D_R,by=c('sample_id','Accession'))
    df2<-df2 %>% dplyr::group_split(Accession,sample_id)
    df2<-purrr::map(df2,function(x){ x %>%  dplyr::mutate(missing_pct=100*sum(is.na(value))/length(value))})
    df2<-dplyr::bind_rows(df2)
    
  }
  
  
  if(isTRUE(Batch)){
    
    i=1
    file.list<-f
    df<-vector("list",length(file.list))
    df[[i]]<-data.frame()
    for( i in seq(file.list)){
      df.raw <- readxl::read_excel(file.list[i])
      
      check<-any(stringr::str_detect(names(df.raw),'Grouped'))
      
      if (isTRUE(check)){
        df[[i]] <- df.raw %>%
          dplyr::select(Accession,tidyselect::starts_with('Abundance'),-tidyselect::contains('Grouped')) %>%
          tidyr::gather('id', 'value', -Accession) %>%
          dplyr::mutate(sample_id = as.factor(paste("F",as.character(i),sep="")),
                        temp_ref = stringr::str_extract_all(id,find),
                        missing=ifelse(is.na(value),1,0))
      }else{
        df[[i]] <- df.raw %>%
          dplyr::select(Accession,tidyselect::starts_with('Abundance')) %>%
          tidyr::gather('id', 'value', -Accession) %>%
          dplyr::mutate(sample_id = ifelse(!is.na(str_extract(id, str_c('F',"[:digit:][:digit:]"))),str_extract(id, str_c('F',"[:digit:][:digit:]")),stringr::str_extract_all(id,str_c('F','[:digit:]'))),
                        temp_ref = as.factor(unlist(stringr::str_extract_all(id,find))),
                        missing=ifelse(is.na(value),1,0))
      }
    }
    df<-lapply(df,function(x) x %>% tidyr::unnest(cols="sample_id"))
    df<-dplyr::bind_rows(df)
    
    df_raw_D_R<-df %>% dplyr::filter(temp_ref=="126") %>% dplyr::mutate(rank=dplyr::ntile(value,3))%>% dplyr::select(sample_id,Accession,rank)
    df<-df %>% dplyr::left_join(df_raw_D_R,by=c('sample_id','Accession'))
    df<-df %>% dplyr::group_split(Accession,sample_id)
    df<-purrr::map(df,function(x){ x %>%  dplyr::mutate(missing_pct=100*sum(is.na(value))/length(value))})
    df<-dplyr::bind_rows(df)
    
    return(df)
    
  }else{
    
    return(df2)
    
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
clean_cetsa <- function(df, temperatures = NULL, samples = NULL,PSM=FALSE) {
  if (is.null(temperatures)) {return(warning('No temperature data'))}
  if (is.null(samples)) {return(warning('No sample data'))}
  
  df <- df %>%
    dplyr::left_join(samples, by = c('sample_id' = 'sample_id'))
  
  
  df <- df %>%
    dplyr::left_join(temperatures, by = 'temp_ref')
 
  if(isTRUE(PSM) & any(str_detect(names(df),"Fraction")=="TRUE")){
    df <- df %>%
      dplyr::group_split(Accession, sample_id, temperature)
    df<-lapply(df,function(x) x%>% dplyr::mutate(value=sum(value)))
    
    df<-dplyr::bind_rows(df)
    
    df <- df %>%
      dplyr::select(Accession, sample_id, temperature,value,rank,Fraction,missing,missing_pct,sample_name) %>% 
      dplyr::filter(!is.na(temperature),!is.na(value)) %>%
      dplyr::group_by(Accession,sample_id) %>%
      dplyr::mutate(rank=rank,missing=missing,value = value / value[temperature == min(temperature)]) %>% unique(.)%>% 
      dplyr::rename("sample"="sample_id")
    df<-df[!is.na(df$Accession),]
    
  }else{
    df <- df %>%
      dplyr::select(Accession, sample_id, temperature,value,rank,missing,missing_pct,sample_name) %>% 
      dplyr::filter(!is.na(temperature),!is.na(value)) %>%
      dplyr::group_by(Accession,sample_id) %>%
      dplyr::mutate(rank=rank,missing=missing,value = value / value[temperature == min(temperature)]) %>% unique(.) %>% 
      dplyr::rename("sample"="sample_id")
    df<-df[!is.na(df$Accession),]
    
  }
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
normalize_cetsa <- function(df, temperatures,PSM=FALSE,Batch=TRUE,filters=FALSE) {
  temperatures <- sort(temperatures)
  df$Accession<-as.factor(df$Accession)
  
  df.jointP <- suppressWarnings(df %>%
                                  dplyr::group_split(Accession,sample) %>% 
                                  purrr::map(function(x) x %>% dplyr::mutate(n=dplyr::n()) %>% 
                                               dplyr::mutate(.,T7 = try(mean(value[x$temperature == temperatures[7]]/value[temperature == temperatures[1]])),
                                                             T9 = try(mean(value[temperature == temperatures[9]]/value[temperature == temperatures[1]])),
                                                             T10 = try(mean(value[temperature == temperatures[10]]/value[temperature == temperatures[1]])))))
  df.jointP<- dplyr::bind_rows(df.jointP)
  if(isTRUE(filters)){
    df.jointP<-df.jointP %>% dplyr::filter(T7 >= 0.4 & T7 <= 0.6)
    df.jointP<-df.jointP %>% dplyr::filter(T9 < 0.3)%>% dplyr::select(-T7,-T9,-n)
    if(any(names(df.jointP)=="T10")){
      df.jointP<- df.jointP %>% dplyr::filter(T10 < 0.2)%>% dplyr::select(-T10)#normalization from TPP
    }
  }
  if(nrow(df)==0){
    return(warning("Please disable filters, all data was filtered out."))
  }

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
    dplyr::rowwise(.) %>% 
    dplyr::mutate(fitted_values = ifelse(!is.logical(fit),list(data.frame(fitted_values=predict(fit))),NA),
                  temperature=ifelse(!is.logical(fit),list(data.frame(temperature=fit$data$t)),NA)) %>% 
    dplyr::select(sample,fitted_values,temperature) 
  
  
  
  ## calculate the fitted values
  d<-length(df.fit$fitted_values[[1]])
  #unnest fitted values from list and name value column
  check<-df.fit %>% unnest(c(fitted_values,temperature))
  
  check<-check %>% unique(.) 
  check$sample<-as.factor(check$sample)
  
  test<-df.median %>% dplyr::group_by(sample,temperature) %>% dplyr::right_join(check,c('sample','temperature'))
  ## calculate ratios between the fitted curves and the median values
  df.out <- test %>%
    dplyr::mutate(correction = ifelse(is.na(fitted_values / value),NA,fitted_values / value)) %>%
    dplyr::select('sample','temperature','value','fitted_values','correction')
  df.out<-df.out %>% dplyr::select(-fitted_values,-value)
  ## apply normalization factor to data
  
  df$correction<-df.out$correction
  df <- df %>% 
    dplyr::mutate(norm_value = ifelse(is.na(correction),value,value * correction)) %>% dplyr::ungroup(.)
  df <- df %>% 
    dplyr::mutate(norm_value = ifelse(is.na(correction),value,value * correction)) %>% dplyr::ungroup(.)
  df <- df %>% 
    dplyr::rename("uniqueID"="Accession", "C"="temperature","I"="norm_value")
  df<-df %>% dplyr::ungroup(.) %>% dplyr::select(-value,-correction)
  if(isTRUE(filters)){
    df<-df %>% dplyr::select(-T10,-T7,-T9)
  }

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
  pars<-f$fit[[1]]$m$getPars()
  Tm<-pars['a']/(pars['b'] - log(1-pars['Pl'])/(1/2 - pars['Pl']-1))
  return(round(Tm,1)[['a']]) 
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
###Upset plots#################
#To generate information on missing values
################################
#df_ is the input from df_norm which includes missing percentages
#returns segmented data for vehicle and treated
upMV <- function(df_,condition,N,plot_multiple=FALSE){
  if(isTRUE(plot_multiple)){
    df_1<-dplyr::bind_rows(df_)
    df_1<-df_1%>% dplyr::rename("uniqueID"="Accession","I"="value","C"="temperature")
    df_1$rank<-ifelse(df_1$rank==1,"high",ifelse(df_1$rank==2,"medium","low"))
    df_1$rank<-factor(df_1$rank,levels=c("high","medium","low"))
    df_1$uniqueID<-as.factor(df_1$uniqueID)
    df_1$dataset<-as.factor(df_1$dataset)
    df_1$sample_name<-as.factor(df_1$sample_name)
    df_1<-df_1%>% 
      dplyr::group_split(uniqueID,sample_id)
    df_1<-lapply(df_1,function(x) x[1,]) 
    df_1<-dplyr::bind_rows(df_1)
   
    df_1$missing_pct<-as.numeric(df_1$missing_pct)
    df_1$missing_pct<-round(df_1$missing_pct,0)

    
    df_1<-df_1 %>% dplyr::group_split(sample_name)
    df_1<-purrr::map(df_1,function(x)
      pivot_wider(x,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA) %>%
        dplyr::select(-C,-I,-CC,-missing,-dataset,-sample_id,-temp_ref,-id))
    
    df_2<-purrr::map(df_1,function(x) x %>% dplyr::select(rank,sample_name))
    df_1<-purrr::map(df_1,function(x) 1+x[,lapply(x,class)=="numeric"])
    df<-lapply(df_1,function(x) colnames(x))
    df<-lapply(df,function(x) str_replace(x,x,
                                          paste0("missing ",x,"%", sep = " ")))
    
    
    
    
    df_1<-purrr::map2(df_1,df,function(x,y){setNames(x,y)})
    df_1<-purrr::map2(df_1,df_2,function(x,y)cbind(x,y))
    df_1<-purrr::map(df_1,function(x)x[!is.na(x$rank),])
    
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    check<-list()
    check<- purrr::map(df_1,function(x)upset(x,names(x),
                                             min_degree=1,
                                             set_sizes=FALSE,
                                             guides='collect',
                                             n_intersections=N,
                                             height_ratio = 0.7,
                                             stripes='white',
                                             base_annotations=list(
                                               '# of Replicates'=intersection_size(
                                                 counts=TRUE,
                                               )
                                             ),
                                             annotations =list(
                                               "Bin %"=list(
                                                 aes=aes(x=intersection, fill=x$rank),
                                                 
                                                 geom=list(
                                                   geom_bar(stat='count', position='fill', na.rm=TRUE,show.legend=FALSE),
                                                   
                                                   geom_text(
                                                     aes(
                                                       label=!!aes_percentage(relative_to='intersection'),
                                                       color=ifelse(!is.na(rank), 'show', 'hide')
                                                     ),
                                                     stat='count',
                                                     position=position_fill(vjust = .5),
                                                     
                                                   ),
                                                   
                                                   scale_y_continuous(labels=scales::percent_format()),
                                                   show_hide_scale,
                                                   rating_scale
                                                   
                                                 )
                                               )
                                             )
                                             
    )+ggtitle(paste0(x$sample_name[1])))
    check1<-upset(df_1[[1]],names(df_1[[1]]),min_degree=1,
                                             set_sizes=FALSE,
                                             guides='collect',
                                             n_intersections=N,
                                             
                                             stripes='white',
                                             base_annotations=list(
                                               '# of Replicates'=intersection_size(
                                                 counts=TRUE,
                                               )
                                             ),
                                             annotations =list(
                                               "Ranked Intensity  "=list(
                                                 aes=aes(x=intersection, fill=df_1[[1]]$rank),
                                                 
                                                 geom=list(
                                                   geom_bar(stat='count', position='fill', na.rm=TRUE),
                                                   themes=theme(legend.position="bottom", legend.box = "horizontal"),
                                                   geom_text(
                                                     aes(
                                                       label=!!aes_percentage(relative_to='intersection'),
                                                       color=ifelse(!is.na(rank), 'show', 'hide')
                                                     ),
                                                     stat='count',
                                                     position=position_fill(vjust = .5),
                                                     
                                                   ),
                                                   
                                                   scale_y_continuous(labels=scales::percent_format()),
                                                   show_hide_scale,
                                                   rating_scale
                                                   
                                                 )
                                               )
                                             )
                                             
    )+ggtitle(paste0(df_1[[1]]$sample_name[1]))
   y<-get_legend(check1$patches$plots[[1]])
   
    P<-ggarrange(plotlist=check,ncol=4,nrow=2,font.label = list(size = 14, color = "black", face = "bold"),labels = "AUTO",legend.grob = y)
    print(P)
  }else{
    ########Plot separately
    
    df_<-df_ %>% subset(sample_name== condition) 
    
    df_<-df_ %>% dplyr::rename("uniqueID"="Accession","C"="temperature","I"="value")
    
    df_$rank<-ifelse(df_$rank==1,"high",ifelse(df_$rank==2,"medium","low"))
    df_$rank<-factor(df_$rank,levels=c("high","medium","low"))
    df_<-df_ %>% dplyr::filter(sample_name==condition)%>%
      dplyr::group_split(uniqueID) 
    df_<-lapply(df_,function(x) x[1,]) %>% dplyr::bind_rows(.)
    df_$missing_pct<-as.numeric(df_$missing_pct)
    df_$missing_pct<-round(df_$missing_pct,0)
    df_$dataset<-as.factor(df_$dataset)
    
    d2<-dplyr::bind_rows(df_)%>% dplyr::filter(dataset=="vehicle") 
    d2$uniqueID<-as.factor(d2$uniqueID)
    d2$sample_name<-as.factor(d2$sample_name)
    d2$C<-as.factor(d2$C)
    d2$rank<-as.factor(d2$rank)
    
    
    d1<-dplyr::bind_rows(df_)%>% dplyr::filter(dataset=="treated") 
    d1$uniqueID<-as.factor(d1$uniqueID)
    d1$sample_name<-as.factor(d1$sample_name)
    d1$C<-as.factor(d1$C)
    d1$rank<-as.factor(d1$rank)
    
    
    d3<-rbind(d1,d2)
    
    d3<-tidyr::pivot_wider(d3,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d1<-pivot_wider(d1,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    d2<-pivot_wider(d2,names_from=c(missing_pct),values_from=missing_pct,values_fill=NA)
    
    d1<-d1%>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    d2<-d2 %>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    d3<-d3 %>% dplyr::select(-uniqueID,-C,-I,-CC,-missing,-dataset,-sample_id,-sample_name)
    
    rd1<-d1$rank
    rd2<-d2$rank
    rd3<-d3$rank
    
    d1<-1+d1[,lapply(d1,class)=="numeric"]
    d2<-1+d2[,lapply(d2,class)=="numeric"]
    d3<-1+d3[,lapply(d3,class)=="numeric"]
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    colnames(d1) <- paste("missing ", colnames(d1),"%", sep = " ")
    colnames(d2) <- paste("missing ", colnames(d2),"%", sep = " ")
    colnames(d3) <- paste("missing ", colnames(d3),"%", sep = " ")
    
    rating_scale = scale_fill_manual(name="Ranked Intensity",
                                     values=c(
                                       'high'='#2ca25f', 'medium'='#99d8c9',
                                       'low'='#e5f5f9'
                                     ))
    show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide=FALSE)
    
    d1$rank<-rd1
    d2$rank<-rd2
    d3$rank<-rd3
    
    p<-upset(d2,names(d2),
             min_degree=1,
             set_sizes=FALSE,
             n_intersections=N,
             stripes='white',
             base_annotations=list(
               '# of Replicates'=intersection_size(
                 counts=TRUE,
               )
             ),
             annotations =list(
               "Ranked Intensity %"=list(
                 aes=aes(x=intersection, fill=d2$rank),
                 geom=list(
                   geom_bar(stat='count', position='fill', na.rm=TRUE),
                   geom_text(
                     aes(
                       label=!!aes_percentage(relative_to='intersection'),
                       color=ifelse(!is.na(rank), 'show', 'hide')
                     ),
                     stat='count',
                     position=position_fill(vjust = .5)
                   ),
                   scale_y_continuous(labels=scales::percent_format()),
                   
                   show_hide_scale,
                   rating_scale
                 )
               )
             )
             
    )+ggtitle(paste0( condition," vehicle"))
    
    p2<-upset(d1,names(d1),name="Sample bins",
              min_degree=1,
              n_intersections=N,
              set_sizes=FALSE,
              stripes='white',
              base_annotations=list(
                '# of Replicates'=intersection_size(
                  counts=TRUE,
                )
              ),
              annotations =list(
                "Ranked Intensity %"=list(
                  aes=aes(x=intersection, fill=d1$rank),
                  geom=list(
                    geom_bar(stat='count', position='fill', na.rm=TRUE),
                    
                    geom_text(
                      aes(
                        label=!!aes_percentage(relative_to='intersection'),
                        color=ifelse(!is.na(rank), 'show', 'hide')
                      ),
                      stat='count',
                      position=position_fill(vjust = .5)
                    ),
                    scale_y_continuous(labels=scales::percent_format()),
                    
                    show_hide_scale,
                    rating_scale
                  )
                )
              )
              
    )+ggtitle(paste0( condition," treated"))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    
    d<-names(d3)[which(!names(d3)=="rank")]
    dr<-which(names(d3)=="rank")
    d3<-d3 %>% dplyr::select(-rank)
    d<-as.numeric(unlist(str_extract_all(d,"[[:digit:]]+")))
    ints<-names(d3)[c(order(d,decreasing=FALSE))] %>% as.list(.)
    d3$rank<-rd3
    
    p3<-ComplexUpset::upset(d3,names(d3), 
                            n_intersections=N,
                            set_sizes=FALSE,
                            stripes='white',
                            base_annotations=list(
                              '# of Replicates'=intersection_size(
                                counts=TRUE,
                              )
                            ),
                            annotations =list(
                              "Ranked Intensity %"=list(
                                aes=aes(x=intersection, fill=d3$rank),
                                geom=list(
                                  geom_bar(stat='count', position='fill', na.rm=TRUE),
                                  
                                  geom_text(
                                    aes(
                                      label=!!aes_percentage(relative_to='intersection'),
                                      color=ifelse(!is.na(rank), 'show', 'hide')
                                    ),
                                    stat='count',
                                    position=position_fill(vjust = .5)
                                  ),
                                  scale_y_continuous(labels=scales::percent_format()),
                                  
                                  show_hide_scale,
                                  rating_scale
                                )
                              )
                            ) 
                            
    )+ggtitle(paste0(condition))+ guides(color=guide_legend(title="Ranked Intensity"))
    
    return(list(plot(p),plot(p2),plot(p3)))
  }
}

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
  
  df_1<-purrr::map(d, function(x) { 
    x %>%
      dplyr::group_by(C) %>%
      dplyr::mutate(I=mean(I,na.rm=TRUE)) %>% 
      dplyr::ungroup(.)
    
  })
  #rank intensity values using 3 regions,  rename column as LineRegion
  LR<-purrr::map(df_1, function(x) { 
    dplyr::ntile(dplyr::desc(x$I),3)%>%
      as.data.frame(.) %>% dplyr::rename("LineRegion"=".")})
  df_1<-purrr::map(df_1,function(x){x %>% dplyr::select(-LineRegion)})#remove Line Region column from one dataset before merging
  
  #Add LR to the list
  df_1<-purrr::map2(df_1,LR, function(x,y) {c(x,y) %>% as.data.frame(.)})
  df_1 <-purrr::map(df_1,function(x){x %>% dplyr::mutate(C = C,I=I,CC=as.factor(CC))})
  
  #separate by Line Regions
  df1<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==1) %>% as.data.frame(.)})
  df2<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==2) %>% as.data.frame(.)})
  df3<-purrr::map(df_1,function(x){x %>% dplyr::filter(LineRegion==3) %>% as.data.frame(.)})
  
  #preallocate model data per line region
  LM1<-list(NA)
  LM2<-list(NA)
  LM3<-list(NA)
  df1<-purrr::map(df1,function(x) x[order(x$C),])
  df2<-purrr::map(df2,function(x) x[order(x$C),])
  df3<-purrr::map(df3,function(x) x[order(x$C),])
  # #Flag NA values
  df1<-purrr::map(df1,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df2<-purrr::map(df2,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  df3<-purrr::map(df3,function(x) x %>% dplyr::mutate(missing=is.na(x$I)))
  # #remove NA values
  df1<-purrr::map(df1,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df2<-purrr::map(df2,function(x) x %>% dplyr::filter(!is.na(x$I)))
  df3<-purrr::map(df3,function(x) x %>% dplyr::filter(!is.na(x$I)))
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
  L1<-purrr::map(df1,function(x){ 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM1<-purrr::map2(df1,L1,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L2<-purrr::map(df2,function(x){ 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM2<-purrr::map2(df2,L2,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  L3<-purrr::map2(df3,seq(df3),function(x,y) { 
    tryCatch(lm(formula = I~C,data = x ,na.action='na.omit'), error = function(e){NA})})
  LM3<-purrr::map2(df3,L3,function(x,y) x %>% purrr::keep(function(x) any(!is.na(y))))
  
  #linear fit per line region
  LM1<-purrr::map2(df1,L1,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM2<-purrr::map2(df2,L2,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  LM3<-purrr::map2(df3,L3,function(x,y)x %>% dplyr::mutate(M1 = list(y)))
  
  
  #fitted curves
  x1<-purrr::map(LM1, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x2<-purrr::map(LM2, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  x3<-purrr::map(LM3, function(x) try(ifelse(class(x$M1[[1]])=="lm",TRUE,NA)))
  
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
  df_0<-purrr::map(df_0, function(x) x %>% dplyr::select(-missing,-CI))
  
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
  df_0<-purrr::map(df_0,function(x)x %>% dplyr::arrange(C) %>% dplyr::group_by(C,LineRegion)%>%dplyr::mutate(n=dplyr::n()) %>% dplyr::ungroup())
  
  #subset of the data with shared line region values, using purrr::map to keep size constant
  Split<-purrr::map(df_0,function(x)x %>%subset(n<max(n)) %>% data.frame(.) %>% dplyr::mutate(LineRegion=as.numeric((.$LineRegion))))
  #remove NA values
  Split<-dplyr::bind_rows(Split)
  df_0<-dplyr::bind_rows(df_0)
  #Join df_0 with the subset of values
  dap<-df_0 %>% dplyr::left_join(Split,by=c("C"="C","uniqueID"="uniqueID","dataset"="dataset","I"="I","sample_name"="sample_name","CC"="CC","n"="n","sample"="sample","missing_pct"="missing_pct","rank"="rank"))
  dap<-dap %>% dplyr::group_split(uniqueID)
  dap<-purrr::map(dap,function(x)x %>% dplyr::mutate(LineRegion=ifelse(.$C<=tail(.$C[.$LineRegion.x==1],1),1,ifelse(.$C<=tail(.$C[.$LineRegion.x==2],1),2,3))) %>% dplyr::select(-LineRegion.x,-LineRegion.y,-missing.y))
  
  return(dap)
}
#Trilinear functions
tlstat<-function(DF,df,df1,norm=FALSE,Filters=FALSE,Ftest=FALSE){
  i<-1
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
  
  if(!isTRUE(norm)){
    mean1<-list()
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                         dplyr::group_nest(LineRegion,uniqueID) %>%
                         dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=FALSE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=purrr::map(M1,function(x){deviance(x)}),
                                       Rsq=purrr::map(M1,function(x){summary(x)$r.squared}), 
                                       dataset="vehicle",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=min(x$I,na.rm=FALSE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}), 
                                         dataset="treated",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
    
    
    
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                         dplyr::group_nest(LineRegion,uniqueID) %>% 
                         dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx( x$I,x$C, xout=min(x$I,na.rm=FALSE)+(0.5*(max(x$I, na.rm=TRUE)-min(x$I, na.rm=TRUE))))$y),
                                       slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="null",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    if (isTRUE(Filters)){
      #Apply lax Rsq and negative slope filter to remove flat melt curves
      mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      
      mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    }
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
    if(isTRUE(Ftest)){
      #Calculate rss0 and rss1 null vs alt
      rss0<-purrr::map(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
      rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
                                                              Tm = y$Tm[[1]]-x$Tm[[1]]))
      #params for null and alternative models
      pN<-purrr::map(mean3,function(x)x %>% dplyr::summarise(pN = 4))
      pA<-purrr::map(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
      
      #sum residuals
      n1<-purrr::map2(mean1,mean1_1,function(x,y) data.frame(n1 = as.numeric(nrow(dplyr::bind_rows(x$data))) + as.numeric(nrow(dplyr::bind_rows(y$data)))))
      #degrees of freedom before
      d1<-purrr::map2(pA,pN,function(x,y)data.frame(d1=x$pA-y$pN))
      d2<-purrr::map2(n1,pA,function(x,y)data.frame(d2=x$n1-y$pA)) 
      #delta RSS
      rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS-y$RSS %>% as.data.frame(.))
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
    }
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    return(results)
  }else if (isTRUE(norm)){
    mean1<-list()
    mean1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df[[i]]$uniqueID[1],Tm=rep(0,1))
    mean1<- purrr::map(df,function(x) x %>% as.data.frame(.) %>% 
                         dplyr::group_nest(LineRegion,uniqueID) %>%
                         dplyr::mutate(M1=purrr::map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=0.5))$y,
                                       slope=purrr::map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=purrr::map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="vehicle",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    #define linear models with outputs
    
    mean1_1<-list()
    mean1_1[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="treated",uniqueID=df1[[i]]$uniqueID[1],Tm=rep(0,1))
    
    mean1_1<- purrr::map(df1,function(x) x %>% as.data.frame(.) %>%
                           dplyr::group_nest(LineRegion,uniqueID) %>% 
                           dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                         CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                         Tm=with(x, stats::approx(x$I,x$C, xout=0.5))$y,
                                         slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                         intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                         rss=map(M1,function(x){deviance(x)}),
                                         Rsq=map(M1,function(x){summary(x)$r.squared}),
                                         dataset="treated",
                                         uniqueID=x$uniqueID[1],
                                         n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
    
    
    
    # null hypothesis
    #null
    mean3<-list()
    mean3[[1]]<-data.frame(slope=rep(0,1),intercept=rep(0,1),rss=rep(0,1),Rsq=rep(0,1),AUC = rep(0,1),dataset="null",uniqueID=DF[[i]]$uniqueID[1],Tm=rep(0,1))
    
    
    mean3<- purrr::map(DF,function(x) x %>% as.data.frame(.) %>%
                         dplyr::group_nest(LineRegion,uniqueID) %>% 
                         dplyr::mutate(M1=map(data,function(x){stats::lm(x$I ~ x$C)}),
                                       CI=purrr::map(M1,function(x){predict(x,interval="confidence")}),
                                       Tm=with(x, stats::approx(x$I,x$C, xout=0.5))$y,
                                       slope=map(M1,function(x){as.numeric(coef(x)[2])}),
                                       intercept=map(M1,function(x){as.numeric(coef(x)[1])}),
                                       rss=map(M1,function(x){deviance(x)}),
                                       Rsq=map(M1,function(x){summary(x)$r.squared}),
                                       dataset="null",
                                       uniqueID=x$uniqueID[1],
                                       n=ifelse(class(M1)=="lm",1,0)))
    
    
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(AUC = pracma::trapz(x$M1[[1]]$fitted.values)))
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
    
    
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    
    if(isTRUE(Filters)){
      
      #Apply lax Rsq and negative slope filter to remove flat melt curves
      mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) all(unlist(x$Rsq)>0.5)))
      
      mean1<-suppressWarnings(mean1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      mean1_1<-suppressWarnings(mean1_1 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
      mean3<-suppressWarnings(mean3 %>% purrr::keep(function(x) any(unlist(x$slope)<0)))
    }
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
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    if(isTRUE(Ftest)){
      #Calculate rss0 and rss1 null vs alt
      rss0<-purrr::map(mean3,function(x)data.frame(RSS = sum(as.numeric(x$rss))))
      rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = sum(as.numeric(x$rss))+sum(as.numeric(y$rss)),
                                                              Tm = y$Tm[[1]]-x$Tm[[1]]))
      #params for null and alternative models
      pN<-purrr::map(mean3,function(x)x %>% dplyr::summarise(pN = 4))
      pA<-purrr::map(mean1_1,function(x)x %>% dplyr::summarise(pA = 8))
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
        coord_cartesian(xlim=c(0,100))+
        ggplot2::xlab("F-values")
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
      test<-testScaled %>% dplyr::filter(.$pAdj<0.01)
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
    }
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
    #apply Tm data from LR 2
    mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
    
    results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
    return(results)
  }
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
  #apply Tm data from LR 2
  mean1<-purrr::map(mean1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
  mean1_1<-purrr::map(mean1_1,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
  mean3<-purrr::map(mean3,function(x) x %>% dplyr::mutate(Tm=x[which(x$LineRegion==2),"Tm"]))
  
  results<-dplyr::bind_rows(mean1,mean1_1,mean3) %>% dplyr::group_split(uniqueID)
  
  return(results)
}
tlf<-function(tlresults,DFN,APfilt=TRUE,PF=TRUE){
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
  Nsum<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "null")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
  
  #get the summed rss values for vehicle
  Rssv<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "vehicle")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>%head(.,1))
  #get the summed rss values for treated
  Rsst<-purrr::map(tlresults, function(x) x %>% subset(stringr::str_detect(tolower(dataset), pattern = "treated")) %>% 
                     dplyr::rowwise(.) %>%  dplyr::mutate(RSS=sum(unlist(.$rss)))%>% dplyr::select(RSS,Tm,dataset,uniqueID)%>% head(.,1))
  #find the rss difference between treated and vehicle 
  
  Rssv<-purrr::map(Rssv,function(x)na.omit(x))
  Rsst<-purrr::map(Rsst,function(x)na.omit(x))
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
    Df1<-Dsum %>% dplyr::left_join(tlresults,by=c("uniqueID")) %>% as.data.frame(.) %>% dplyr::rename("dataset"="dataset.y") %>% 
      dplyr::select(-dataset.x)
    Df1<-Df1 %>% dplyr::group_split(uniqueID)
    
  }
  
  df1<-list()
  #get uniqueID and dataset for stable proteins with decreasing RSS differences
  df1<-purrr::map(Df1,function(x) x %>% dplyr::select(uniqueID,dataset) %>% head(.,1))
  df1<-data.frame(dplyr::bind_rows(df1))
  
  #unlist to data.frame
  #order the original data by RSS differences
  #
  DFN<- dplyr::bind_rows(DFN)
  DFN$uniqueID<-as.vector(DFN$uniqueID)
  df1$uniqueID<-as.vector(df1$uniqueID)
  
  
  df2<-df1 %>% dplyr::right_join(DFN,by=c("uniqueID")) %>% dplyr::rename("dataset"="dataset.y") %>% 
    dplyr::select(-dataset.x)
  
  
  return(list(df1,df2,Df1))
}
tlCI<-function(i,df1,df2,Df1,overlay=TRUE,residuals=FALSE){
  null<-data.frame()
  i<-i
  df1<-df1
  Df1<-Df1[[i]]
  DF1<-data.frame(NA)
  DF1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) 
  
  null<-Df1 %>% subset(dataset == "null")
  
  pred1<-predict(null$M1[[1]], interval="confidence") %>% as.data.frame(.)
  if(nrow(null)==2|nrow(null)==3){
    pred2<-predict(null$M1[[2]], interval="confidence")%>% as.data.frame(.)
    pred2<-na.omit(pred2)
  }else{
    pred2<-data.frame()
  }
  if(nrow(null)==3){
    pred3<-predict(null$M1[[3]], interval="confidence")%>% as.data.frame(.)
    pred3<-na.omit(pred3)
  }else{
    pred3<-data.frame()
  }
  Pred1<-NA
  pred1<-na.omit(pred1)
  
  
  
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
    ggplot2::geom_point(ggplot2::aes(x=C,y=I))+ ggplot2::ggtitle(paste(Df1$uniqueID[1],DF1$sample_name[1]))+
    ggplot2::geom_ribbon(data=Pred,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+ 
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ 
    annotate("text", x=60, y=min(Pred$I),label=paste("RSS= ",round(sum(unlist(null$rss)),3)))+
    annotate("text",
             x = Pred[which(round(Pred$fit,1)==0.5)[1],]$C,
             y = 0.45,
             label=paste(Pred[which(round(Pred$fit,1)==0.5)[1],]$C),
             colour="red"
    )
  
  
  DF_f<-df2 %>%subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated')) %>% subset(dataset=="vehicle")
  
  vehicle<-Df1 %>% subset(dataset == "vehicle")
  
  pred1<-predict(vehicle$M1[[1]], interval="confidence")%>% as.data.frame(.)
  if(nrow(vehicle)==2|nrow(vehicle)==3){
    pred2<-predict(vehicle$M1[[2]], interval="confidence")%>% as.data.frame(.)
    pred2<-na.omit(pred2)
  }else{
    pred2<-data.frame()
  }
  if(nrow(vehicle)==3){
    pred3<-predict(vehicle$M1[[3]], interval="confidence")%>% as.data.frame(.)
    pred3<-na.omit(pred3)
  }else{
    pred3<-data.frame()
  }
  Pred1<-NA
  pred1<-na.omit(pred1)
  
  
  
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
  
  
  DF_f1<-data.frame()
  DF_f1<-df2 %>% subset(uniqueID == df1$uniqueID[i]) %>% dplyr::mutate(dataset=ifelse(CC==0,'vehicle','treated'))
  DF_f1<-DF_f1 %>% subset(dataset =="treated")
  
  
  treated<-data.frame()
  treated<-Df1 %>% subset(dataset == "treated")
  
  pred1<-predict(treated$M1[[1]], interval="confidence")%>% as.data.frame(.)
  if(nrow(treated)==2|nrow(treated)==3){
    pred2<-predict(treated$M1[[2]], interval="confidence")%>% as.data.frame(.)
    pred2<-na.omit(pred2)
  }else{
    pred2<-data.frame()
  }
  if(nrow(treated)==3){
    pred3<-predict(treated$M1[[3]], interval="confidence")%>% as.data.frame(.)
    pred3<-na.omit(pred3)
  }else{
    pred3<-data.frame()
  }
  
  pred1<-na.omit(pred1)
  
  
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
  rn<-data.frame(residuals(null$M1[[1]]))
  if(nrow(null)==3){
    rn<-data.frame(residuals=c(residuals(null$M1[[1]]),residuals(null$M1[[2]]),residuals(null$M1[[3]])))
  }else{
    rn<-data.frame(residuals=c(residuals(null$M1[[1]]),residuals(null$M1[[2]])))
  }
  Pred<-cbind(Pred,rn[1:nrow(Pred),])
  names(Pred)<- names(Pred)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(residuals(vehicle$M1[[1]]))
  if(nrow(vehicle)==3){
    rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]]),residuals(vehicle$M1[[3]])))
  }else{
    rn<-data.frame(c(residuals(vehicle$M1[[1]]),residuals(vehicle$M1[[2]])))
  }
  Pred1<-cbind(Pred1,rn[1:nrow(Pred1),])
  names(Pred1)<- c("fit","lower","upper","C","I","Treatment",'residuals')
  rn<-data.frame(residuals(treated$M1[[1]]))
  if(nrow(treated)==3){
    rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]]),residuals(treated$M1[[3]])))
  }else{
    rn<-data.frame(c(residuals(treated$M1[[1]]),residuals(treated$M1[[2]])))
  }
  Pred2<-cbind(Pred2,rn[1:nrow(Pred2),])
  names(Pred2)<- names(Pred2)<-c("fit","lower","upper","C","I","Treatment",'residuals')
  Preds<-rbind(Pred1,Pred2)
  Preds$C<-as.numeric(as.vector(Preds$C))
  Preds$I<-as.numeric(as.vector(Preds$I))
  Pred2$C<-as.numeric(as.vector(Pred2$C))
  Pred2$I<-as.numeric(as.vector(Pred2$I))
  DF1$dataset<-as.factor(DF1$dataset)
  
  PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = residuals,color=Treatment)) +
    ggplot2::geom_point()+ ggplot2::ggtitle(paste(Df1$uniqueID[1],DF1$sample_name[1]))+
    ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
  if(isTRUE(residuals)){
  print(PLrs)
  }
  Tm_d<-round(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),1)
  PLR_P1<-ggplot2::ggplot(Pred1, ggplot2::aes(x = C,y = fit,color=Treatment))+ggplot2::geom_point(Pred1, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_ribbon(data=Pred1,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
    annotate("text",
             x = 2+round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
             y = 0.55,
             label=paste0(round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1)),
             colour="blue"
    )+
    annotate("segment", x = min(Pred1$C), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
             colour = "blue",linetype=2)+
    annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
             colour = "blue",linetype=2)
  
  
  PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
    ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+
    annotate("text",
             x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
             y = 0.55,
             label=paste0(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
             colour="red"
    )+
    annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
             colour = "red",linetype=2)+
    annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
             colour = "red",linetype=2)
  if(overlay=="TRUE"){
    AUCd<-round(P2_AUC-P1_AUC,2)
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
    miss_v<-data.frame(NA)
    miss_t<-data.frame(NA)
    miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
    miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
    Pred2$missing_v<-rep(round(miss_v$missing_pct[1],0),nrow(Pred2))
    Pred2$missing_t<-rep(round(miss_t$missing_pct[1],0),nrow(Pred2))
    
    Pred1$missing_v<-rep(miss_v$missing_pct[1],nrow(Pred1))
    Pred1$missing_t<-rep(miss_t$missing_pct[1],nrow(Pred1))
    
    
    PLR_P2<-PLR_P1+ggplot2::geom_point(Pred2, mapping=ggplot2::aes(x = C,y = I,color=Treatment)) +
      ggplot2::geom_ribbon(data=Pred2,ggplot2::aes(x=C,ymin=lower,ymax=upper,fill=Treatment),alpha=0.2)+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Df1$uniqueID[1],DF1$sample_name[1]))+
      ggplot2::annotate("text", x=min(Pred2$C)+5, y= 0.55, label= paste("\u03A3","RSS= ",round(sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "vehicle"),'rss']))+
                                                                                                 sum(unlist(Df1[stringr::str_detect(tolower(Df1$dataset), pattern = "treated"),'rss'])),3)))+
      ggplot2::annotate("text", x=min(Pred2$C)+5, y= 0.45, label=  paste("\u0394", "AUC = ",AUCd))+ 
      ggplot2::annotate("text", x=min(Pred2$C)+5, y= 0.35, label= paste("\u0394","Tm = ",round(Tm_d,1),"\u00B0C"))+ 
      ggplot2::annotate("text", x=min(Pred2$C)+5, y= 0.25, label= paste("missing vehicle",Pred2$missing_v[1],"%"))+ 
      ggplot2::annotate("text", x=min(Pred2$C)+5, y= 0.15, label= paste("missing treated",Pred2$missing_t[1],"%"))+
      annotate("text",
               x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1),
               y = 0.55,
               label=paste0(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
               colour="red"
      )+
      annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    par(mfrow=c(2,2))
    print(PLR_P2)
  }else if(overlay=="FALSE"){
    miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
    miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
    Pred2$missing_v<-rep(round(miss_v$missing_pct[1],0),nrow(Pred2))
    Pred2$missing_t<-rep(round(miss_t$missing_pct[1],0),nrow(Pred2))
    PLR<-PLR_P2+
      facet_wrap("Treatment") + 
      ggplot2::annotate("text", x=min(Pred$C), y=min(Pred2$I)+0.45, label= paste("missing % v",Pred2$missing_v[1]))+ 
      ggplot2::annotate("text", x=min(Pred$C), y=min(Pred2$I)+0.35, label= paste("missing % t",Pred2$missing_t[1]))+
      annotate("text",
               x = 2+round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1),
               y = 0.55,
               label=paste(round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1)),
               colour="red"
      )+
      annotate("segment", x = round(with(Pred1, stats::approx(Pred1$fit,Pred1$C,xout=max(Pred1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(Pred2, stats::approx(Pred2$fit,Pred2$C,xout=max(Pred2$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    print(PLR)
    
    if(bootstrap==TRUE){
      set.seed(233)
      n<-length(Pred$I)
      mean_orig_v<-mean(Pred$I)
      mean_orig_t<-mean(Pred2$I)
      N<-1000
      for (i in 1:N){
        bs_vehicle[i]<-sample(Pred$I,n,replace=TRUE)
        mean_v[i]<-mean(bs_vehicle[i])
        
        bs_treated[i]<-sample(Pred2$I,n,replace=TRUE)
        mean_t[i]<-mean(bs_treated[i])
      }
     mean_bs_v<-mean(mean_v)
     mean_bs_t<-mean(mean_t)
     bias_v<-mean_orig_v-mean_bs_v
     bias_t<-mean_orig_t-mean_bs_t
     
     summary<-data.frame(uniqueID=rep(Pred$uniqueID[1],2),
                         dataset=c(Pred$dataset[1],Pred2$dataset[1]),
                         Orig_mean=c(mean_orig_v,mean_orig_t),
                         BS_mean=c(mean_bs_v,mean_bs_t),
                         bias_bs=c(bias_v,bias_t),
                         sd_bs_v=sd(mean_v),
                         sd_bs_t=sd(mean_t),
                         CI_2.5_v=quantile(mean_v,0.025),
                         CI_97.5_v=quantile(mean_v,0.975),
                         CI_2.5_t=quantile(mean_t,0.025),
                         CI_97.5_t=quantile(mean_t,0.975))
                         
      
    }
  }
}


#Spline functions
spstat<-function(DF,df,df1,norm=FALSE,Ftest=TRUE){
  
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
  DF<- purrr::map(DF,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df<-df %>% dplyr::group_split(C,uniqueID,dataset) 
  df<- purrr::map(df,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  df1<-df1 %>% dplyr::group_split(C,uniqueID,dataset) 
  df1<- purrr::map(df1,function(x) x %>% dplyr::mutate(CV_pct = 100*sd(.$I,na.rm=TRUE)/mean(.$I,na.rm=TRUE)))
  
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
  DF<-purrr::map(DF,function(x)na.omit(x))
  df<-purrr::map(df,function(x)na.omit(x))
  df1<-purrr::map(df1,function(x)na.omit(x))
  if(!isTRUE(norm)){
    
    #alternative spline fit method : Generalized Additive Models
    #fit penalized splines
    m <- purrr::map(df,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(x$I ~ s(x$C,k=5), data = x , method = "REML")))))
    m<-m %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m<-purrr::map(m,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                   sum = list(summary(.$M1[[1]])),
                                                   Tm=ifelse(any(class(dplyr::first(.$M1))=="gam"),try(with(x, stats::approx(x$I,x$C,xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                                   rss=deviance(.$M1[[1]]),
                                                   CV_pct = .$CV_pct,
                                                   AUC = pracma::trapz(.$M1[[1]]$fit),
                                                   rsq=summary(x$M1[[1]])$r.sq,
                                                   n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0)))
    #m<-lapply(m,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    m1 <- purrr::map(df1,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data = x , method = "REML")))))
    m1<-m1 %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    m1<-purrr::map(m1,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                     sum = list(summary(.$M1[[1]])),
                                                     Tm=ifelse(any(class(dplyr::first(.$M1))=="gam"),try(with(x, stats::approx(x$I,x$C,xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                                     rss=deviance(.$M1[[1]]),
                                                     CV_pct = .$CV_pct,
                                                     AUC = pracma::trapz(.$M1[[1]]$fit),
                                                     rsq=summary(x$M1[[1]])$r.sq,
                                                     n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0)))
    #m1<-lapply(df1,function(x) x %>% dplyr::mutate(sig = ifelse(sum[[1]]$p.pv[[1]]<0.05,list(mgcv::gam(I ~ s(C,k=k_[[1]]-1), data = x , method = "REML")),"ns")))
    
    
    mn<- purrr::map(DF,function(x)x %>% dplyr::mutate(M1 = list(try(mgcv::gam(I ~ s(C,k=5), data =x, method = "REML")))))
    mn<-mn %>% purrr::keep(function(x)any(class(dplyr::first(x$M1))=="gam"))
    #check significance and refit data with more k 
    mn<-purrr::map(mn,function(x)x %>% dplyr::mutate(k_ = .$M1[[1]]$rank,
                                                     sum = list(summary(.$M1[[1]])),
                                                     Tm=ifelse(any(class(dplyr::first(.$M1))=="gam"),try(with(x, stats::approx(x$I,x$C,xout=max(x$I, na.rm=TRUE)-0.5))$y),NA),
                                                     rss=deviance(.$M1[[1]]),
                                                     CV_pct=.$CV_pct,
                                                     AUC = pracma::trapz(.$M1[[1]]$fit),
                                                     rsq=summary(.$M1[[1]])$r.sq,
                                                     n = ifelse(any(class(dplyr::first(.$M1))=="gam"),1,0)))
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
    return(results)
    if(isTRUE(Ftest)){
      #Calculate rss0 and rss1 null vs alt
      rss0<-purrr::map(mean3,function(x)data.frame(RSS = deviance(x$M1[[1]]),
                                                   RSS1 = deviance(x$sig[[1]])))
      rss1<-purrr::map2(mean1,mean1_1,function(x,y)data.frame(RSS = (deviance(x$M1[[1]])+deviance(y$M1[[1]])),
                                                              RSS1 =(deviance(x$sig[[1]])+deviance(y$sig[[1]])),
                                                              Tm = y$Tm[[1]]-x$Tm[[1]]))
      #params for null and alternative models
      pN<-purrr::map(mean3,function(x)x %>% dplyr::summarise(pN = x$M1[[1]]$rank,
                                                             n=x$M1[[1]]$df.null,
                                                             pN1 = x$sig[[1]]$rank,
                                                             n1 = x$sig[[1]]$df.null))
      pA<-purrr::map(mean1_1,function(x)x %>% dplyr::summarise(pA = 2*(x$M1[[1]]$rank),
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
      rssDiff<-purrr::map2(rss0,rss1,function(x,y) x$RSS1-y$RSS1 %>% as.data.frame(.))
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
      return(results)  
    }
    return(results)
  }
  return(results)
}
spf<-function(spresults,DFN,filters = TRUE){
  spresults<-spresults %>% dplyr::group_split(uniqueID)
  
  if(!isTRUE(filters))
  {
    sl<-purrr::map(seq_len(length(spresults)),function(x) as.numeric({paste(x)})) 
    sp<-purrr::map2(spresults,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))  
    sp<-dplyr::bind_rows(sp)  
    df1<-data.frame(uniqueID = unique(sp$uniqueID))  
    df2<-dplyr::bind_rows(DFN)  
    df2$C<-as.numeric(as.vector(df2$C)) 
    df2$I<-as.numeric(as.vector(df2$I))  
    df2<-sp %>% left_join(df2, by = c("uniqueID"="uniqueID","dataset"="dataset","C"="C","I"="I","CC"="CC","sample_name"="sample_name","LineRegion"="LineRegion","missing.x"="missing.x","missing_pct"="missing_pct","rank"="rank","N"="N","n"="n")) 
    Df1<-spresults
  }else{
    #Apply filters 
    #keep the positive AUC differences
    
    spresults<-spresults %>% keep(function(x) mean(x$AUC[x$dataset=="treated"],na.rm=TRUE)>mean(x$AUC[!x$dataset=="vehicle"],na.rm=TRUE))
    
    spresults<-spresults %>% keep(function(x) max(x$lambda)<1)
    
    if (is.null(nrow(spresults))){
      return(warning("all proteins filtered out by AUC and lambda value"))
    }
    #get Tm and RSS differences
    sp<-purrr::map(spresults, function(x) x %>% dplyr::mutate(Tmd= x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'Tm'][[1]] - x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'Tm'][[1]],
                                                              RSSd = sum(x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']) - sum(x[!stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "null"),'rss']),
                                                              AUCd = x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "treated"),'AUC'])[[1]]- x[stringr::str_detect(tolower(data.frame(x)$dataset), pattern = "vehicle"),'AUC'][[1]])
    #conserve list indexes
    sl<-purrr::map(seq_along(length(sp)),function(x) as.numeric({paste(x)}))
    
    #insert list index column
    sp<-purrr::map2(sp,sl,~.x %>% dplyr::mutate(id = as.numeric(.y)))
    sp<-dplyr::bind_rows(sp)
    
    sp<-dplyr::arrange(sp,dplyr::desc(AUCd),dplyr::desc(RSSd),dplyr::desc(Tmd)) %>% dplyr::select(uniqueID,id) %>% unique(.) 
    #arrange results by decreasing AUCd, RSSd and Tmd and standardize the order in spresults
    #Df1 holds the model results and stats for splines 
    
    df1<-data.frame(uniqueID = unique(sp$uniqueID))
    df2<-dplyr::bind_rows(DFN) 
    df2$C<-as.numeric(as.vector(df2$C))
    df2$I<-as.numeric(as.vector(df2$I))
    df2<-sp %>% 
      left_join(df2, by = c("uniqueID"="uniqueID","dataset"="dataset","C"="C","I"="I","CC"="CC","sample_name"="sample_name","LineRegion"="LineRegion","missing.x"="missing.x")) %>% 
      dplyr::rename("missing"="missing.x")
    Df1<-spresults[sp$id]
  }
  ret<-list()
  ret[[1]]<-df1
  ret[[2]]<-df2
  ret[[3]]<-Df1
  return(ret)
}
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}
spCI<-function(i,df1,df2,Df1,overlay=TRUE,alpha){
  null<-data.frame()
  i<-i
  
  #set C and I as numeric
  df2$C<-as.numeric(as.vector(df2$C))
  df2$I<-as.numeric(as.vector(df2$I))
  df2<-df2  %>%  mutate_if(is.logical,as.numeric) 
  df2$uniqueID<-as.character(df2$uniqueID)
  
  
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
  pred1$AUC<-pracma::trapz(pred1$fit)-pracma::trapz(pred$fit)
  pred1$AUC<-round(pred1$AUC[1],3)
  if( pred1$AUC[1] > 5){
    pred1$AUC<-pracma::trapz(pred1$lwrP)-pracma::trapz(pred$uprP)#AUC diff in stabilized CI
    pred1$AUC<-round(pred1$AUC,3)
    pred1$AUC<-pracma::trapz(pred1$lwrP)-pracma::trapz(pred$uprP)#AUC diff in stabilized CI
    pred1$AUC<-round(pred1$AUC,3)
  }else if ( pred1$AUC[1]< -5){
    pred1$AUC<-pracma::trapz(pred$lwrP)-pracma::trapz(pred1$uprP)#AUC diff in destabilized CI
    pred1$AUC<-round(pred1$AUC,3)
  }else{
    pred1$AUC<-round(pracma::trapz(pred1$fit)-pracma::trapz(pred$fit),2) #AUC diff in fit
  }
  pred1$RSS<- deviance(m1)-deviance(m)#RSS diff(+ stabilized)
  
  pred1$RSS<- round(pred1$RSS,3)
  #Residuals
  
  pred1$Tm<-round(treated$Tm[1]-vehicle$Tm[1],1)
  #missing values
  miss_v<-data.frame(NA)
  miss_t<-data.frame(NA)
  miss_v<-DF1%>% dplyr::filter(dataset=="vehicle")
  miss_t<-DF1 %>% dplyr::filter(dataset=="treated")
  
  Pred<-data.frame(m$fitted.values,m$residuals)
  names(Pred)<- c("fit","rn")
  Pred$dataset<-as.factor("vehicle")
  BSVar$dataset<-as.factor("vehicle")
  Pred1<-data.frame(m1$fitted.values,m1$residuals)
  names(Pred1)<-c("fit","rn")
  Pred1$dataset<-as.factor("treated")
  Preds<-rbind(Pred,Pred1)
  BSVar1$dataset<-as.factor("treated")
  #get fitted value data
  fitted.values<-data.frame(C=BSVar$M1[[1]]$model$`x$C`,fit=predict.gam(BSVar$M1[[1]],se.fit=TRUE))
  names(fitted.values)<-c("C","fit","se.fit")
  fitted.values1<-data.frame(C=BSVar1$M1[[1]]$model$C,fit=predict.gam(BSVar1$M1[[1]],se.fit=TRUE))
  names(fitted.values1)<-c("C","fit","se.fit")
  #append missing value data
  BSVar$missing_v<-rep(max(miss_v$missing_pct,na.rm=TRUE),nrow(BSVar))
  BSVar$missing_t<-rep(max(miss_t$missing_pct,na.rm=TRUE),nrow(BSVar))
  #append Tm values on predicted data
  pred1$Tm<-round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)-round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)
  
  PLrs<-ggplot2::ggplot(Preds, ggplot2::aes(x =fit,y = rn,color=dataset)) +ggplot2::geom_point()+ 
    ggplot2::ggtitle(paste(Df1[[i]]$uniqueID[1]," ",Df1[[i]]$dataset[1]))+ggplot2::xlab("Fitted Intensities")+ggplot2::ylab("Residuals")
  print(PLrs)
  plot1<-ggplot2::ggplot(BSVar,ggplot2::aes(x =C,y = I,color=dataset))+
    ggplot2::geom_point(BSVar,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
    ggplot2::geom_ribbon(data.frame(pred),mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
    ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(c(as.character(df1),"alternative"))+
    # ggplot2::annotate("text", x=50, y=1, label= paste("missing values: vehicle",vmissing[1]))+
    # ggplot2::annotate("text", x=50, y=0.9, label= paste("missing values treated",BSvar1$Dataset.x[1],":",tmissing[1]))                    
    ggplot2::annotate("text", x=min(BSVar$C)+5, y=min(BSVar$I,na.rm=TRUE)+0.55, label= paste("\u03A3","RSS= ", abs(pred1$RSS[1])))+
    ggplot2::annotate("text", x=min(BSVar$C)+5, y=min(BSVar$I,na.rm=TRUE)+0.45, label=  paste("\u0394", "AUC = ",pred1$AUC[1]))+
    ggplot2::annotate("text", x=min(BSVar$C)+5, y=min(BSVar$I,na.rm=TRUE)+0.35, label= paste("\u0394","Tm = ",round(pred1$Tm[1],1),"\u00B0C"))+ 
    ggplot2::annotate("text", x=min(BSVar$C)+5, y=min(BSVar$I)+0.25, label= paste("missing vehicle",round(BSVar$missing_v[1],0),"%"))+ 
    ggplot2::annotate("text", x=min(BSVar$C)+5, y=min(BSVar$I)+0.15, label= paste("missing treated",round(BSVar$missing_t[1],0),"%"))+
    annotate("text",
             x = 2+round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1),
             y = 0.55,
             label=paste0(round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1)),
             colour="blue"
    )+
    annotate("segment", x = min(fitted.values$C), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
             colour = "blue",linetype=2)+
    annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
             colour = "blue",linetype=2)
  
  
  plot<-plot1+
    ggplot2::geom_point(BSvar1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
    ggplot2::geom_ribbon(pred1,mapping=ggplot2::aes(x=C,y=fit,ymin = lwrP, ymax = uprP ,fill=CI), alpha = 0.2 ) +
    ggplot2::labs(y = "Relative Solubility",
                  x = "Temperature (\u00B0C)")+
    coord_cartesian(xlim = c(37,67))+annotate("text",
                                              x = 2+round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1),
                                              y = 0.55,
                                              label=paste0(round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1)),
                                              colour="red"
    )+
    annotate("segment", x = round(with(fitted.values, stats::approx(fitted.values$fit,fitted.values$C,xout=max(fitted.values$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0.5, yend = 0.5,
             colour = "red",linetype=2)+
    annotate("segment", x = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), xend = round(with(fitted.values1, stats::approx(fitted.values1$fit,fitted.values1$C,xout=max(fitted.values1$fit, na.rm=TRUE)-0.5))$y,1), y = 0, yend = 0.5,
             colour = "red",linetype=2)
  
  print(plot)  
}

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

sigC<-function(vehicle,treated,null){
  df_<-vehicle
  df_1<-treated
  DFN<-null
  
  nlm2<-df_  %>% 
    dplyr::group_by(uniqueID) %>%  
    dplyr::mutate(fit=list(
      tryCatch(fitSingleSigmoid(  C,  I))))
  nlm2$Tm<-Tm(nlm2)
  if(class(nlm2$fit[[1]])=='try-error'){
    return(warning("the function could not fit the data"))
  }
  #remove proteins with a plateau value not =  zero 
  CT<-nlm2 %>%
    dplyr::filter( !coef(fit[[1]])[[1]]==0)#find values where Pl = 0 
  if(nrow(CT)==0){
    return(warning("error with the plateau"))
  }
  
  #get confidence intervals for sigmoidal function
  dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95)) #predicted values
  #obtain sigmoidal equation parameters
  nlm2<- CT  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
  #ready confidence intervals for plot
  result<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                        HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                        nV=length(predict(.$fit[[1]])))# lower CI
  result$AUC<-round(pracma::trapz(result$C,result$I),2)
  
  result <-result %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]]))
  
  
  nlm1<-list()
  CT<-list()
  dfc<-list()
  
  #sigmoidal fit for treated
  nlm2<-df_1 %>% 
    dplyr::group_by(uniqueID) %>%  
    dplyr::mutate(fit=list(
      tryCatch(fitSingleSigmoid(  C,  I))))
  nlm2$Tm<-Tm(nlm2)
  
  if(class(nlm2$fit[[1]])=='try-error'){
    return(warning("the sigmoidal function could not fit the treated data"))
  }
  #remove proteins with a plateau value not =  zero 
  CT<-nlm2 %>%
    dplyr::filter( !coef(fit[[1]])[[1]]==0)#find values where Pl = 0 
  if(nrow(CT)==0){
    return(warning("error with the plateau"))
  }
  
  #get confidence intervals for sigmoidal function
  dfc <-tryCatch(nlstools::confint2(nlm2$fit[[1]],level=0.95)) #predicted values
  #obtain sigmoidal equation parameters
  nlm2<- CT  %>% dplyr::mutate(Pl=dfc[1],a=dfc[2],b=dfc[3],Pl1=dfc[4],a1=dfc[5],b1=dfc[6])
  #ready confidence intervals for plot
  result1<-  nlm2 %>%dplyr::rowwise(.) %>% dplyr::mutate(LOW = list(((1-.$Pl[1])/(1+exp(.$b[1]-(.$a[1]/.$C))))+.$Pl[1]),
                                                         HI = list(((1-.$Pl1[1])/(1+exp(.$b1[1]-(.$a1[1]/.$C))))+.$Pl1[1]),
                                                         nV=length(predict(.$fit[[1]])))# lower CI
  result1$AUC<-round(pracma::trapz(result1$C,result1$I),2)
  
  result1<-result1 %>% dplyr::rowwise() %>% dplyr::mutate(rss=deviance(.$fit[[1]]))
  
  #remove fit column
  Pred<- result %>% dplyr::select(-fit)
  Pred1<- result1 %>% dplyr::select(-fit)
  
  return(rbind(Pred,Pred1))
}

sigfit<-function(SigF,i,MD=FALSE){
  i<-i
  #
  if(isTRUE(MD)){
    Pred<-SigF[[i]] %>%
      subset(dataset=="vehicle") %>% 
      dplyr::select(uniqueID,dataset ,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,N,sample_name,sample,missing_pct,rank)
    Pred1<-SigF[[i]]%>%
      subset(dataset=="treated") %>% 
      dplyr::select(uniqueID,dataset ,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,N,sample_name,sample,missing_pct,rank)
    Pred$LOW<-Pred$LOW[[1]]
    Pred$HI<-Pred$HI[[1]]
    
    Pred1$LOW<-Pred1$LOW[[1]]
    Pred1$HI<-Pred1$HI[[1]]
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    Pred1$dAUC<-as.double(round(Pred1$AUC[1]-Pred$AUC[1],2))
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=dataset))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.55, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.35, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.25, label= paste("missing vehicle",round(Pred$missing_pct[1],0),"%"))+ 
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.15, label= paste("missing treated",round(Pred1$missing_pct[1],0),"%"))+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = 0.55,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = round(Pred$Tm[1],1), xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 )+
      annotate("text",
               x = 2+round(Pred$Tm[1],1),
               y = 0.55,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = round(min(Pred$C),1), xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2) 
    
    
    print(P1)
  }else{
    Pred<-SigF[[i]] %>%
      subset(dataset=="vehicle") %>% 
      dplyr::select(uniqueID,dataset,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,N,sample_name,sample,missing_pct,rank)
    Pred1<-SigF[[i]]%>%
      subset(dataset=="treated") %>% 
      dplyr::select(uniqueID,dataset,C,I,CC,LineRegion,Pl,a,Pl1,a1,b1,Tm,rss,
                    AUC ,LOW,HI,N,sample_name,sample,missing_pct,rank)
    Pred$LOW<-Pred$LOW[[1]]
    Pred$HI<-Pred$HI[[1]]
    
    Pred1$LOW<-Pred1$LOW[[1]]
    Pred1$HI<-Pred1$HI[[1]]
    Pred1$dTm<-round(Pred1$Tm[1]-Pred$Tm[1],1)
    Pred1$dAUC<-as.double(round(Pred1$AUC[1]-Pred$AUC[1],2))
    Pred1$RSS<-round(sum(Pred1$rss[1]+Pred$rss[1]),3)
    #Check sigmoidal fit
    P<-ggplot2::ggplot(Pred1, ggplot2::aes(x =C,y =I,color=dataset))+
      ggplot2::geom_point(Pred1,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred1,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 ) +
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.55, label= paste("\u03A3","RSS = ",Pred1$RSS[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.45, label=  paste("\u0394", "AUC = ",Pred1$dAUC[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.35, label= paste("\u0394","Tm = ",Pred1$dTm[1],"\u00B0C"))+
      ggplot2::xlab("Temperature (\u00B0C)")+ggplot2::ylab("Relative Intensity")+ ggplot2::ggtitle(paste(Pred1$uniqueID[1]))+
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.25, label= paste("missing vehicle",round(Pred1$missing_pct[1],0),"%"))+ 
      ggplot2::annotate("text", x=min(Pred1$C)+5, y=0.15, label= paste("missing treated",round(Pred1$missing_pct[1],0),"%"))+
      annotate("text",
               x = 2+round(Pred1$Tm[1],1),
               y = 0.55,
               label=paste0(round(Pred1$Tm[1],1)),
               colour="red"
      )+
      annotate("segment", x = round(Pred$Tm[1],1), xend =round(Pred1$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "red",linetype=2)+
      annotate("segment", x = round(Pred1$Tm[1],1), xend = round(Pred1$Tm[1],1), y = 0, yend = 0.5,
               colour = "red",linetype=2)
    
    
    P1<- P +ggplot2::geom_point(Pred,mapping=ggplot2::aes(x=C,y=I,color = dataset))+
      ggplot2::geom_ribbon(Pred,mapping=ggplot2::aes(ymin = LOW, ymax = HI ,fill=dataset), alpha = 0.2 )+
      annotate("text",
               x = 2+round(Pred$Tm[1],1),
               y = 0.55,
               label=paste0(round(Pred$Tm[1],1)),
               colour="blue"
      )+
      annotate("segment", x = round(min(Pred$C),1), xend =round(Pred$Tm[1],1), y = 0.5, yend = 0.5,
               colour = "blue",linetype=2)+
      annotate("segment", x = round(Pred$Tm[1],1), xend = round(Pred$Tm[1],1), y = 0, yend = 0.5,
               colour = "blue",linetype=2) 
    
    
    print(P1)
    
    
  }
  
}

###############################

plan(multisession,workers=8)

##################################

# #prepare a list of proteins
# setwd("~/Cliff prot pep")

#rename to the folder where your PSM file is located
f<- list.files(pattern='*PEPTIDES2.xlsx')
f<- list.files(pattern='*PROTEINS.xlsx')

#New
f<- list.files(pattern='*PSMs.xlsx')
f<- list.files(pattern='*_0.xlsx')

#Covid
# f<- list.files(pattern='*Proteins.xlsx')
# f<- list.files(pattern='*PSMs.xlsx')

# df.samples<-read_csv(list.files(pattern="mapping.csv")) %>% 
#   dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample"="Sample","sample_name"="MS_sample_number")
# df_raw<-df_raw %>% dplyr::left_join(df.samples,by=c("temp_ref",))


# PSMs<-read_excel(f)
# PSMs<-PSMs %>% dplyr::rename("uniqueID"="Protein","sample_name"="Mixture","dataset"="BioReplicate","temp_ref"="Channel","I"="Abundance")
# PSMs<-PSMs %>% dplyr::select(uniqueID,sample_name,dataset,temp_ref,I)
# PSMs<-PSMs %>% dplyr::mutate(sample_name=ifelse(!is.na(str_match(PSMs$sample_name,'[:punct:][:digit:][:digit:][:digit:][:punct:][:digit:]')),str_match(PSMs$sample_name,'[:punct:][:digit:][:digit:][:digit:][:punct:][:digit:]'),PSMs$sample_name)) 
# PSMs<-PSMs %>% dplyr::mutate(sample_name=str_replace(PSMs$sample_name," ","_")) %>%
#   dplyr::left_join(df.temps,by="temp_ref") %>%
#   dplyr::rename("C"="temperature") %>% 
#   dplyr::select(-temp_ref)
# d<-PSMs
# d<-d %>% dplyr::mutate(CC=ifelse(stringr::str_detect(d$sample_name,"DMSO"),0,1),
#                        rank= dplyr::ntile(I,3),
#                        sample=sample_name)
#rename to the folder where your Protein file is located
# f<-"~/Cliff prot pep/Proteins.xlsx"
# f<-"C:/Users/figue/OneDrive - Northeastern University/CETSA R/CP_Exploris_20200811_DMSOvsMEKi_carrier_FAIMS_PhiSDM_PEPTIDES.xlsx"

df_raw <- read_cetsa("~/Cliff_new","~/Cliff_new/PSMs","_0",PSM=TRUE,Batch=FALSE)
#new for PSMs
#df_raw<-df_raw %>% dplyr::rename("sample_name"="Spectrum_File")
#annotate protein data with missing values
MID<-df_raw[is.na(df_raw$value),]
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 40.1, 43.5, 47.5, 50.4, 54, 57, 60.8, 65, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = unique(df_raw$temp_ref),temperature = c(40, 42.1, 43.8, 46.5, 50, 54, 57.3, 60.1, 62, 64), stringsAsFactors = FALSE)

df.t <- function(n){
  TMT<-data.frame(NA)
  if(n==10){
    TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67), stringsAsFactors = FALSE)
  }else{
    TMT <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), temperature = c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67,68), stringsAsFactors = FALSE)
  }
  return(TMT)
}

df.temps<-df.t(10)
# df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67), stringsAsFactors = FALSE)
#df.temps <- data.frame(temp_ref = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = c(67, 64, 60.4, 57.1, 53.8, 50.5, 47.2, 43.9, 40.6, 37.3), stringsAsFactors = FALSE)
#df.samples <- data.frame(sample_id = c('F1', 'F2', 'F3','F4'), sample_name = c('DMSO_1','DMSO_2', '655_1','655_2'), stringsAsFactors = FALSE)

df.s <- function(data_path,n,rep_,bio_,vehicle_name,treated_name,Batch=FALSE,PSM=TRUE){#n is df_raw, rep is tech rep
  
  if(any(names(n)=="sample_name") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,sample_name) %>% unique(.)
    return(n)
  }
  if(any(names(n)=="Spectrum_File") & any(names(n)=="sample_id")){
    n<-n %>% dplyr::select(sample_id,Spectrum_File) %>% unique(.)
    return(n)
  }
  if(isTRUE(Batch)){
    if (!isTRUE(PSM)){
    f<-data_path
    find<-c('[:upper:][[:digit:]]+')
    check<-list()
    check<-purrr::map(seq(f),function(x){
      data.frame(sample_name= str_extract_all(names(read_xlsx(f[x]))[str_detect(names(read_xlsx(f[x])),"Abundance")],find)[[1]])
    })
    
    df.samples<-data.frame(sample_name=data.frame(f),sample_id=dplyr::bind_rows(check) )
    names(df.samples)<-c("sample_name","sample_id")
    return(df.samples)
    }
      
  }
  b<-2*rep_*bio_
  samples<-data.frame(sample_name=c(paste0(rep(as.factor(vehicle_name),rep_*bio_)),paste0(rep(as.factor(treated_name),rep_*bio_))),
                      sample_id=NA)
  samples<-rownames_to_column(samples)
  samples$sample_id<-paste0("F",samples$rowname)
  samples<-samples %>% dplyr::select(-rowname)
  samples$sample_id<-as.factor(samples$sample_id)
  samples$sample_name<-as.factor(samples$sample_name)
  return(samples)
}
df.samples<-df.s(f,df_raw,3,2,"DMSO","TREATED",Batch=TRUE)

#Convert to MSStatsTMT

editPSMs2<-data.frame()
editPSMs2<-df_raw %>% dplyr::rename("ProteinName"="Accession",
                                    "PeptideSequence"="Annotated_Sequence",
                                    "Run"="Spectrum_File",
                                    "PSM"="PSMs_Peptide_ID",
                                    "Channel"="temp_ref",
                                    "Intensity"="value")
#Condition, Bioreplicate and TechRepMixture need to be filled

editPSMs2$Condition<-ifelse(editPSMs2$Channel=="126","Norm",0)
editPSMs2$BioReplicate<-ifelse(str_detect(editPSMs2$Run,"DMSO")=="TRUE","vehicle","treated")
editPSMs2<-editPSMs2 %>% 
  dplyr::mutate(Mixture=paste0(ifelse(str_detect(Run,"NOcarrier"),"nC",ifelse(str_detect(Run,"carrier"),"C",NA)),'_',
                                   ifelse(str_detect(Run,"NO_FAIMS"),"nF",ifelse(str_detect(Run,"r_FAIMS"),"F",NA)),'_',
                                   ifelse(str_detect(Run,"S_eFT"),"E",ifelse(str_detect(Run,"S_Phi"),"S",NA))))

editPSMs2$TechRepMixture<-1
editPSMs2<-editPSMs2 %>% dplyr::select(ProteinName,PeptideSequence,Charge,PSM,Mixture,TechRepMixture,Run,Channel,Condition,BioReplicate,Intensity)

Annotation<-editPSMs2 %>% dplyr::select(Run,TechRepMixture,Channel,Condition,Mixture,BioReplicate)
Annotation$Fraction<-1
Annotation<-Annotation %>% dplyr::group_split(Run)


#method development 
####
#read filenames from separate excel sheets
###



#assign TMT channel and temperature data
df_clean <- clean_cetsa(df_raw, temperatures = df.temps, samples = df.samples,PSM=FALSE)#ssigns temperature and replicate values


df_clean <-df_clean %>% dplyr::mutate(CC=ifelse(stringr::str_detect(sample_name,"DMSO")==TRUE,0,1))#concentration values are defined in uM

df_clean$dataset<-ifelse(df_clean$CC==0,"vehicle","treated")


df_clean$sample_name<-paste0(ifelse(str_detect(df_clean$sample_name,"NOcarrier")==TRUE,"nC",ifelse(str_detect(df_clean$sample_name,"carrier")==TRUE,"C",NA)),'_',
                                                      ifelse(str_detect(df_clean$sample_name,"NO_FAIMS")==TRUE,"nF",ifelse(str_detect(df_clean$sample_name,"r_FAIMS")==TRUE,"F",NA)),'_',
                                                      ifelse(str_detect(df_clean$sample_name,"S_eFT")==TRUE,"E",ifelse(str_detect(df_clean$sample_name,"S_Phi")==TRUE,"S",NA)))
df_<-df_clean %>%
  dplyr::ungroup() %>% 
  dplyr::group_split(sample_name)#split by method development parameters before curve fitting
df_norm <- purrr::map(df_,function(x) normalize_cetsa(x, df.temps$temperature,Batch=TRUE,filters=FALSE)) #normalizes according to Franken et. al. without R-squared filter

df_norm1<-df_norm


df_v<-df_norm %>% dplyr::filter(dataset=="vehicle")
df_t<-df_norm %>% dplyr::filter(dataset=="treated")
###Generate upset plots for missing value data###
df_<-df_clean %>% dplyr::rename("sample_id"="sample")%>%
  dplyr::select(-missing_pct,-value,-missing,-rank)
df_<-df_raw%>% dplyr::right_join(df_,by=c("Accession","sample_id"))

listUP<-upMV(df_,"C_F_S",5,plot_multiple=TRUE)
pdf("UpsetMV.pdf",pointsize= 14)
listUP[[1]]

dev.off()


df_norm<-purrr::map(df_norm1,function(x)x %>% dplyr::filter(uniqueID %in% c("P36507","Q02750")))


PlotTrilinear<-function(df_norm,target){
  ##SCRIPT STARTS HERE
  DF<-df_norm %>% dplyr::group_split(uniqueID) #split null dataset only by protein ID
  d_<-df_norm %>% dplyr::filter(CC == 0) %>% dplyr::group_split(uniqueID,dataset) #split vehicle dataset
  d_1<-df_norm %>% dplyr::filter(CC > 0) %>% dplyr::group_split(uniqueID,dataset) #split treated dataset
  
  
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
  
  #preallocate list
  results<-vector(mode = "list", length(d_))
  results_t<-vector(mode = "list",length(d_1))
  results_n<-vector(mode = "list",length(DF))
  
  
  results<-suppressWarnings(DLR(d_))#First guess at line regions
  results_t<-suppressWarnings(DLR(d_1))
  results_n<-suppressWarnings(DLR(DF))
  
  
  
  #reassign shared points between line regions
  df_<-suppressWarnings(CP(results,d_))
  df_1<-suppressWarnings(CP(results_t,d_1))
  DFN<-suppressWarnings(CP(results_n,DF))# 
  #remove results to save space 
  rm(results,results_t,results_n,d_,d_1,DF)#10
  
  
  
  df_<-lapply(df_ ,function(x)x[order(x$C),])
  df_1<-lapply(df_1,function(x)x[order(x$C),])
  DFN<-lapply(DFN,function(x)x[order(x$C),])
  
  df_<-purrr::map2(df_,seq(df_),function(x,y)x %>% dplyr::mutate(N=y))
  df_1<-purrr::map2(df_1,seq(df_1),function(x,y)x %>% dplyr::mutate(N=y))
  DFN<-purrr::map2(DFN,seq(DFN),function(x,y)x %>% dplyr::mutate(N=y))
  
  #gettrilinear results
  tlresults<-list()
  tlresults_PI<-list()
  #confidence intervals
  tlresults<-tlstat(DFN,df_,df_1,norm=FALSE,Filters=FALSE,Ftest=TRUE)
  
  #return filtered lists
  res<-tlf(tlresults,DFN,APfilt=FALSE,PF=FALSE)
  i=which(res[[1]]$uniqueID %in% target)
  plotTL1<-tlCI(i,res[[1]],res[[2]],res[[3]],overlay=TRUE,residuals=FALSE)

  return(list(plotTL1))
}

plot<-purrr::map(df_norm,function(x) PlotTrilinear(x,"P36507"))

plot<-purrr::map(df_norm,function(x) PlotTrilinear(x,"Q02750"))
#saveIDs filtered
saveRDS(res[[1]]$uniqueID,"proteins_trilinear_filters_Ftest_CFE.rds")
pdf("Target_curves_trilinear_CnFE.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
plotTL1
plotTL2
dev.off()


#df1 <- only IDs in order desc(stability)
#df2<-original data in order  
#Df1 <- ordered spline results 
###############################
#get spline results
spresults<-list()
spresults_PI<-list()

spresults<-spstat(DFN,df_,df_1,Ftest=FALSE)

res_sp<-spf(spresults,DFN,filters=FALSE)
#saveIDs filtered
saveRDS(res_sp[[1]]$uniqueID,"proteins_splines_no_filters_noFtest_CFE.rds")
saveRDS(res_sp[[3]],paste0("protein_data_splines",as.character(unique(res_sp[[3]][[1]]$sample_name[[1]])),".rds"))
i=which(res_sp[[1]]$uniqueID %in% "P36507")
#generate 95%CI for splines
Pred1<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],overlay=TRUE,alpha=0.05)
i=which(res_sp[[1]]$uniqueID %in% "Q02750")
#generate 95%CI for splines
Pred2<-spCI(i,res_sp[[1]],res_sp[[2]],res_sp[[3]],overlay=TRUE,alpha=0.05)
pdf("Target_curves_splines_CnFE.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
Pred1
Pred2
dev.off()


###################################################                                                                                                                                                                                                                                                                                                                            ##################################################
#Sigmoidal function with confidence intervals
###################################################
CID<-dplyr::bind_rows(df_)$uniqueID
CID<-intersect(dplyr::bind_rows(df_1)$uniqueID,CID)
df_<-df_ %>% purrr::keep(function(x) x$uniqueID[1] %in% CID)
df_1<-df_1%>% purrr::keep(function(x) x$uniqueID[1] %in% CID)


PlSig<-purrr::map2(df_,df_1,function(x,y)try(sigC(x,y,NA)))
SigF<-PlSig %>% purrr::keep(function(x)any(class(x)=="data.frame"))

ID<-unique(dplyr::bind_rows(SigF)$uniqueID)
#saveIDs filtered
saveRDS(ID,"proteins_sigmoidal_no_filters_noFtest.rds")
i=13
sig<-sigfit(SigF,i,MD=FALSE)


