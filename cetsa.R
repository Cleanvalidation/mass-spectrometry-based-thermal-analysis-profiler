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
find_pat <- function(pat, x) 
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


