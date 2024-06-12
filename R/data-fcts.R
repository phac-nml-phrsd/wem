#' @title Add '.obs' suffix to strings.
#' 
#' @description Function appends strings with '.obs'.
#' 
#' @param x String. 
#' 
#' @return String object with '.obs' suffix.
add_obs_string <- function(x) {
  paste0(x,'.obs')
}




#' @title Create a data object from dataframes.
#'
#' @description User provides dataframes that records the number of cases,
#'  (optionally) hospital admissions and viral concentration in wastewater for a
#'  specific location. The various data source do \emph{not} have to have
#'  matching date.
#' 
#' @param cases Dataframe with two columns named \code{date} and \code{value} 
#' that records the number of cases with respect to time.
#' 
#' @param ww Dataframe with two columns named \code{date} and \code{value} 
#' that records the viral concentration in wastewater with respect to time.
#' 
#' @param hosp Dataframe with two columns named \code{date} and \code{value} 
#' that records the number of hospitalizations with respect to time.
#' 
#' @param hosp.type String. Type of hospitalization provided in \code{path.hosp}:
#'  \itemize{
#'   \item \code{NULL},
#'   \item \code{'hosp.adm'} for hospital admissions,
#'   \item \code{'hosp.occ'} for hospital occupancy.
#'  } 
#'
#' @param case.date.type String. Type of date which cases are based on:
#'  \itemize{
#'   \item \code{'report'} for reported date 
#'   \item \code{'episode'} for episode date (date of symptoms onset). 
#'  }
#'  
#' 
#' @return A list of dataframes.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' 
#' # Load data sets examples from `wem` package:
#' data('cases')
#' data('hosp')
#' data('wwviralconc')
#' 
#' # Build the data object:
#' dat = build_data(cases = cases, 
#'                  hosp = hosp, 
#'                  ww = wwviralconc, 
#'                  hosp.type = 'hosp.adm', 
#'                  case.date.type = 'report')
#'
build_data <- function(cases, hosp, ww, hosp.type, case.date.type){
  
  # --- Checks
  
  hosp.types      = c('hosp.adm', 'hosp.occ')
  case.date.types = c('report', 'episode')
  
  if(!is.null(hosp.type)) 
    check.h = hosp.type %in% hosp.types
  
  if(!check.h) stop(paste('parameter `host.type` must be NULL (:ignored) or',
                          'take these values:',
                          paste(hosp.types, collapse=', ')))
  
  check.c = case.date.type %in% case.date.types
  
  if(!check.c) stop(paste('parameter `case.date.type` must take these values:',
                          paste(case.date.types, collapse=', ')))
  
  
  # Make sure dates are in the `Date` format
  cases$date = as.Date(cases$date)
  ww$date    = as.Date(ww$date)
  if(!is.null(hosp)) { hosp$date = as.Date(hosp$date) }
  
  # "origin" date to calculate simulation time:
  date.vec = c(cases$date, ww$date)
  if(!is.null(hosp)) date.vec = c(date.vec, hosp$date)
  d0 = as.Date(min(date.vec, na.rm=TRUE)) - 1
  
  #---- Join all data sources
  
  # Clinical case reports only:
  dat.cl = cases %>%
    distinct() %>% 
    rename(clin.obs = value)
  
  # Wastewater concentrations:
  dat.ww = ww %>%
    distinct() %>%
    rename(ww.obs = value)
  
  # Hospitalizations
  if(!is.null(hosp)){
    dat.hp = hosp %>%
      distinct() %>%
      rename(hosp = value)
  }
  
  dtrng = range(date.vec, na.rm = TRUE)
  tmp = data.frame(date = seq.Date(dtrng[1], dtrng[2], by=1))
  
  tmp2 = tmp |> 
    left_join(dat.cl, by = 'date') |>
    left_join(dat.ww, by = 'date')
  
  if(!is.null(hosp)) {
    tmp2 = tmp2 |> 
      left_join(dat.hp, by = 'date') |> 
      rename(hosp.obs = hosp)
  }
  if(is.null(hosp)) tmp2$hosp.obs = NA
  
  obs = tmp2 |> 
    # Remove dates with no observations at all
    mutate(chk = is.na(clin.obs) & 
             is.na(ww.obs) & 
             is.na(hosp.obs)) |>
    filter(!chk) |> 
    select(-chk) |>
    # Calculate times:
    mutate(time = as.numeric(as.Date(date) - d0)) |> 
    arrange(date) |> 
    ungroup()
  
  # Long format:
  obs.long = obs %>%
    pivot_longer(cols=-c(time, date)) %>%
    filter(!is.na(value)) %>%
    mutate(type = ifelse(grepl('ww',name),'ww','clinical'))
  
  return(list(obs = obs, 
              obs.long = obs.long,
              hosp.var = hosp.type,
              case.var = case.date.type))
}


#' @title Create a data object from CSV files.
#'
#' @description User provides CSV files that records the number of cases,
#'  (optionally) hospital admissions and viral concentration in wastewater for a
#'  specific location. The various data source do \emph{not} have to have
#'  matching date.
#' 
#' @param path.cases String. Path to the CSV file with two columns named
#'  \code{date} and \code{value} that records the number of cases with respect
#'  to time.
#' 
#' @param path.hosp String. Path to the CSV file with two columns named
#'  \code{date} and \code{value} that records the number of hospitalizations
#'  with respect to time.
#' 
#' @param path.ww String. Path to the CSV file with two columns named
#'  \code{date} and \code{value} that records the viral concentration in
#'  wastewater with respect to time.
#' 
#' @param hosp.type String. Type of hospitalization provided in \code{path.hosp}:
#'  \itemify{
#'   \item \code{NULL}
#'   \item \code{'hosp.adm'} for hospital admissions
#'   \item \code{'hosp.occ'} for hospital occupancy.
#'  }
#'   
#' @param case.date.type String. Type of date which cases are based on:
#'  \itemify{
#'   \item \code{'report'} for reported date
#'   \item \code{'episode'} for episode date (date of symptoms onset).
#'  }
#' 
#' @return A list of dataframes.
#' @export
#'
build_data_csv <- function(path.cases,
                           path.hosp,
                           path.ww,
                           hosp.type,
                           case.date.type) {
  
  # load csv data files 
  cases = read.csv(path.cases)
  if(!is.null(path.hosp)){
    hosp = read.csv(path.hosp)  
  }else{
    hosp = NULL
  } 
  ww = read.csv(path.ww)
  
  res = build_data(cases = cases, 
                   hosp = hosp, 
                   ww = ww, 
                   hosp.type = hosp.type,
                   case.date.type = case.date.type)
  return(res)
}
