

#' helper function
add_obs_string <- function(x) {
    paste0(x,'.obs')
}




#' @title Create a data object from dataframes.
#'
#' @description User provides dataframes that records the number of cases, (optionally) hospital admissions and viral concentration in wastewater for a specific location.
#' The various data source do \emph{not} have to have matching date.
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
#' @param hosp.type String. Type of hospitalization provided in \code{path.hosp}: \code{NULL}, \code{'hosp.adm'} for hospital admissions, \code{'hosp.occ'} for hospital occupancy.
#'
#' @param case.date.type String. Type of date which cases are based on: \code{'report'} for reported date and \code{'episode'} for episode date 
#' 
#' @return A list of dataframes.
#' @export
#'
build_data <- function(cases, hosp, ww, hosp.type, case.date.type){
    
    # --- Checks
    check = vector()
    
    check[1] = TRUE  # if hosp.type=NULL
    check[2] = TRUE  # if case.date.type=
    if(!is.null(hosp.type)) check[1] = hosp.type %in% c('hosp.adm', 'hosp.occ')
    check[2] = case.date.type %in% c('report', 'episode')
    
    stopifnot(all(check))  
    
    # Make sure dates are in the `Date` format
    cases$date = as.Date(cases$date)
    if(!is.null(hosp)){
    hosp$date  = as.Date(hosp$date)
    }
    ww$date    = as.Date(ww$date)
    
    # "origin" date to calculate simulation time:
    d0 = as.Date(min(c(cases$date,hosp$date,ww$date), na.rm=TRUE)) - 1
    
    
    #---- Join wastewater and clinical observations
    
    # Clinical case reports only:
    dat.cl = cases %>%
        distinct() %>% 
        rename(clin.obs = value)%>%
        mutate(time = as.numeric(as.Date(date) - d0))
    
    # Hospitalizations
    if(!is.null(hosp)){
        dat.hp = hosp %>%
            distinct() %>%
            rename(hosp = value)%>%
            mutate(time = as.numeric(as.Date(date) - d0))
    }
    
    # Wastewater concentrations:
    dat.ww = ww %>%
        distinct() %>%
        rename(ww.obs = value)%>%
        mutate(time = as.numeric(as.Date(date) - d0))
    
    # Join clinical reports and wastewater
    obs.cl.ww = left_join(dat.cl, dat.ww, by='time') 
    
    # Join hospital (optional):
    if(is.null(hosp))   obs.tmp = obs.cl.ww
    if(!is.null(hosp))  obs.tmp = left_join(obs.cl.ww, dat.hp, by='time')
    
    # Final reformating:
    obs = obs.tmp %>%
        ungroup() %>%
        arrange(time) %>% 
        mutate(date = d0 + time)
    
    
    if(is.null(hosp)){
        obs = obs %>% 
            select(date, time, clin.obs, ww.obs)
    }
    
    if(!is.null(hosp)) {  
        obs = obs %>% 
            select(date, time, clin.obs, ww.obs, hosp)
        
       
        obs = obs %>%  
            rename_at(vars(matches("^hosp")), .funs = add_obs_string) %>% 
            select(date, time, clin.obs, ww.obs, starts_with('hosp'))
    }
    
    # Long format:
    
    obs.long = obs %>%
        pivot_longer(cols=-c(time, date)) %>%
        filter(!is.na(value))%>%
        mutate(type = ifelse(grepl('ww',name),'ww','clinical'))
    
    
    return(list(obs = obs, 
                obs.long = obs.long,
                hosp.var = hosp.type,
                case.var = case.date.type))
}


#' @title Create a data object from CSV files.
#'
#' @description User provides CSV files that records the number of cases, (optionally) hospital admissions and viral concentration in wastewater for a specific location.
#' The various data source do \emph{not} have to have matching date.
#' 
#' @param path.cases String. Path to the CSV file with two columns named \code{date} and \code{value} 
#' that records the number of cases with respect to time.
#' 
#' @param path.hosp String. Path to the CSV file with two columns named \code{date} and \code{value} 
#' that records the number of hospitalizations with respect to time.
#' 
#' @param path.ww String. Path to the CSV file with two columns named \code{date} and \code{value} 
#' that records the viral concentration in wastewater with respect to time.
#' 
#' @param hosp.type String. Type of hospitalization provided in \code{path.hosp}: \code{NULL}, \code{'hosp.adm'} for hospital admissions, \code{'hosp.occ'} for hospital occupancy.
#'
#' @return A list of dataframes.
#' @export
#'
build_data_csv <- function(path.cases,
                           path.hosp,
                           path.ww,
                           hosp.type) {
    
    # load csv data files 
    cases = read.csv(path.cases)
    if(!is.null(path.hosp)) hosp = read.csv(path.hosp)
    ww = read.csv(path.ww)
    
    res = build_data(cases = cases, 
                     hosp = hosp, 
                     ww = ww, 
                     hosp.type = hosp.type)
    return(res)
}
