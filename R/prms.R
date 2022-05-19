
#' @title Identify scalar and vector pramaters 
#'
#' @param x dataframe of parameters
#'
#' @return Dataframe with additional columns indicating the type of parameter (scalar vs vector)
prm_scalar_vec <- function(x) {
    res = x %>% 
        mutate(isnum = grepl('^\\d\\.*e*[+-]*\\d*$', value),
               isvec = grepl('^\\d+\\.*\\d*.*;', value))
    return(res)
}



#' @title Check if prior's name exist
#' @description Checks if the name of all prior parameters are consistent 
#' with the model parameters already defined. 
#' The model parameters and the priors parameters are separately defined by the user,
#' hencethis function checks for typos and inconsistent naming.
#'
#' @param lhs Dataframe of prior samples.
#' @param prm List of model parameters.
#'
#' @return Stop the program if at least one inconsistency in the prior parameter names is found.
#'
check_prior_name_scalar <- function(lhs, prm) {
    idx.prm.scalar = which(!grepl('\\w+_\\d+', names(prm)))
    idx.lhs.scalar = which(!grepl('\\w+_\\d+', names(lhs)))
    chk.scalar.names = names(lhs)[idx.lhs.scalar] %in% names(prm)[idx.prm.scalar]
    if(!all(chk.scalar.names)) {
        ii = idx.lhs.scalar[!chk.scalar.names]
        msg = paste0('ERROR: the name of the prior parameter `',
                     names(lhs)[ii],
                     '` is not found in model parameters.')
        msg
        stop(msg)
    }   
}


#' Read parameters from a dataframe and convert into a list
read_prm_list <- function(p) {
    x = list()
    for(i in 1:nrow(p)){
        x[[i]] <- p$value[i]
        if(p$isnum[i]) x[[i]] = as.numeric(x[[i]])
        if(p$isvec[i]) 
            x[[i]] = as.numeric(unlist(strsplit(x[[i]],split = ';')))
    }
    names(x) <- p$name
    return(x)
}


digest_prm <- function(prm) {
    # prm = prm.tmp
    prm[['nE']]  = length(prm$shed.E)
    prm[['nEv']] = length(prm$shed.Ev)
    prm[['nA']]  = length(prm$shed.A)
    prm[['nH']]  = length(prm$shed.H)
    prm[['nI']]  = length(prm$shed.I)
    prm[['nIH']] = length(prm$shed.IH)
    prm[['nZ']]  = length(prm$shed.Z)
    
    # Viral shedding values were inputed as log values:
    idx.log = which(grepl('^shed', names(prm)))
    for(i in idx.log) prm[[i]] = 10^prm[[i]]
    
    return(prm)
}

#' @title Load all baseline model parameters
#' 
#' @description Load all model parameters defined in a CSV file. 
#' This is an alternative (that may be more practical) to defining the 
#' full list of parameters in a R script. 
#' Parameters are related to the simulation (population size, time horizon, etc.),
#'  the epidemic (reproduction number, etc.) and wastewater (decay, delay, etc.) 
#' The expected format of the CSV file consists of three columns named 
#' \code{name} for the parameter name, \code{value} for the corresponding value
#'  and \code{comment} for a brief explanation of the parameter.
#' 
#' @param path String. Path to the model parameters. 
#' 
#' @return A list of all model parameters.
#' 
#' @seealso The function \code{model_prm_example()} provides a set of 
#' parameters ready to use. It is a helpful template from which to start
#' customizing parameter values.
#' 
#' @export
#' 
load_prm <- function(path) {
    
    p = read.csv(path, 
                 stringsAsFactors = F, 
                 strip.white = TRUE) %>% 
        prm_scalar_vec()
    
    prm.tmp = read_prm_list(p)
    prm = digest_prm(prm.tmp)
    return(prm)
}




#' Example of model parameters
#'
#' @description This is a "helper" function that provides an 
#' example of the list of all model parameters.
#'
#' @return List of all model parameters
#' @export
#'
model_prm_example <- function() {
    
    x = list(
        dur.latent.mean = 2,  # Latent mean duration in days
        dur.latent.mean.t = NULL,  # Break time for latent mean duration in days
        dur.latent.mean.v = NULL,  # Break values for latent mean duration in days
        dur.inf.symp.mean = 12,  #Infectiousness duration for symptomatic individual in days
        dur.inf.sympHosp.mean = 8,  # Infectiousness duration for symptomatic individual before admission to hospital in days
        dur.inf.asymp.mean = 10,  # Infectiousness duration for asymptomatic individual in days
        dur.shed.recov = 24,  # fecal shedding duration after infectious period in days
        hosp.length.mean = 11,  # Length of hospital stay (obtained from DAD by Health Canada)
        hosp.length.mean.t = NULL,  # Break time for length of hospital stay
        hosp.length.mean.v = NULL,  # Break values for length of hospital stay
        dur.immunity.R = 550, # duration of full immunity for Recovered in days
        dur.immunity.R.t = NULL, # break times for time-dependent duration of full immunity for Recovered in days
        dur.immunity.R.v = NULL, # break values for time-dependent duration of full immunity for Recovered in days
        dur.immunity.V = 550, # duration of full immunity in days
        dur.immunity.V.t = NULL, # break times for time-dependent duration duration of full immunity for Vaccinated in days
        dur.immunity.V.v = NULL, # break values for time-dependent duration duration of full immunity for Vaccinated in days
        vacc.rate = 0.00001, # proportional rate of individuals get vaccinated every day (number/population)
        vacc.rate.t = NULL,  # break times for the time-dependent vaccination rate
        vacc.rate.v = NULL,   # break values for the time-dependent vaccination rate (number/population)
        vacc.eff.infection = 0.9, # vaccine effectiveness against infection given exposure 
        vacc.eff.symptomatic = 0.95, # vaccine effectiveness against symptomatic disease given exposure
        vacc.eff.hospitalization = 0.97, # vaccine effectiveness against hospitalization given exposure
        vacc.eff.t      = NULL, # break times for change in vacc. effectiveness values
        vacc.eff.inf.v  = NULL, # break values for change in vacc. effectiveness in infection given exposure 
        vacc.eff.symp.v = NULL, # break values for change in vacc. effectiveness in symptomatic given exposure 
        vacc.eff.hosp.v = NULL, # break values for change in vacc. effectiveness in hospitalization given exposure 
        dur.build.immun = 40, # days take to build immunity after receiving two-dose vaccine 
        R0 = 3.0,  # basic effective reproduction number
        transm.t =  '45 ; 55',  # break times for change in contact rate
        transm.v =  '1; 0.20' ,  # break values (multiplier for transmission rate estimated from R0) for change in contact rate
        init.I1 = 20,  # initial number of symptomatic infections introduced in the population
        init.V = 0, #initial number of vaccinated introduced in the population
        asymp.prop = 0.316,    # asymptomatic proportion
        death.prop = 0.19,  # in-hospital death rate(CIHI report https://www.cihi.ca/en/covid-19-hospitalization-and-emergency-department-statistics)
        death.prop.t = NULL,  # break time for in-hospital death rate
        death.prop.v = NULL,  # break value for in-hospital death rate
        hospital.prop = 0.02,  #proportion of hospitalized cases from symptomatic cases
        hospital.prop.t = NULL,  # break times for change in hospital proportion out of symptomatic infection
        hospital.prop.v = NULL,  # break values for change in hospital proportion out of symptomatic infection
        asymp.prop.t = NULL,  # break times for change in asymptomatic proportion 
        asymp.prop.v = NULL,  # break values for change in asymptomatic proportion (multiplier for asymp.prop)
        rel.inf.asymp = 0.8,  # relative infectiousness of asymptomatic compared to symptomatic states
        inf.A = '3;6;5;4;3;2' ,  # relative infectiousness during infectious period for asymptomatic infections (proportional to logVL)
        inf.I = '3;6;5;4;3;2' ,  # relative infectiousness during infectious period for symptomatic infections (proportional to logVL)
        inf.IH = '3;6;5;4' ,  # relative infectiousness during infectious period for symptomatic infections before admission to hospital (proportional to logVL)
        shed.E = 0.0001,  # log10 fecal shedding kinetics for exposed
        shed.Ev = 0.0001, #log10 fecal shedding kinetics for exposed
        shed.H = 0.0001,  # log10 fecal shedding kinetics for Hospital
        shed.A = '5;7.30103;7.083333;6.677121;6.30;5.9',  # log10 fecal shedding kinetics for asymptomatic
        shed.I = '5;7.30103;7.083333;6.677121;6.30;5.9',  # log10 fecal shedding kinetics for symptomatic
        shed.IH = '5;7.30103;7.083333;6.677121',  # log10 fecal shedding kinetics for symptomatic before admission to hospital
        shed.Z = '5.40103;4.60103;3.90103;3.10103;2.10103;1.2',  # log10 fecal shedding kinetics for shedding Not infectious 
        mult.shed.t= '1;50;100' ,  # time for the multiplicative factor of all fecal shedding
        mult.shed.v= '1;1;1' ,  # value for the multiplicative factor of all fecal shedding
        horizon = 200,  #horizon of the simulation
        sim.steps = 1,  #time steps per time unit
        pop.size = 50000,  #population size of the catchment area
        report.prop = 0.45,  #proportion of symptomatic reported cases (explained how drived from parameter excel sheet)
        report.lag = 10,  #lag between incident date and reported date of clinical cases in days
        episode.lag = 3,  #lag between incident date and symptom onset date for clinical cases in days
        report.lag.ww = 2,  #reporting lag between sampling date and reporting date of viral concentration in ww (in days)
        decay.rate = 0.18,  #decay rate of genetic materials of SARS-CoV-2 in ww
        transit.time.mean=1,  #mean transit time between deposited viral concentration in the sewer system and sampling sites (in days)
        transit.time.cv=0.3,  #std dev transit time between shedding and sampling sites (in days)
        ww.scale=0.0003  #scaling factor for viral concentration
    )
    
    n = length(x)
    q = data.frame(name = rep(NA,n), value = rep(NA,n))
    
    for(i in 1:n){ 
        q$name[i] = names(x)[i]
        if(is.null(x[[i]])){
            q$value[i] = 'NULL'
        }
        else{
            q$value[i] = x[[i]]
        }
    }
    
    p = prm_scalar_vec(q)
    
    prm.tmp = read_prm_list(p)
    prm = digest_prm(prm.tmp)
    
    return(prm)
}

