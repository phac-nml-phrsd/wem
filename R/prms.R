
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

#' @title Load all baseline parameters for the wastewater epidemic model.
#' 
#' @description Load all model parameters defined in a CSV file.
#' Parameters are related to the simulation (population size, time horizon, etc.),
#'  the epidemic (reproduction number, etc.) and wastewater (decay, delay, etc.) 
#' The expected format of the CSV file consists of three columns named 
#' \code{name} for the parameter name, \code{value} for the corresponding value
#'  and \code{comment} for a brief explanation of the parameter.
#' 
#' 
#' @param path String. Path to the model parameters. 
#' 
#' @return A list of all parameters.
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
        dur.inf.symp.mean = 12,  #Infectiousness duration for symptomatic individual in days
        dur.inf.sympHosp.mean = 8,  # Infectiousness duration for symptomatic individual before admission to hospital in days
        dur.inf.asymp.mean = 10,  # Infectiousness duration for asymptomatic individual in days
        dur.shed.recov = 24,  # fecal sheding duration after infectious period in days
        hosp.length.mean = 11,  # Length of hospital stay for AB (obtained from DAD by Health Canada)
        dur.immunity = 365, # duration of full immunity in days
        vacc.rate = 0.00001, # number of indiviudals get vaccinated every day
        vacc.eff.infection = 0.9999, # vaccine effectiveness against infection
        dur.build.immun = 40, # days take to build immunity after reciveing two-dose vaccine 
        R0 = 3.0,  # basic effective reproduction number
        init.I1 = 20,  # initial number of symptomatic infections introduced in the population
        init.V = 0, #initial number of vaccinated introduced in the population
        asymp.prop = 0.316,    # asymptomatic proportion
        asymp.prop.vacc = 0.70, # asymptomatic proportion for vaccinated individuals
        death.prop = 0.19,  # proportion of death from hospitalized (CIHI report https://www.cihi.ca/en/covid-19-hospitalization-and-emergency-department-statistics)
        hospital.prop = 0.02,  #proportion of hospitalized cases from symptomatic cases
        hospital.prop.vacc = 0.01, #proportion of hospitalized cases from symptomatic cases
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
        mult.shed.val= '1;1;1' ,  # value for the multiplicative factor of all fecal shedding
        transm.t =  '45 ; 55',  # break points dates for change in contact rate
        transm.v =  '1; 0.20' ,  # break points values for change in contact rate
        hosp.rate.t = '30; 40',    # break times for change in hospitalization rate 
        hosp.rate.v = '1.0; 1.0',  # break values for change in hospitalization rate 
        vacc.rate.t = '1 ; 10',  # break times for the time-dependent vaccination rate
        vacc.rate.v = '0 ; 0',   # break values for the time-dependent vaccination rate
        asymp.prop.t = '30; 40',  # break times for change in asymptomatic proportion
        asymp.prop.v = '1.0; 1.0',  # break values for change in asymptomatic proportion
        hosp.rate.vacc.t =  '30; 40',#break times for change in hospital rate (dates between waves) 
        hosp.rate.vacc.v = '1; 1.001', #break values for change in hospital rate (value for each wave) 
        asymp.prop.vacc.t = '30; 40', #break times for change in asymptomatic proportion
        asymp.prop.vacc.v = '1.0; 1.0', #break values for change in asymptomatic proportion
        horizon = 200,  #?horizon of the simulation
        sim.steps = 1,  #?time steps per time unit
        pop.size = 50000,  #?population size ofthe catchment area
        report.prop = 0.45,  #?proportion of symptomatic reported cases (explained how drived from parameter excel sheet)
        report.lag = 10,  #?lag between infection and report in days
        report.lag.ww = 2,  #?reporting lag between sampling date and reporting date in days
        decay.rate = 0.18,  #?decay rate of RNA in ww
        transit.time.mean=1,  #?mean transit time between shedding and sampling sites (in days)
        transit.time.cv=0.3,  #?std dev transit time between shedding and sampling sites (in days)
        ww.scale=0.0003  #?scaling factor for viral concentration
    )
    
    n = length(x)
    q = data.frame(name = rep(NA,n), value = rep(NA,n))
    
    for(i in 1:n){ 
        q[i,1] = names(x)[i]
        q[i,2] = x[[i]]
    }
    
    p = prm_scalar_vec(q)
    
    prm.tmp = read_prm_list(p)
    prm = digest_prm(prm.tmp)
    
    return(prm)
}

