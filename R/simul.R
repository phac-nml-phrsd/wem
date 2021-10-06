

#' Reporting of clinical cases
#' @param ts Dataframe. Time series of the epidemiological variables. 
#' @param report.prop Numerical. Proportion of cases reported (between 0 and 1).
#' @param report.lag Numerical. Lag (delay) of reporting in days. 
#' @param sim.steps Numerical. simulation time steps between every day
report_cases <- function(ts, report.prop, report.lag, sim.steps) {
    if(0){ # DEBUG
        report.prop = 0.44
        report.lag  = 3
    }
    sim.lag = sim.steps * report.lag # simulation lag for report cases 
    tmp = ts$sympinc * report.prop  
    res = c( rep(0, sim.lag-1), tmp[1:(length(tmp)-sim.lag+1)] )
    return(res)
}


#' calculate average daily deposited RNA concentration
calc_concen <- function(ts, 
                        lambda_I,
                        lambda_IH,
                        lambda_A,
                        lambda_Z,
                        mult.shed.t,
                        mult.shed.val,
                        popSize)
{
    # Identify the columns of the variables
    # associated with fecal shedding:
    indx.I  <- find_col(ts,"I")
    indx.IH <- find_col(ts,"IH")
    indx.A  <- find_col(ts,"A")
    indx.Z  <- find_col(ts,"Z")
    
    mult = sapply(ts$time, FUN = broken_line, 
                  b = mult.shed.t, v = mult.shed.val)
    
    # Calculate the product element-wise
    lambda_Is  = mapply(`*`, ts[,indx.I], lambda_I)
    lambda_IHs = mapply(`*`, ts[,indx.IH],lambda_IH)
    lambda_As  = mapply(`*`, ts[,indx.A], lambda_A)
    lambda_Zs  = mapply(`*`, ts[,indx.Z], lambda_Z)
    
    # Dataframe of all shedding contributions
    tmp = data.frame(lambda_Is,lambda_IHs,lambda_As,lambda_Zs)
    
    # Total shedding per capita
    concen = rowSums(tmp) /popSize  
     
    # Apply the multiplicative factor for variable shedding
    concen = concen * mult
    
    return(concen)
}  


#' Estimate RNA/day received at sampling site after decay and delay. 
#' Using O. Levenspiel, The Dispersion Model, pp. 47-70. 2012. 
#'
#' @param ts Dataframe. Time series from the simulation.
#' @param kappa Numeric. Daily decay rate.
#' @param sim.steps Integer. Number of time steps per unit of time.
#' @param transit.time.mean Numeric. Mean transit time of RNA particles in wastewater from shedding to sampling sites.
#' @param transit.time.cv Numeric. Coeeficient ot variation for the transit time of RNA particles in wastewater from shedding to sampling sites.
calc_delayed_concen <- function(ts, 
                                kappa, 
                                sim.steps, 
                                transit.time.mean, 
                                transit.time.cv)
{
    con   = ts$concen
    mu    = transit.time.mean
    sigma = transit.time.cv * transit.time.mean
    
    # Taking 5 sd to make sure we cover 99.99...% of the possible range
    tmax  = max(1, round(mu + 5 * sigma) )
    
    mat = matrix(0, nrow = length(con), ncol=tmax)
    tmp =  1 / sigma / sqrt(2*pi)
    for(i in 1:tmax)
    {
        # estimate transfer curve
        g = tmp * exp(-(i-mu)^2 / sigma^2 / 2)
        m = g * con * exp(-kappa*i)
        
        # apply lag in receiving m at the site
        # 'i-1' shows the lagging day
        sim.lag = (i-1) * sim.steps
        if(i==1) mat[,i] = m
        if(i>1)  mat[,i] = c(rep(0,sim.lag), m[1:(length(m)-sim.lag)])
    }
    return(mat)
}



#' turnaround time from sampling to result
calc_WWreport <- function(ts, WWreport.lag, sim.steps) 
{
    tmp = ts$samp.concen.daily
    lag = WWreport.lag * sim.steps  
    res = c( rep(0, lag), tmp[1:(length(tmp)-lag)] )
    return(res)
}

#' @title Simulate an epidemic
#' 
#' @description Simulate an epidemic by solving the system of
#' ordinary differential equations that define the SEIR-like 
#' compartmental model.
#' 
#' @param prm List of all model parameters.
#' 
#' @return A dataframe of the time series of all epidemiological variables.
#' 
#' @export
#' 
#' @examples
#' prm = model_prm_example() ; sim = simul(prm)
#' 
simul <- function(prm){
    
    # Overwrite vector elements, if specified:
    prm <- overwrite_vectors(prm, verbose = FALSE)
    
    # Ingest all model parameters:
    horizon          <- prm[["horizon"]]
    sim.steps        <- prm[["sim.steps"]]
    nE               <- prm[["nE"]]
    nEv              <- prm[["nEv"]]
    nA               <- prm[["nA"]]
    nI               <- prm[["nI"]]
    nIH              <- prm[["nIH"]]
    nZ               <- prm[["nZ"]]
    nH               <- prm[["nH"]]
    rel.inf.a        <- prm[["rel.inf.asymp"]]
    latent_mean      <- prm[["dur.latent.mean"]]
    inf_symp_mean    <- prm[["dur.inf.symp.mean"]]
    inf_sympH_mean   <- prm[["dur.inf.sympHosp.mean"]]
    inf_asymp_mean   <- prm[["dur.inf.asymp.mean"]]
    immunity         <- prm[["dur.immunity"]]
    vac.rate         <- prm[["vacc.rate"]] 
    dur.build.immun  <- prm[["dur.build.immun"]]
    eff.inf          <- prm[["vacc.eff.infection"]]
    popSize          <- prm[["pop.size"]] 
    R0               <- prm[["R0"]]
    I.init           <- prm[["init.I1"]] 
    V.init           <- prm[["init.V"]] 
    inf.A            <- prm[["inf.A"]]
    inf.I            <- prm[["inf.I"]]
    inf.IH           <- prm[["inf.IH"]]
    vload_E          <- prm[["shed.E"]]
    vload_I          <- prm[["shed.I"]]
    vload_IH         <- prm[["shed.IH"]]
    vload_A          <- prm[["shed.A"]]
    vload_H          <- prm[["shed.H"]]
    vload_Z          <- prm[["shed.Z"]] 
    mult.shed.t      <- prm[["mult.shed.t"]]
    mult.shed.val    <- prm[["mult.shed.val"]]
    alpha            <- prm[["asymp.prop"]]
    h                <- prm[["hospital.prop"]]
    alpha.vac        <- prm[["asymp.prop.vacc"]]
    h.vac            <- prm[["hospital.prop.vacc"]]
    delta            <- prm[["death.prop"]]
    shedNotInf       <- prm[["dur.shed.recov"]]
    hosp.stay        <- prm[["hosp.length.mean"]]
    kappa            <- prm[["decay.rate"]]
    report.prop      <- prm[["report.prop"]]
    report.lag       <- prm[["report.lag"]]
    WWreport.lag     <- prm[["report.lag.ww"]]
    ww.scale         <- prm[["ww.scale"]]
    transm.t         <- prm[["transm.t"]]
    transm.v         <- prm[["transm.v"]]
    vacc.rate.t      <- prm[["vacc.rate.t"]]
    vacc.rate.v      <- prm[["vacc.rate.v"]]
    hosp.rate.t      <- prm[["hosp.rate.t"]]
    hosp.rate.v      <- prm[["hosp.rate.v"]]
    hosp.rate.vacc.t <- prm[["hosp.rate.vacc.t"]]
    hosp.rate.vacc.v <- prm[["hosp.rate.vacc.v"]]
    asymp.prop.t     <- prm[["asymp.prop.t"]]
    asymp.prop.v     <- prm[["asymp.prop.v"]]
    asymp.prop.vacc.t <- prm[["asymp.prop.vacc.t"]]
    asymp.prop.vacc.v <- prm[["asymp.prop.vacc.v"]]
    transit.time.mean <- prm[['transit.time.mean']] 
    transit.time.cv   <- prm[['transit.time.cv']] 
    
    # Checks for inputs
    if(length(inf.I)!=length(vload_I)){
        stop('Inputs inconsistent: Length for `inf.I` must be the same as `shed.I`. Aborting.')
    }
    if(length(inf.A)!=length(vload_A)){
        stop('Inputs inconsistent: Length for `inf.A` must be the same as `shed.A`. Aborting.')
    }
    
    # Define simulation parameters
    n.time.steps <- horizon * sim.steps
    dt       <- seq(0, horizon, length.out = n.time.steps+1)
    epsilon  <- 1 / latent_mean
    nepsilon <- epsilon * nE
    nepsilon.vac <- epsilon * nEv
    tau      <- 1 / inf_symp_mean
    ntau     <- tau * nI
    tau.immu <- 1 / immunity
    mu       <- 1 / inf_sympH_mean
    nmu      <- mu * nIH
    theta    <- 1/inf_asymp_mean
    ntheta   <- theta * nA
    ell      <- 1/ hosp.stay
    nell     <- ell * nH
    eta      <- 1/shedNotInf
    neta     <- eta * nZ
    
    # define vaccine parameters 
    r <- vac.rate         
    d <- 1 / dur.build.immun 
    
    
    # RNA copies concentration
    lambda_E  <- vload_E 
    lambda_I  <- vload_I
    lambda_IH <- vload_IH
    lambda_A  <- rel.inf.a * vload_A
    lambda_H  <- vload_H
    lambda_Z  <- vload_Z
    
    #----- beta from R0
    S0 = popSize - V.init - I.init
    V0 = V.init
        
    # Unvaccinatd part
    adj.J  = sum(inf.I[1:nIH]) / sum(inf.I)
    
    sA  =  alpha / theta * rel.inf.a
    sI  =  (1-h) * (1 - alpha) / tau
    sJ  =      h * (1 - alpha) / mu
    s.S   = (sA + sI + sJ * adj.J)*S0/popSize  
    
    # Vaccinated part 
    sA.V  =  alpha.vac / theta * rel.inf.a
    sI.V  =  (1-h.vac) * (1 - alpha.vac) / tau
    sJ.V  =      h.vac * (1 - alpha.vac) / mu
    s.V   = (sA.V + sI.V + sJ.V * adj.J) * (1 - eff.inf)*V0/popSize
    
    # Overall R0
    s = s.S +s.V
    beta = R0 / s
    
    params.SEIR <- list(  
        h = h,
        alpha = alpha,
        h.vac = h.vac,
        alpha.vac = alpha.vac,
        delta = delta,
        beta = beta,
        nepsilon = nepsilon,
        nepsilon.vac = nepsilon.vac,
        ntau = ntau,
        nmu = nmu,
        ntheta = ntheta,
        neta = neta,
        nell = nell,
        nE=nE, nEv=nEv, nI=nI, nIH=nIH, 
        nA=nA, nH=nH, nZ=nZ,
        popSize = popSize,
        transm.t = transm.t, 
        transm.v = transm.v,
        vacc.rate.t = vacc.rate.t, 
        vacc.rate.v = vacc.rate.v,
        hosp.rate.t = hosp.rate.t, 
        hosp.rate.v = hosp.rate.v,
        asymp.prop.t = asymp.prop.t, 
        asymp.prop.v = asymp.prop.v,
        hosp.rate.vacc.t = hosp.rate.vacc.t, 
        hosp.rate.vacc.v = hosp.rate.vacc.v,
        asymp.prop.vacc.t = asymp.prop.vacc.t, 
        asymp.prop.vacc.v = asymp.prop.vacc.v,
        inf.A = inf.A, 
        inf.I = inf.I,
        inf.IH = inf.IH,
        rel.inf.a = rel.inf.a,
        eff.inf = eff.inf,
        r=r, d=d, tau.immu=tau.immu)
    
    ### Initial conditions
    ###
    ### * * * WARNING * * *
    ### MUST BE IN THE SAME ORDER AS 
    ### THE VARIABLES IN THE ODE DEFINED 
    ### BY THE `seir` FUNCTION
    
    EEvIAHZ_vec <- c(
        E  = rep(0, ifelse(nE==Inf,0,nE)),
        Ev = rep(0, ifelse(nEv==Inf,0,nEv)),
        I  = c(I.init, rep(0, ifelse(nI==Inf, 0, nI-1))),
        IH = rep(0, ifelse(nIH==Inf,0,nIH)),
        A  = rep(0, ifelse(nA==Inf,0,nA)),
        H  = rep(0, ifelse(nH==Inf,0,nH)),
        Z  = rep(0, ifelse(nZ==Inf,0,nZ)))
    
    inits.SEIR <- c(S = popSize -V.init -I.init,
                    Vw = 0,
                    V  = V.init,
                    EEvIAHZ_vec,
                    R=0, D=0, 
                    cuminc = I.init, 
                    cumincsymp = I.init, 
                    cumHospAdm = 0)
    
    #### Simulation (i.e., solutions of the ODEs)
    ts <- as.data.frame(
        lsoda(y      = inits.SEIR, 
              times  = dt, 
              func   = seir, 
              parms  = params.SEIR))
    
    ts$Eall   = calc.all(ts,"E")  # <- all exposed indiv.
    ts$Evall  = calc.all(ts,"Ev") # <- all vaccinated exposed indiv.
    ts$Iall   = calc.all(ts,"I")  # <- all symptomatic indiv.
    ts$IHall  = calc.all(ts,"IH") # <- all symptomatic goes to hospital 
    ts$Aall   = calc.all(ts,"A")  # <- all asymptomatic indiv.
    ts$Hall   = calc.all(ts,"H")  # <- all hospital occupancy
    ts$Zall   = calc.all(ts,"Z")  # <- all recovered yet shedding indiv.
    ts$prev   = ts$Aall + ts$Iall + ts$IHall + ts$Hall   # <- global prevalence
    ts$inc    = c(I.init, diff(ts$cuminc))               # <- global incidence
    ts$sympinc  = c(I.init, diff(ts$cumincsymp))         # <- symptomatic incidence
    ts$hosp.admission = c(0, diff(ts$cumHospAdm))        # <- hospital admissions
    
    ts$report = report_cases(ts, report.prop, report.lag, sim.steps)
    
    # Calculate the daily concentration 
    # deposited in the sewer system:
    ts$concen   = calc_concen(ts,
                              lambda_I,lambda_IH,
                              lambda_A,lambda_Z,
                              mult.shed.t,
                              mult.shed.val,
                              popSize)
    
    # Calculate the concentration arriving 
    # at the sampling site after decay and delay:
    samp.concen = calc_delayed_concen(ts,
                                      kappa,
                                      sim.steps,
                                      transit.time.mean, 
                                      transit.time.cv)
    
    # ww.scale is a constant for uncertainties 
    # involves in RNA transfer, dilution, hydraulic degradation
    # and experimental extraction/quantification 
    samp.concen.daily = rowSums(samp.concen) * ww.scale
    
    # Add wastewater sampling variables to the human epi dataframe:
    ts = data.frame(ts, samp.concen.daily)
    
    # Wastewater concentration reported (ie, including reporting delay):
    ts$WWreport = calc_WWreport(ts, WWreport.lag, sim.steps)
    
    return(list(ts=ts, R0=R0))
}




simul_post_unit <- function(i, post.abc, prm) {
    
    np = ncol(post.abc)
    for(k in 1:np){
        prm[[names(post.abc)[k]]] <- post.abc[i,k]
    }
    sim = simul(prm)
    df  = sim$ts %>% 
        mutate(iter = i)
    return(df)
}


#' Simulate epidemic from posterior distribution
#' 
#' @param post.abc Dataframe of posterior values.
#' @param prm List of (not fitted) model parameters. 
#' @param ci Numeric. Credible interval.
#' @param n.cores Integer. Number of cores for parallel computing
#' 
#' @return Dataframe of summary statistics time series.
#' 
#' @export
#' 
simul_from_post <- function(post.abc, prm, ci=0.95, n.cores = 4) {
    
    n.abc.post = nrow(post.abc)
    message(paste('Starting simulation using',n.abc.post,'posterior values...'))
    
    sfInit(parallel = n.cores > 1, cpus = n.cores)
    sfExportAll()
    suppressMessages({
        sfLibrary(deSolve)
        sfLibrary(stringr)
        sfLibrary(dplyr)
    })
    tmp = sfLapply(x   = 1:n.abc.post, 
                   fun = simul_post_unit, 
                   post.abc = post.abc, prm = prm)
    sfStop()
    
    simp = do.call('rbind', tmp)
    
    ss = simp %>%
        group_by(time) %>% 
        summarise(prev.m  = mean(prev),
                  prev.lo = quantile(prev, probs = 0.5 - ci/2),
                  prev.hi = quantile(prev, probs = 0.5 + ci/2),
                  inc.m   = mean(inc),
                  inc.lo  = quantile(inc, probs = 0.5 - ci/2),
                  inc.hi  = quantile(inc, probs = 0.5 + ci/2),
                  hosp.adm.m  = mean(hosp.admission),
                  hosp.adm.lo = quantile(hosp.admission, probs = 0.5 - ci/2),
                  hosp.adm.hi = quantile(hosp.admission, probs = 0.5 + ci/2),
                  report.m  = mean(report),
                  report.lo = quantile(report, probs = 0.5 - ci/2),
                  report.hi = quantile(report, probs = 0.5 + ci/2),
                  ww.m  = mean(WWreport),
                  ww.lo = quantile(WWreport, probs = 0.5 - ci/2),
                  ww.hi = quantile(WWreport, probs = 0.5 + ci/2),
                  S.m   = mean(S),
                  S.lo  = quantile(S, probs = 0.5 - ci/2),
                  S.hi  = quantile(S, probs = 0.5 + ci/2)) %>% 
        # Calculate cumulative incidence
        mutate(cuminc.m  = 1 - S.m/S.m[1],
               cuminc.lo = 1 - S.hi/S.hi[1],
               cuminc.hi = 1 - S.lo/S.lo[1])
    
    return(ss)
}




#' Add an observation process distributed as a negative binomial.
#'
#' @param x Numerical vector. Mean of observations.
#' @param disp Numerical. Dispersion parameter (the larger, the less dispersed)
#'
#' @return A vector the same length as x where each element is a 
#' random sample with mean x[i].
#' 
obs_nbinom <- function(x, disp) {
    # Note: variance = m + m^2/disp
    # where m is the mean.
    rnbinom(n=length(x), size = disp, mu = x)
    return(y)
}

#' Add an observation process distributed as normal.
#'
#' @param x Numerical vector. Mean of observations.
#' @param cv Numerical. Coefficient of variation.
#' @param do.round Logical. Round random samples (default = FALSE).
#'
#' @return A vector the same length as x where each element is a 
#' random sample with mean x[i].
#' 
obs_norm <- function(x, cv, do.round = FALSE) {
    # Note: CV = stddev / mean
    y = rnorm(n=length(x), mean = x, sd = cv*x)
    if(do.round) y = round(y)
    return(y)
}

draw_obs <- function(x, prms) {
    
    if(prms$distribution == 'norm'){
        
        do.round = FALSE
        if(!is.null(prms$do.round)){
            do.round = prms$do.round
        }
        y = obs_norm(x, cv = prms$cv, do.round = do.round)
    }
    if(prms$distribution == 'nbinom'){
        y = obs_binom(x, disp = prms$disp)
    }
    
    res = y 
    
    if(prms$freq < 1){
        # draw the times of observations
        idx = rbinom(n = length(x), size = 1, prob = prms$freq)
        idx
        res[!idx] <- NA
        res
    }
    return(res)
}

