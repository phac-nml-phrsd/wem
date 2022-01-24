

#' Reporting of clinical cases
#' @param ts Dataframe. Time series of the epidemiological variables. 
#' @param report.prop Numerical. Proportion of cases reported (between 0 and 1).
#' @param lag Numerical. Lag (delay) btw. incident date and reported date or episode date in days. 
#' @param sim.steps Numerical. simulation time steps between every day
report_cases <- function(ts, report.prop, lag, sim.steps) {
    if(0){ # DEBUG
        report.prop = 0.44
        lag  = 3
    }
    sim.lag = round(sim.steps * lag) # simulation lag for report cases 
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
                        mult.shed.v,
                        popSize)
{
    # Identify the columns of the variables
    # associated with fecal shedding:
    indx.I  <- find_col(ts,"I")
    indx.IH <- find_col(ts,"IH")
    indx.A  <- find_col(ts,"A")
    indx.Z  <- find_col(ts,"Z")
    
    mult = sapply(ts$time, FUN = broken_line, 
                  b = mult.shed.t, v = mult.shed.v)
    
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
#' @param transit.time.cv Numeric. Coefficient ot variation for the transit time of RNA particles in wastewater from shedding to sampling sites.
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
        sim.lag = round((i-1) * sim.steps)
        if(i==1) mat[,i] = m
        if(i>1)  {
          foo = c(rep(0,sim.lag), m[1:(length(m)-sim.lag)])
          mat[,i] = foo
        }
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
    eff.symp         <- prm[["vacc.eff.symptomatic"]]
    eff.hosp         <- prm[["vacc.eff.hospitalization"]]
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
    mult.shed.v      <- prm[["mult.shed.v"]]
    alpha            <- prm[["asymp.prop"]]
    hosp.prop        <- prm[["hospital.prop"]]
    delta            <- prm[["death.prop"]]
    shedNotInf       <- prm[["dur.shed.recov"]]
    hosp.stay        <- prm[["hosp.length.mean"]]
    kappa            <- prm[["decay.rate"]]
    report.prop      <- prm[["report.prop"]]
    report.lag       <- prm[["report.lag"]]
    episode.lag      <- prm[["episode.lag"]]
    WWreport.lag     <- prm[["report.lag.ww"]]
    ww.scale         <- prm[["ww.scale"]]
    transm.t         <- prm[["transm.t"]]
    transm.v         <- prm[["transm.v"]]
    vacc.rate.t      <- prm[["vacc.rate.t"]]
    vacc.rate.v      <- prm[["vacc.rate.v"]]
    hosp.prop.t      <- prm[["hospital.prop.t"]]
    hosp.prop.v      <- prm[["hospital.prop.v"]]
    asymp.prop.t     <- prm[["asymp.prop.t"]]
    asymp.prop.v     <- prm[["asymp.prop.v"]]
    
    # the  time vector `eff.t` (below) applies to 
    # all other time dependent parameters
    # associated with vaccine effectiveness:
    eff.t            <- prm[["vacc.eff.t"]]    
    
    eff.inf.v        <- prm[["vacc.eff.inf.v"]]
    eff.symp.v       <- prm[["vacc.eff.symp.v"]]
    eff.hosp.v       <- prm[["vacc.eff.hosp.v"]]
    
    transit.time.mean <- prm[['transit.time.mean']] 
    transit.time.cv   <- prm[['transit.time.cv']] 
    
    # --- Checks for inputs
    
    if(length(inf.I)!=length(vload_I)){
        stop('Inputs inconsistent: Length for `inf.I` must be the same as `shed.I`. Aborting.')
    }
    if(length(inf.A)!=length(vload_A)){
        stop('Inputs inconsistent: Length for `inf.A` must be the same as `shed.A`. Aborting.')
    }
    
    # vaccine effectiveness inputs
    if(length(eff.inf.v)!=length(eff.symp.v)){
      stop('Inputs inconsistent: Length for `vacc.eff.inf.v` must be the same as `vacc.eff.symp.v`. Aborting.')
    }
    if(length(eff.hosp.v)!=length(eff.symp.v)){
      stop('Inputs inconsistent: Length for `vacc.eff.hosp.v` must be the same as `vacc.eff.symp.v`. Aborting.')
    }
    if(length(eff.t)!=length(eff.symp.v)){
      stop('Inputs inconsistent: Length for `vacc.eff.t` must be the same as `vacc.eff.symp.v`. Aborting.')
    }
    if(length(eff.t)==1 & eff.t!='NULL'){
      stop('ERROR: Length for `eff.t` must be of size 2 or more. For a constant values for vaccine effectiveness, use those parameters instead: `vacc.eff.infection`, `vacc.eff.symptomatic`, `vacc.eff.hospitalization` . Aborting.')
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
    
    #--- Define vaccine parameters before passing for simulation
    
    r <- vac.rate         
    d <- 1 / dur.build.immun 
    
    
    # Check if the input are time-dependent
    is.time.dep.vacc.eff = length(eff.t)>1
    
    if(!is.time.dep.vacc.eff){
      #message('Vacine effectiveness parameters are constant.')
      
      eff.symp.inf  = 1 - (1-eff.symp)/(1-eff.inf)
      eff.hosp.symp = 1 - (1-eff.hosp)/(1-eff.symp)
      alpha.vac     = eff.symp.inf
      h.vac         = 1 - eff.hosp.symp 
      
      asymp.prop.vacc.v = NULL
      hosp.rate.vacc.v  = NULL
      asymp.prop.vacc.t = NULL
      hosp.rate.vacc.t  = NULL
    }
    
    if(is.time.dep.vacc.eff){
      
      message('Vacine effectiveness parameters are time-dependent.')
      
      # constant vacc. variables
      # even in case of using time-dependent h.vac and alpha.vac, 
      # we still use const. alpha.vac and h.vac for R0 and beta 
      eff.symp.inf  = 1 - (1 - eff.symp.v[1])/(1 - eff.inf.v[1])
      eff.hosp.symp = 1 - (1 - eff.hosp.v[1])/(1 - eff.symp.v[1])
      alpha.vac     = eff.symp.inf
      h.vac         = 1 - eff.hosp.symp 
      
      # time-dependent vacc. variables
      asymp.prop.vacc.v = 1 - (1-eff.symp.v)/(1-eff.inf.v)
      hosp.rate.vacc.v  = 1 - (1-eff.hosp.v)/(1-eff.symp.v) 
      asymp.prop.vacc.t = eff.t
      hosp.rate.vacc.t  = eff.t
    }
    
    #--- RNA copies concentration
    lambda_E  <- vload_E 
    lambda_I  <- vload_I
    lambda_IH <- vload_IH
    lambda_A  <- rel.inf.a * vload_A
    lambda_H  <- vload_H
    lambda_Z  <- vload_Z
    
    #----- beta from R0
    S0 = popSize - V.init - I.init
    V0 = V.init
        
    # Unvaccinated part
    adj.J  = sum(inf.I[1:nIH]) / sum(inf.I)
    
    sA  =  alpha / theta * rel.inf.a
    sI  =  (1-hosp.prop) * (1 - alpha) / tau
    sJ  =      hosp.prop * (1 - alpha) / mu
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
        hosp.prop = hosp.prop,
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
        hosp.prop.t = hosp.prop.t, 
        hosp.prop.v = hosp.prop.v,
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
      deSolve::lsode(
        y      = inits.SEIR,
        times  = dt,
        func   = seir,
        parms  = params.SEIR,
        mf     = 10,
        rtol   = 1e-2,
        atol   = 1e-2)
    )
    
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
    
    ts$report         = report_cases(ts, report.prop, lag = report.lag, sim.steps)
    ts$report.episode = report_cases(ts, report.prop, lag = episode.lag, sim.steps)
    
    # Calculate the daily concentration 
    # deposited in the sewer system:
    ts$concen   = calc_concen(ts,
                              lambda_I,lambda_IH,
                              lambda_A,lambda_Z,
                              mult.shed.t,
                              mult.shed.v,
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


#' @title Simulate epidemic from posterior distributions.
#' 
#' @description Use posterior distributions of a fitted object to simulate an epidemic.
#' 
#' @param post.abc Dataframe of posterior values.
#' @param prm List of (not fitted) model parameters.
#' @param hosp.var String. Type of hospitalization; \code{NULL}, \code{'hosp.adm'} and \code{'hosp.occ'}
#' @param case.var String. Type of date for clinical cases; \code{'report'} and \code{'episode'}
#' @param ci Numeric. Credible interval level. Default = 0.95.
#' @param n.cores Integer. Number of cores used for parallel computing.
#' 
#' @return Dataframe of summary statistics time series.
#' 
#' @export
#' 
simul_from_post <- function(post.abc, prm, 
                            hosp.var, case.var, 
                            ci = 0.95, n.cores = 4) {
  
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
  
  # Determine type of date for clinical cases (reported date or episode date)
  if(case.var=='report'){
    simp = simp %>%
      mutate(clin.case = report)
  }
  if(case.var=='episode'){
    simp = simp %>%
      mutate(clin.case = report.episode)
  }
  
  if(!is.null(hosp.var)){
    
    # Determine hospital type (new admissions or occupancy)
    if(hosp.var=='hosp.adm'){
      simp = simp %>%
        mutate(hospital = hosp.admission)
    }
    if(hosp.var=='hosp.occ'){
      simp = simp %>%
        mutate(hospital = Hall)
    }
  }  
  #---- Summary stats of posterior simulations:
  
  varnames = c('prev', 'inc',  'clin.case', 
               'WWreport', 'S', 'Aall', 'V')
  
  if(!is.null(hosp.var)) varnames = c(varnames, 'hospital') 
  
  ss = simp %>% 
    group_by(time) %>%
    # Calculate mean and quantiles for 
    # each variables defined in `vv` :
    summarise(
      across(.cols = varnames, 
             .fns = list(m  = mean, 
                         lo = ~quantile(.x, probs= 0.5 - ci/2),
                         hi = ~quantile(.x, probs= 0.5 + ci/2)), 
             .names = "{.col}.{.fn}")) %>% 
    #
    # Calculate cumulative incidence
    mutate(cuminc.m  = 1 - S.m/S.m[1],
           cuminc.lo = 1 - S.hi/S.hi[1],
           cuminc.hi = 1 - S.lo/S.lo[1]) %>%
    #
    # Rename variables.  
    # TODO: do something smarter! 
    # Maybe do not rename and use original var names, 
    # as this would be more informative and consistent (?)
    #
    rename(ww.lo = WWreport.lo, 
           ww.m  = WWreport.m, 
           ww.hi = WWreport.hi,
           report.lo = clin.case.lo,
           report.m  = clin.case.m,
           report.hi = clin.case.hi,
           asymp.lo = Aall.lo, 
           asymp.m  = Aall.m, 
           asymp.hi = Aall.hi)
  
  if(!is.null(hosp.var)){
    ss = rename(ss, 
                hosp.lo = hospital.lo, 
                hosp.m  = hospital.m, 
                hosp.hi = hospital.hi)
  }
  
  return(ss)
}



#' @title Generate observation noise
#' 
#' @description Generate observation noise on the reported values for clinical cases
#' (variable `report`) and concentration in wastewater (variable `WWreport`). 
#' The observations are generated by drawing samples from a negative binomial distribution
#' with mean the value from the deterministic model and a variance parametrized
#' using the coefficient of variation (std. dev / mean). 
#'
#' @param df Dataframe. The time series of the epidemiological variables generated
#' by the function \code{wem::simul(...)$ts}.
#' @param prms List of parameters. Must contain elements named \code{cv} and \code{cv.ww} 
#' to inform the desired coefficient of variation of the added noise for clinical 
#' reports and concentration in wastewater, respectively. 
#'
#' @return The input dataframe with three additional variables, 
#' \code{report.obs}, \code{hosp.admission.obs} and \code{WWreport.obs}, representing the 
#' simulated observations of the variables  \code{report}, \code{hosp.admission}  and \code{WWreport}.
#' 
#' @export
#'
#' @examples 
#' prm = model_prm_example() 
#' sim = simul(prm)
#' df = generate_obs_noise(sim$ts) 
#' plot(x=df$time, y=df$report.obs) ; lines(x=df$time, y=df$report)
#' 
generate_obs_noise <- function(df, prms = list(cv = 0.1, cv.ww = 0.1)) {
    
    mu    = sapply(df$report, mean)    # for clinical cases
    mu.ww = sapply(df$WWreport, mean)  # for ww concentration
    mu.ha = sapply(df$hosp.admission, mean)  # for hospital admissions
    
    # See `rbinom` documentation for the interpretation of the "variance" term.
    # Here we reparameterize with the coefficient of variation (=sd/mean)
    a    = mu / (mu * prms$cv -1) 
    a.ww = mu.ww / (mu.ww * prms$cv.ww -1)
    a.ha = mu.ha / (mu.ha * prms$cv -1)
    
    # Avoid NaNs:
    tiny = 1e-6
    a[a <= 0]       <- tiny
    a.ww[a.ww <= 0] <- tiny
    a.ha[a.ha <= 0] <- tiny
    
    df$report.obs   = rnbinom(n=nrow(df), mu = df$report, size = a)
    df$WWreport.obs = rnbinom(n=nrow(df), mu = df$WWreport, size = a.ww)
    df$hosp.admission.obs = rnbinom(n=nrow(df), 
                                    mu = df$hosp.admission, size = a.ha)
    return(df)
}

