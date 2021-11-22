
#------ Calculate Rt 
#' @title Estimate Effective Reproduction Number.
#'
#' @description Calculate the effective reproduction number from a model 
#' fitted on observations. The fitted parameters are sampled from 
#' their posterior distribution to simulate epidemic curves. 
#' The effective reproduction number, Rt, is estimated from those simulations. 
#' 
#'
#' @param fitobj List. Results of the fitted model. Output of function \code{fit()}. 
#' @param ci Numerical. Confidence interval (default = 0.95).
#' @param n.cores Numerical. Number of computer CPU cores used for computation.
#' (default = 1)
#'
#' @return A dataframe of estimates for the mean Rt and its confidence interval bounds. 
#' @export
#'
estimate_Rt <- function(fitobj, ci = 0.95, n.cores = 1){
  
  print('--> Estimating effective reproduction number')   
  
  # ---- Rt from mechanistic model
  sim.Rt = calc_R_from_model(fitobj = fitobj, ci=ci, n.cores=n.cores)
  
  return(sim.Rt)
}


#' @title Estimate Effective Reproduction Number 
#'
#' @param fitobj List. Results of the fitted model. Output of function \code{fit()}. 
#' @param ci Numerical. Confidence Interval.
#' @param n.cores Numerical. Number of computer CPU cores used for computation.
#'
#' @return
#' @export
#'
calc_R_from_model <- function(fitobj,ci,n.cores){
  
  # UNPACK
  prm = fitobj$prm
  prm.abc  = fitobj$prm.abc
  breaks.val  = prm[['transm.v']]  
  breaks.time = prm[['transm.t']]
  
  
  # run simulation with posterior distribution
  runpost  = run_from_posteriors(fitobj = fitobj, 
                                 ci = ci,
                                 n.cores = n.cores)
  
  
  # -- Rt
  df.R = calc_Reff_from_post(runpost, breaks.time, breaks.val)
  
  df.R = df.R %>%
    mutate(width.ci = Reff.hi - Reff.lo)
  
  return(df.R)
}

#' @title Run Simulation from Porsterios
#'
#' @param fitobj List. Results of the fitted model. Output of function \code{fit()}. 
#' @param ci Numerical. Confidence Interval.
#' @param time.horizon Numerical. Maximum time for simulation.
#' @param n.cores Numerical. Number of computer CPU cores used for computation.
#'
#' @return
#' @export
#'
run_from_posteriors <- function(fitobj, 
                                ci=0.95,
                                time.horizon = NULL,
                                n.cores = 4) {

  
  # Load parameters from fitted object
  obs = fitobj$obs
  prm = fitobj$prm
  prm.abc = fitobj$prm.abc
  post.abc = fitobj$post.abc
  last.date = fitobj$last.date
  hosp.var = fitobj$hosp.var
  case.var = fitobj$case.var
  
  
  if(!is.null(time.horizon)) prm$horizon <- time.horizon
  
  # Run simulations
  ss = simul_from_post(post.abc, prm, hosp.var, case.var, ci=ci, n.cores=n.cores)
  d0 = min(obs$date)
  sim.post = ss %>% 
    mutate(date = d0 + time)
  
  
  return(list(
    sim.post = sim.post, 
    post.abc = post.abc,
    hosp.var = hosp.var,
    case.var = case.var,
    ci = ci,
    prm = prm, 
    prm.abc = prm.abc,
    last.date = last.date,
    date.first.obs = d0))
}

#' @title Calculate Rt from Posteriors 
#'
#' @param runpost Output of function \code{run_from_posteriors}
#' @param breaks.time Numerical. Time when transmission rate changes.
#' @param breaks.val Numerical. Value for the change of transmission rate at \code{breaks.time}
#'
#' @return
#' @export
#'
calc_Reff_from_post <- function(runpost,
                                breaks.time,
                                breaks.val) {
  
  # Unpack:
  sim.post = runpost$sim.post
  post.abc = runpost$post.abc
  pop.size = sim.post$S.m[1]
  d0       = runpost$date.first.obs
  ci       = runpost$ci
  
  # Stats on the transmission multiplicative factor: 
  m  = apply(post.abc, MARGIN=2, FUN=mean)
  lo = apply(post.abc, MARGIN=2, FUN=quantile, probs = 0.5 - ci/2)
  hi = apply(post.abc, MARGIN=2, FUN=quantile, probs = 0.5 + ci/2)
  
  tmax = max(sim.post$time)
  tvec = 0:tmax
  
  v.m  = m[grepl('^v_', names(m))]
  v.lo = lo[grepl('^v_', names(lo))]
  v.hi = hi[grepl('^v_', names(hi))]
  idx  = as.integer(stringr::str_extract( names(v.m), '\\d+'))
  
  bv.m = bv.lo = bv.hi = breaks.val
  
  bv.m[idx]  <- v.m
  bv.lo[idx] <- v.lo
  bv.hi[idx] <- v.hi
  
  y.m  = sapply(X=tvec, FUN = broken_line, b = breaks.time, v = bv.m)
  y.lo = sapply(X=tvec, FUN = broken_line, b = breaks.time, v = bv.lo)
  y.hi = sapply(X=tvec, FUN = broken_line, b = breaks.time, v = bv.hi)
  
  # Retrieve R0:
  R0.m  = m['R0']
  R0.lo = lo['R0']
  R0.hi = hi['R0']
  
  # Calculate effective reproduction number::
  res = data.frame(
    time    = sim.post$time,
    date    = d0 + sim.post$time,
    Reff.m  = R0.m  * y.m  * (sim.post$S.m  + sim.post$V.m)  / pop.size,
    Reff.lo = R0.lo * y.lo * (sim.post$S.lo + sim.post$V.lo) / pop.size,
    Reff.hi = R0.hi * y.hi * (sim.post$S.hi + sim.post$V.hi) / pop.size
  )
  return(res)
}

