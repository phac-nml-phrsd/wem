


#' @title Define parameters for Appoximate Bayesian Computation (ABC) fitting method. 
#' 
#' @description Determine the required parameters (iteration, acceptance and weight)
#' for the ABC fitting process. 
#' The weights specify the importance of different data sources utilized simultaneously 
#'  in the fitting process.  There are three data sources 
#' for simultaneous fitting: reported cases, viral concentration in wastewater 
#' and hospitalization (admission or occupancy).  
#'
#' @param iter Numeric. Number of prior iterations for ABC fitting. 
#' @param accept Numeric. Acceptance ratio (so the number of posterior 
#' samples is \code{iter * accept}). 
#' @param case.weight Numeric, float. Relative weight for clinical cases 
#' @param ww.weight Numeric, float. Relative weigth for viral concentration in wastewater. 
#' @param hosp.weight Numeric, float. Relative weight for hospitalization.
#' @param hosp.type String. Type of hospitalization data: 
#' \code{"hosp.adm"} for hospital admissions, \code{"hosp.occ"} for hospital occupancy and
#' \code{NULL} for no hospital data.
#'
#' @return Nested list of ABC parameters.
#' @export
#' 
#' @seealso \code{define_fit_priors()}, \code{fit()}
#'
#' @examples 
#' 
#' prm.abc = define_abc_prms(
#' iter        = 1e4,
#' accept      = 1e-2,
#' case.weight = 1.0,
#' ww.weight   = 2.5,
#' hosp.weight = 1.0,
#' hosp.type   = 'hosp.adm')
#'
#'
define_abc_prms <- function(iter,
                            accept,
                            case.weight,
                            ww.weight,
                            hosp.weight,
                            hosp.type) {
    
    
    if(!is.null(hosp.type)){
        if(!hosp.type %in% c('hosp.adm','hosp.occ')){
            stop(paste0('Hospital type `hosp.type=',hosp.type,'` unknown in function call `define_abc_prms()`. Aborting.'))
        }
        
        if(hosp.type == 'hosp.adm'){
            w = c(cl = case.weight,
                  ww = ww.weight,
                  hosp.adm = hosp.weight)
        }
        if(hosp.type == 'hosp.occ'){
            w = c(cl = case.weight,
                  ww = ww.weight,
                  hosp.occ = hosp.weight)
        }
    }
    
    if(is.null(hosp.type)){
        # hosp.type is NULL and hosp.weight must be zero
        if(hosp.weight != 0){
            msg = 'hosp.type is NULL. For fitting process, hosp.weight must be 0.0'
            stop(msg)
        }
        w = c(cl = case.weight,
              ww = ww.weight,
              hosp.adm = hosp.weight,
              hosp.occ = hosp.weight)
    }
    
    
    
    
    res = list(n      = iter,
               accept = accept,
               weight = w)
    return(res)
}



#' @title Randomly sample from priors for fitting process 
#' 
#' @description  This function uses the dataframe output of \code{define_priors()} that defines the prior distributions for all parameter that must be fitted.
#'
#' @param prior Dataframe. Priors' parameters and distributions. Output of \code{define_fit_priors()}.
#' @param prm.abc List. Varibales required for ABC fitting method created by function define_abc_prms(). 
#' @param all.positive Logical. Default is \code{TRUE}, which means sampled values for priors will be truncated to be positive.
#'
#' @return Long dataframe. All the prior sampled values for the parameters to be fitted.

#' @export
#'
#' 
sample_priors <- function(prior, prm.abc, all.positive = TRUE) {
    
    # Draw the prior from their specified distribution:
    n.abc = prm.abc[['n']]
    x = list()
    for(i in 1:nrow(prior) ){
        prms.i = strsplit(prior$prms[i], split = ';')[[1]] %>%
            as.numeric() %>% as.list()
        dist.i = prior$distrib[i]
        x[[i]] = do.call(dist.i, c(n=n.abc, prms.i))
    }
    res = do.call('cbind.data.frame', x)
    names(res) <- prior$name
    
    # TODO: Do this name dependent?
    if(all.positive) res[res<0] <- 1e-9
    
    return(res)
}


#' @title  Prior distributions for fitted model parameters.
#'
#' @description Create a dataframe with the definition of the prior distribution 
#' for all parameters that will be fitted to data. 
#' The definition are written in a CSV file.
#'
#' @param path String. Path to the CSV file. 
#' The CSV file must have a specific format: 
#' \itemize{
#' \item column \code{name}: the variable name must match the name of the 
#' model parameters that will be fitted. 
#' For example, the basic reproduction number is \code{R0} (and not, say, \code{Rzero}).
#' \item column \code{distrib}: the name of the distribution type for the prior variable. Use the standard nomenclature for probability distribution in R (e.g., \code{rnorm()}, \code{rexp()}, ...)
#' \item column \code{prms}: the parameters for the distribution type for the prior variable. Use the standard nomenclature for probability distribution in R (e.g., \code{mean} and \code{sd} for \code{rnorm()}, etc.) separated by a semi-colon \code{;}. 
#' }
#' 
#' @return Dataframe with three columns: name, distrib,prms.  
#' @export
#' 
#' @examples 
#' \dontrun{
#' An example for the CSV file is:
#' 
#' name,        distrib, prms
#' R0,          rnorm,   3.0; 0.2
#' transm.v_3,  rexp,    0.10 
#' 
#' which translates into:
#'  - The parameter R0 has a normal distribution 
#'    with mean 3.0 and std dev 0.2 as a prior distribution. 
#'  - The 3rd element of the vector that represents 
#'    the time-dependent transmission rate `transm.v` 
#'    (beta in the mathematical documentation) has an 
#'    exponential prior distribution with mean 1/0.10.
#' }
#' 
#' 
define_fit_priors <- function(path) {
    
    # Retrieve the distribution definitions
    res = read.csv(path, strip.white = TRUE)
    return(res)
}



#' @title Estimate statistical values after fiiting
#'
#' @param post.abc Results of abc fitting.
#' @param ci.tight Credible interval 
#' @param ci.broad Credible interval
#' @param df.true 
#'
#' @return
#' @export
#'
post_ss_abc <- function(post.abc, ci.tight=0.50, ci.broad=0.95,
                        df.true = NULL) {
    
    df = post.abc %>%
        pivot_longer(cols = 1:ncol(post.abc)) %>%
        group_by(name) %>%
        summarise(m = mean(value),
                  qvhi = quantile(value, probs = 0.5 + ci.broad/2), 
                  qhi  = quantile(value, probs = 0.5 + ci.tight/2), 
                  qlo  = quantile(value, probs = 0.5 - ci.tight/2), 
                  qvlo = quantile(value, probs = 0.5 - ci.broad/2)
        ) 
    
    col.post = 'steelblue3'
    sz.seg = 4
    alpha.seg = 0.6
    
    g = df %>% 
        ggplot() + 
        geom_segment(aes(x=name,xend=name, y=qvlo, yend=qvhi),size=sz.seg,
                     col=col.post, alpha=alpha.seg*0.8) +
        geom_segment(aes(x=name,xend=name, y=qlo, yend=qhi),size=sz.seg, 
                     col=col.post, alpha=alpha.seg) +
        geom_point(aes(x=name, y=m), shape=21, size=3,stroke=1, 
                   fill='white', col=col.post) + 
        facet_wrap(~name, scales = 'free') + 
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              strip.background = element_rect(fill=col.post, linetype = 0),
              strip.text = element_text(face='bold', colour = 'white'),
              axis.ticks = element_blank(),
              axis.text.x = element_blank())+
        xlab('') + ylab('') + 
        ggtitle('Posterior Summary Stats',
                paste0('Mean, ',ci.tight*100,'% and ',ci.broad*100,'% CrI'))
    
    if(!is.null(df.true)){
        g = g + geom_hline(data = df.true, 
                           aes(yintercept = true.value),
                           linetype = 'dashed')
    }
    
    return(list(plot = g, df = df))
    
}



#' @title Calculate error for ABC fitting
#'
#' @param x Observed data
#' @param y Simulated data
#'
#' @return Numerical value for error's fit
#'
err_fct <- function(x,y) {
    res = sum( (x-y)^2 , na.rm = TRUE)
    return(sqrt(res))
}

#' @title Perform fitting for one iteration
#'
#' @param i Numerical. The number of iterated simulation for fitting.  
#' @param prm.abc List. Parameters for ABC fitting.
#' @param prm List. All parameters required for model (\code{wem-prm.csv})
#' @param lhs Dataframe. Sampled for priors.
#' @param obs Datafrmae. Paired clinical-ww data created by \code{build_data()}
#' @param hosp.var String. Type of hospitalization data. Output of function \code{build_data()}.   
#' @param case.var String. Type of date which cases are reported. Output of function \code{build_data()}.   
#'
#' @return
#'
fit_abc_unit <- function(i, 
                         prm.abc, 
                         prm, 
                         lhs, 
                         obs,
                         hosp.var,
                         case.var) {
    
    # Update parameter values and simulate:
    np = ncol(lhs)
    
    # Check that the name of scalar parameter exist:
    check_prior_name_scalar(lhs, prm)
    
    # Add vector elements and replace scalar elements:
    for(k in 1:np){
        prm[[ names(lhs)[k] ]] = lhs[i,k]
    }
    
    # Simulate with updated parameters:
    sim.i = simul(prm)$ts 
    
    # Normalize by max value of observations:
    mx.obs.ww = max(obs$ww.obs, na.rm = TRUE)
    mx.obs.cl = max(obs$clin.obs, na.rm = TRUE)
    
    if(case.var == 'report'){
        sim.norm = sim.i %>% 
            mutate(time = round(time, 1)) %>%
            mutate(norm.ww.sim = WWreport / mx.obs.ww,
                   norm.cl.sim = report / mx.obs.cl) %>% 
            select(time, starts_with('norm'))
    }
    
    if(case.var == 'episode'){
        sim.norm = sim.i %>% 
            mutate(time = round(time, 1)) %>%
            mutate(norm.ww.sim = WWreport / mx.obs.ww,
                   norm.cl.sim = report.episode / mx.obs.cl) %>% 
            select(time, starts_with('norm'))
    }
    
    obs.norm = obs %>% 
        mutate(norm.ww = obs$ww.obs / mx.obs.ww,
               norm.cl = obs$clin.obs /mx.obs.cl)
    
    # If fitting hospitalization data:
    if(!is.null(hosp.var)){
        if(hosp.var == 'hosp.adm'){
            mx.obs.hosp.adm = max(obs$hosp.obs, na.rm = TRUE)
            sim.norm$norm.hosp.adm.sim <- sim.i$hosp.admission / mx.obs.hosp.adm
            obs.norm$norm.hosp.adm     <- obs$hosp.obs / mx.obs.hosp.adm
        }
        if(hosp.var == 'hosp.occ'){
            mx.obs.hosp.occ = max(obs$hosp.obs, na.rm = TRUE)
            sim.norm$norm.hosp.occ.sim <- sim.i$Hall / mx.obs.hosp.occ
            obs.norm$norm.hosp.occ     <- obs$hosp.obs / mx.obs.hosp.occ
        }
    }
    
    # Match the times
    # (observations may be irregular 
    # and missing, unlike for simulations)
    dj = left_join(obs.norm, sim.norm, by='time') %>% 
        filter(!is.na(norm.ww) | !is.na(norm.cl))
    
    # Calculate error between simulations and observations:
    err <- vector()
    err['ww'] = err_fct(dj$norm.ww, dj$norm.ww.sim)
    err['cl'] = err_fct(dj$norm.cl, dj$norm.cl.sim)
    
    if(!is.null(hosp.var)){
        if(hosp.var == 'hosp.adm')
            err['hosp.adm'] = err_fct(dj$norm.hosp.adm, dj$norm.hosp.adm.sim)
        if(hosp.var == 'hosp.occ')
            err['hosp.occ'] = err_fct(dj$norm.hosp.occ, dj$norm.hosp.occ.sim)
    }
    
    # Normalized weights
    weightfit = prm.abc$weight
    ws        = sum(prm.abc$weight)
    weightfit = weightfit / ws
    
    # Make sure the sums are between 
    # paired elements of the same type (e.g., ww with ww):
    err.total = 0
    for(q in names(err)){
        err.total = err.total + weightfit[q] * err[q]
    }
    
    abcerr =  c(
        err.total     = unname(err.total),
        err.ww        = unname(weightfit['ww'] * err['ww']), 
        err.clin      = unname(weightfit['cl'] * err['cl']),
        err.hosp.adm  = unname(weightfit['hosp.adm'] * err['hosp.adm']),
        err.hosp.occ  = unname(weightfit['hosp.occ'] * err['hosp.occ']))
    
    res = list(errors = abcerr, 
               ts = sim.i)
    return(res)
}





#' Fitting ABC function (internal use)
#'
#' @param obs Dataframe. Paired clinical and ww data. Output of function \code{build_data()}.
#' @param priors Dataframe. Sampled priors for fitting. Output of function \code{sample_priors()}
#' @param hosp.var String. Output of function \code{build_data()}.
#' @param case.var String. Output of function \code{build_data()}.
#' @param prm.abc List. Parameters for ABC fitting method. Output of function \code{define_abc_prms} 
#' @param prm List. All parameters model requires. It has to be in the \code{wem-prm.csv} format.    
#' @param n.cores Numerical. Number of cores for fitting compoutation. 
#'
#' @return A list of 
#' 1) a dataframe of variables post fitting, 
#' 2) a dataframe of priors used for fitting,
#' 3) a dataframe of fitting errors,
#' 4) a dataframe of the final states of the system for the posterior trajectories only.
#'
#' 
fit_abc <- function(obs,
                    priors,
                    hosp.var,
                    case.var,
                    prm.abc, 
                    prm, 
                    n.cores=2) {
    
    # Sampled Priors
    lhs = priors
    
    # Run simulation for every parameter combination
    sfInit(parallel = n.cores>1, cpus = n.cores)
    sfExportAll()
    suppressMessages({
        sfLibrary(deSolve)
        sfLibrary(stringr)
        sfLibrary(dplyr)
    })
    
    tmp = sfLapply(x   = 1:prm.abc$n,
                   fun = fit_abc_unit,
                   prm.abc = prm.abc,
                   prm = prm,
                   lhs = lhs,
                   obs = obs,
                   hosp.var = hosp.var,
                   case.var = case.var)
    sfStop()
    
    # Sort the fitting errors after retrieving them:
    df.err = lapply(tmp, '[[', 'errors') %>%
        bind_rows() %>%
        mutate(i = 1:prm.abc$n) %>%
        as.data.frame() %>%
        arrange(err.total)
    
    # Keep only the smallest errors:
    n.accept = round(prm.abc$accept * prm.abc$n)
    idx.keep = df.err$i[1:n.accept]
    df.post  = lhs[idx.keep,]
    
    # Record the final states:
    fs = lapply(tmp, '[[', 'ts') %>%
        bind_rows() %>%
        filter(time == max(time)) %>%
        mutate(i = 1:prm.abc$n) %>%
        filter(i %in% idx.keep)

    return(list(df.post = df.post, 
                priors  = lhs,
                err     = df.err,
                finalstates = fs))
}


#' @title Helper function
#'
#' @param data 
#' @param last.date 
#'
unpack_data <- function(data, last.date) {
    obs        = data[['obs']]
    obs.long   = data[['obs.long']]
    hosp.var   = data[['hosp.var']]
    case.var   = data[['case.var']]
    
    if(!is.null(last.date)){
        obs      = filter(obs, date <= last.date)
        obs.long = filter(obs.long, date <= last.date)
    }
    if(is.null(last.date)){
        last.date = max(obs$date)
    }
    return( list(
        obs = obs,
        obs.long = obs.long,
        last.date = last.date,
        hosp.var = hosp.var,
        case.var = case.var
    ))
}




#' @title Fitting Data to Model
#'
#' @param data List. Output of function build_data()
#' @param prm.abc List. Variables for ABC fitting. Output of function define_abc_prms()
#' @param df.priors Dataframe. Parameters (\code{'name'}), their distribution (\code{'distrib'})
#' and range of values (\code{'prms'}) which to be fitted to the model. Output of function define_fit_priors(). 
#' @param prm List. Initial parameters for fitting plot. Output of function load_prm().
#' @param n.cores Numeric. number of cores used for fitting computation.
#' @param last.date Date. Last date to stop fitting (e.g., ymd('2021-02-01') or NULL). NULL = latest date available.
#' @param do.plot Logical. Plot the fitting results and initial parameters  
#' @param save.rdata Logical. Saving the fitting objects in a .rds file.
#'
#' @return A list containing the fitted object and additional information.
#' @export
#'
fit <- function(data,
                prm.abc,
                df.priors,
                prm,
                n.cores,
                last.date = NULL, 
                save.rdata = FALSE) {
    
    d = unpack_data(data, last.date)
    obs        = d[['obs']]
    obs.long   = d[['obs.long']]
    hosp.var   = d[['hosp.var']]
    case.var   = d[['case.var']]
    last.date  = d[['last.date']]
   
    # --- Draw priors
    samp.priors = sample_priors(df.priors,prm.abc)
    
    # --- Run ABC
    elapsed.time = system.time({
        fit = fit_abc(obs = obs, 
                      priors = samp.priors,
                      hosp.var = hosp.var,
                      case.var = case.var,
                      prm.abc = prm.abc,
                      prm = prm, 
                      n.cores = n.cores)
    })
    print(elapsed.time)
    
    # Retrieve the priors and 
    # posterior distributions:
    post.abc = fit$df.post 
    df.prior = fit$priors %>%
        pivot_longer(cols = 1:length(fit$priors))
    post.ss  = post_ss_abc(post.abc)
    
    time.completed = now_string()
    
    res = list(
        prm = prm,
        prm.abc = prm.abc,
        fit = fit,
        hosp.var = hosp.var,
        case.var = case.var,
        post.abc = post.abc,
        post.ss = post.ss,
        df.prior = df.prior,
        obs = obs, 
        obs.long = obs.long,
        last.date = last.date,
        elapsed.time = elapsed.time,
        time.completed = time.completed
    )
    
    if(save.rdata){
        save(list = names(res), 
             file =  paste0('fitted', time.completed, '.RData'))
    }
    return(res)
}


#' @title Fitting model to recent data only, using existing fit on past data. 
#'
#' @param data List. Output of function build_data()
#' @param prm.abc List. Variables for ABC fitting. Output of function define_abc_prms()
#' @param fitobj.past List as returned by the function \code{fit()}. Fit object for past data. 
#' @param df.priors Dataframe. Parameters (\code{'name'}), their distribution (\code{'distrib'})
#' and range of values (\code{'prms'}) which to be fitted to the model. Output of function define_fit_priors(). 
#' @param prm List. Initial parameters for fitting plot. Output of function load_prm().
#' @param n.cores Numeric. number of cores used for fitting computation.
#' @param last.date Date. Last date to stop fitting (e.g., ymd('2021-02-01') or NULL). NULL = latest date available.
#' @param do.plot Logical. Plot the fitting results and initial parameters  
#' @param save.rdata Logical. Saving the fitting objects in a .rds file.
#'
#' @return A list containing the fitted object and additional information.
#' @export
#'
fit_recent <- function(data,
                       prm.abc,
                       fitobj.past,
                       df.priors,
                       prm,
                       n.cores,
                       last.date = NULL, 
                       save.rdata = FALSE) {
    
    d = unpack_data(data, last.date)
    obs        = d[['obs']]
    obs.long   = d[['obs.long']]
    hosp.var   = d[['hosp.var']]
    case.var   = d[['case.var']]
    last.date  = d[['last.date']]
    
    # --- Draw priors
    
    # Retrieve samples from posteriors 
    # fitted on past data:
    post.past = fitobj.past$post.abc
    n.post.past = nrow(post.past)
    
    # make sure number of priors for recent 
    # is multiple of past posteriors
    q = round(prm.abc$n / n.post.past)
    prm.abc$n <- n.post.past * q
    
    # Sample variables fitted to recent data
    samples.recent = sample_priors(prior = df.priors, 
                                   prm.abc = prm.abc)
    
    # Stitch samples from past posteriors
    #  and priors for recent data.
    # Note: the past posteriors are repeated several
    # times, but the recent priors are unique. This is 
    # to ensure consistency between past and recent.
    
    # First, check that the variable names for 
    # "recent" and "past" or not the same:
    if(any(names(samples.recent) %in% names(post.past))){
        stop('In `fit_recent()`, the recent variables cannot be the same as the variables used to fit past data. Correct variable name in `df.priors`. Aborting.')
    }
    
    # Stitch
    samples.full = cbind(post.past, samples.recent)
    
    
    # --- Run ABC
    elapsed.time = system.time({
        fit = fit_abc(obs = obs, 
                      priors = samples.full,
                      hosp.var = hosp.var,
                      case.var = case.var,
                      prm.abc = prm.abc,
                      prm = prm, 
                      n.cores = n.cores)
    })
    print(elapsed.time)
    
    # Retrieve the priors and 
    # posterior distributions:
    post.abc = fit$df.post 
    df.prior = fit$priors %>%
        pivot_longer(cols = 1:length(fit$priors))
    post.ss  = post_ss_abc(post.abc)
    
    time.completed = now_string()
    
    res = list(
        prm = prm,
        prm.abc = prm.abc,
        fit = fit,
        hosp.var = hosp.var,
        case.var = case.var,
        post.abc = post.abc,
        post.ss = post.ss,
        df.prior = df.prior,
        obs = obs, 
        obs.long = obs.long,
        last.date = last.date,
        elapsed.time = elapsed.time,
        time.completed = time.completed
    )
    
    if(save.rdata){
        save(list = names(res), 
             file =  paste0('fitted-recent', time.completed, '.RData'))
    }
    
    return(res)
}

# Thu Dec 30 08:49:08 2021 ------------------------------

find_timedep_prms <- function(prm, t.or.v) {
    # t.or.v = 't'
    idx = which(grepl(paste0('\\.',t.or.v,'$'), names(prm)))
    # Ignore the ones with NULL values:
    ii = which(prm[idx] != 'NULL')
    idx = idx[ii]
    return(idx)
}

rebase_time_prms <- function(prm, t.pivot) {

    # Identify the time-depedent parameters:    
    idx.t = find_timedep_prms(prm, 't')
    prmt = prm[idx.t]
    
    # Compare all times to the pivot time:
    chk = sapply(prmt, FUN = function(x){x>t.pivot})
    # Is at least one greater than pivot time:
    chk2 = sapply(chk, any)        
    
    y = list()
    
    for(i in seq_along(chk2)){
        
        # For those that have no time definition after pivot, 
        # simply convert to a constant set at the last value
        if(!chk2[i]){
            name.t = names(prmt)[i]
            name.v = str_replace(name.t, 't$','v') # find the associated values
            tmp = prm[[name.v]]
            y[[name.t]] <- 0   # reset time at 0
            y[[name.v]] <- tmp[length(tmp)]
        }
        
        # For those that have a time definition after the 
        # pivot time, translate time and values:
        if(chk2[i]){
            tmp <- prmt[[i]] - t.pivot
            tmpt = tmp[tmp>0]
            name.t = names(prmt)[i]
            name.v = str_replace(name.t, 't$','v') # find the associated values
            tmpv = prm[[name.v]]
            # Keep only the values associated 
            # with the times beyond pivot time:
            nt = length(tmpt)
            nv = length(tmpv)
            tmpv = tmpv[c( (nv-nt+1):nv )]
            y[[name.t]] <- tmpt
            y[[name.v]] <- tmpv
        }
    }
    # overwrite the rebased elements
    res = prm
    res[names(y)] <- y
    return(res)
}

rebase_at_pivot <- function(data, date.pivot, prm) {
    
    obs      = filter(data[['obs']],      date >= date.pivot)
    obs.long = filter(data[['obs.long']], date >= date.pivot)

    date.origin = min(data[['obs']]$date)    
    t.pivot = as.numeric(date.pivot - date.origin)
    
    prm.rebased = rebase_time_prms(prm, t.pivot)
    
    return( list(
        obs      = obs,
        obs.long = obs.long,
        prm      = prm.rebased,
        hosp.var = data[['hosp.var']],  # NEED TO DO THAT??
        case.var = data[['case.var']]   # NEED TO DO THAT??
    ))
}


#' @title Fitting model to recent data only, using existing fit on past data. 
#' @return A list containing the fitted object and additional information.
#' @export
#'
fit_recent_NEW <- function(
    data,
    prm.abc,
    fitobj.past, # NEED THIS FOR FINAL STATES
    df.priors,
    prm,
    n.cores,
    date.pivot,
    save.rdata = FALSE) {
    
    # retrieve new initial values 
    # that start at the `pivot.date`:
    initvar = fitobj.past$fit$finalstates
    
    x = rebase_at_pivot(data, date.pivot, prm)
    
    # d = unpack_data(data, last.date) # DELETE WHEN SURE
    obs        = x[['obs']]
    obs.long   = x[['obs.long']]
    hosp.var   = x[['hosp.var']]
    case.var   = x[['case.var']]
    prm.rebased = x[['prm']]
        
    # --- Draw priors
    
    # Thu Dec 30 09:45:53 2021 ------------------------------
    # STOPPED HERE
    # `prm.rebased` should be used with initvar 
    # in the standard `fit()` (??)
    
# ---- MAYBE USELESS NOW....
    # Retrieve samples from posteriors 
    # fitted on past data:
    post.past = fitobj.past$post.abc
    n.post.past = nrow(post.past)
    
    # make sure number of priors for recent 
    # is multiple of past posteriors
    q = round(prm.abc$n / n.post.past)
    prm.abc$n <- n.post.past * q
# -----------  (useless?)
    
    
    # Sample variables fitted to recent data
    priors.recent = sample_priors(prior = df.priors, 
                                   prm.abc = prm.abc)
    
    # Stitch samples from past posteriors
    #  and priors for recent data.
    # Note: the past posteriors are repeated several
    # times, but the recent priors are unique. This is 
    # to ensure consistency between past and recent.
    
    # First, check that the variable names for 
    # "recent" and "past" or not the same:
    if(any(names(priors.recent) %in% names(post.past))){
        stop('In `fit_recent()`, the recent variables cannot be the same as the variables used to fit past data. Correct variable name in `df.priors`. Aborting.')
    }
    
    # Stitch
    # samples.full = cbind(post.past, samples.recent)
    
    
    # Wed Dec 29 15:34:47 2021 ------------------------------
    #
    # problem to solve: 
    # must rebase time 0 to `last.date`
    # and also make sense of the name of 
    # the prior variable.
    # For example  `transm.v_4` must be "translated" 
    # into something like `transm.v_1` once time
    # is rebased at `last.date` . 
    # 
    
    # --- Run ABC
    elapsed.time = system.time({
        fit = fit_abc(obs = obs, 
                      priors = priors.recent,
                      hosp.var = hosp.var,
                      case.var = case.var,
                      prm.abc = prm.abc,
                      prm = prm, 
                      n.cores = n.cores)
    })
    print(elapsed.time)
    
    # Retrieve the priors and 
    # posterior distributions:
    post.abc = fit$df.post 
    df.prior = fit$priors %>%
        pivot_longer(cols = 1:length(fit$priors))
    post.ss  = post_ss_abc(post.abc)
    
    time.completed = now_string()
    
    res = list(
        prm = prm,
        prm.abc = prm.abc,
        fit = fit,
        hosp.var = hosp.var,
        case.var = case.var,
        post.abc = post.abc,
        post.ss = post.ss,
        df.prior = df.prior,
        obs = obs, 
        obs.long = obs.long,
        last.date = last.date,
        elapsed.time = elapsed.time,
        time.completed = time.completed
    )
    
    if(save.rdata){
        save(list = names(res), 
             file =  paste0('fitted-recent', time.completed, '.RData'))
    }
    
    return(res)
}

