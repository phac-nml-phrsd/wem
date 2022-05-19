#' @title Remove future inference in parameters.
#'
#' @description Remove elements beyond last date of observed data 
#' in the vector of time-dependent prms ('variable.t' and 'variable.v') or
#' in a column priors ('variable_t', 'variable_v')
#'
#' @param idx.t 
#' @param idx.future 
#' @param prm 
#' @param post.abc 
#'
#' @return 
#'
remove.fut.inference <- function(idx.t,
                                 idx.future,
                                 prm,post.abc){
    nm.t = names(prm[idx.t])
    nm.v = paste0(stringr::str_remove(nm.t,'.t$'),'.v')
    idx.v = which(names(prm)==nm.v)
    
    prm[[idx.t]] <- prm[[idx.t]][-idx.future] # remove time
    
    if(length(idx.v)>0){
        prm[[idx.v]] <- prm[[idx.v]][-idx.future] # remove value
        
        # Remove associated columns (:fitted parameters) in `post.abc`:
        nm.post = names(post.abc)
        foo.t = paste(nm.t,idx.future,sep='_') #future index in time
        foo.v = paste(nm.v,idx.future,sep='_') #future index in value
        idx.rm = which(nm.post %in% c(foo.t,foo.v))
        if(length(idx.rm)>0) post.abc = post.abc[, -idx.rm]
    }
    return(list(prm=prm,post.abc=post.abc))
}


#' @title Forecast from a fit object.
#'
#' @param fitobj 
#' @param horizon.fcst 
#' @param dat 
#' @param ci 
#' @param n.cores 
#'
#' @return
#' @export
#'
fcst <- function(fitobj, 
                 horizon.fcst,
                 dat, 
                 ci = 0.95, 
                 n.cores = 1) {
    
    last.date = fitobj$last.date
    
    #check if last.date of fitted object and data matches
    if(max(dat$obs$date)!= last.date){
        msg = 'last.date of fitted object does not match last observed data point in dat'
        stop(msg)
    } 
    
    hosp.var  = fitobj$hosp.var 
    case.var  = fitobj$case.var
    post.abc  = fitobj$post.abc
    prm       = fitobj$prm
    
    # Last time of observation and simulation horizon:
    d0 = min(dat$obs$date)
    last.time = as.numeric(last.date - d0)
    
    prm$horizon <- last.time + horizon.fcst
    
    # Remove any inference done beyond last observed point:
    idx.timedep = which(grepl('\\.t$', names(prm)))
    
    if(0){ # DEBUG
        names(prm)[idx.timedep]
        prm[idx.timedep]
    }
    
    for(idx.t in idx.timedep){
        #print(names(prm)[idx.t]) # DEBUG
        idx.future = NULL
        if(prm[[idx.t]][1] !='NULL')
            idx.future = which(prm[[idx.t]] > last.time)
        
        if(length(idx.future) > 0){
            rem = remove.fut.inference(idx.t,
                                       idx.future,
                                       prm,post.abc)
            prm      = rem$prm
            post.abc = rem$post.abc
        }
    }
    
    # Run simulation
    ss = simul_from_post(post.abc = post.abc, 
                         prm      = prm,
                         hosp.var = hosp.var,
                         case.var = case.var,
                         ci       = ci, 
                         n.cores  = n.cores)
    
    ss = ss %>% 
        mutate(date = d0 + time)
    
    res = list(
        sim.post = ss, 
        post.abc = post.abc,
        prm = prm, 
        prm.abc = fitobj$prm.abc,
        hosp.var = hosp.var,
        case.var = case.var,
        ci = ci,
        date.first.obs = d0,
        date.last.obs  = last.date)
    
    return(res)
}


