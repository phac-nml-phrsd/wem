


#' Return current date and time with a friendly file name format. Omits year.
#' @export
#' 
now_string <- function() {
    x  = as.character(lubridate::now())
    md = stringr::str_extract(x, '-\\d{2}-\\d{2}')
    t  = stringr::str_extract(x, '\\s\\d{2}:\\d{2}')
    t  = stringr::str_replace_all(t, ':','')
    t  = stringr::str_replace_all(t, ' ','-')
    return(paste0(md,t))
}



#' FITTING - Load all baseline parameters
load_prm_fit <- function(path = 'parameters/',weight) {
    # determine fitting source based on weight
    cl = weight[['cl']]
    ww = weight[['ww']]
    hosp = weight[['hosp.adm']]
    
    ww.only = (ww > 0 & cl==0 & hosp==0)
    
    prm.simul = read_prm_list(paste0(path,'prm-simul.csv'))
    prm.epi   = read_prm_list(paste0(path,'prm-epi.csv'))
    prm.ww    = read_prm_list(paste0(path,'prm-ww.csv'))
    
    if(ww.only){
        prm.intrv = read_prm_list(paste0(path,'prm-interv-ww.csv'))
    }
    else{
        prm.intrv = read_prm_list(paste0(path,'prm-interv.csv'))
    }
    
    prm.surv  = read_prm_list(paste0(path,'prm-surv.csv'))
    prm.tmp   = c(prm.simul, prm.epi, prm.ww, 
                  prm.intrv, prm.surv)
    
    prm = digest_prm(prm.tmp)
    return(prm)
}
#' find_col(ts, comp) finds columns for a specific compartment
#' input comp can be E,I,A,H and V and input ts is ODE time-series results
#' We expect the variable name for compartments to be either
#' - "comp"  when only there is one compartment: for example "H" 
#' - "compn" when there are several comaprtments: for example "Hn" and n=1,2,..
#'
#' @param ts 
#' @param col 
#'
#' @return column's number of compartment col
#' @export
#'
#' @examples
#' find_col(ts,"E")
find_col <- function(ts,col)
{
    singlecol = (names(ts) == paste0(col))
    multicol = grepl(paste0("^",col,"\\d+"), names(ts))
    
    if(any(singlecol) & any(multicol)){
        message('Cannot handle variable names with default compartment names')
        stop()
    }
    
    res = which(singlecol | multicol)
    return(res)
}





#' Calculate total numbers per epidemiological states 
#' by summing up values in all compartmental series in input variable comp.
#' input comp can be E,I,A, and V compartments
calc.all <- function(ts,comp)
{
    col =  find_col(ts,comp)
    if(length(col)>1)  res = rowSums(ts[,col])
    if(length(col)==1) res = ts[,col]
    return(res)
}



#' @title Piecewise linear function
#' @description For a given \code{x}, returns the value on a piecewise 
#' linear function as defined by 
#' vector \code{b} of its break points and associated 
#' values \code{v}.
#' It is typically used for relative factor.for time-dependent 
#' parameters. For example: if n = 2, broken line returns below v values:
#' t < t1: returns v1
#' t > t2: returns v2
#' t1 < t < t2: returns the line from (t1,v1) to (t2,v2)
#'
broken_line <- function(x, b, v) {
    res = NA
    n = length(b)
    
    if(length(v)!=n){
        message('Error in `broken_line()`: vectors must have the save size!')
        message(paste('length vector b =',length(b)))
        message(paste(b, collapse=' ; '))
        message(paste('length vector v =',length(v)))
        message(paste(v, collapse=' ; '))
        stop()
    }
    
    i = which(x <= b)[1]
    if(i==1 | is.na(i)){
        if(x <= b[1]) res = v[1]
        if(x > b[n]) res = v[n]
    }
    else{
        t1   = b[i-1]
        t2   = b[i]
        tmp1 = (x-t2)/(t1-t2)
        tmp2 = (x-t1)/(t2-t1)
        res  = tmp1 * v[i-1] + tmp2 * v[i]
    }
    return(res)
}




#' Step line function.
#' 
#' For all x < bb[i], y = vv[i]
#'
#' @param x Numeric. Value at which the step line function is evaluated.
#' @param b Numeric vector. Breakpoints when the step occurs.
#' @param v Numeric vector. Value of the step.
#'
#' @return The y-value where the step function is evaluated
#' 
step_line <- function(x, b, v) {
    res = NA
    n = length(b)
    stopifnot(length(v)==n)
    
    i = which(x <= b)[1]
    if(i==1 | is.na(i)){
        if(x <= b[1]) res = v[1]
        if(x > b[n]) res = v[n]
    }
    else{
        t1   = b[i-1]
        t2   = b[i]
        res  = v[i]
    }
    return(res)
}



#' Find the names that define an overwrite of vector elements
#' For example, "v_2" means replace the 2nd element of vector "v"
#' by the value assigned to prm[["v_2"]].
#' @param prm List of paramters.
#' 
overwrite_vectors <- function(prm, verbose = FALSE) {
    prm2 = prm
    idx = which(grepl('\\w+_\\d+', names(prm2)))
    
    if(length(idx) > 0){
        msg = paste('Overwriting vector elements:', names(prm2)[idx],'\n')
        if(verbose) message(msg)
        
        a = names(prm2)[idx]
        for(i in seq_along(idx)){
            ii = idx[i]
            tmp = unlist(str_split(a[i],'_'))
            prm.name = tmp[1]
            prm.pos  = as.integer(tmp[2])
            
            # Check names are consistent
            if(!prm.name %in% names(prm2)){
              msg = paste0('ERROR: User defined the prior parameter name `',prm.name,
                           '` but it is not found in the parameters defining the model')
              stop(msg)
            }
            
            # Replacing value
            prm2[[prm.name]][prm.pos] <- prm[[ii]]        
        }
    }
    return(prm2)
}


#' Forecast based on a fitted object.
#'
#' @param path.fitted.object 
#' @param df.ww.cl 
#' @param ci 
#' @param time.last.obs 
#' @param time.horizon 
#' @param n.cores 
#'
#' @return
#' @export
#'
fcst_from_post <- function(path.fitted.object, 
                           df.ww.cl, 
                           ci = 0.95,
                           time.last.obs = NULL,
                           time.horizon = NULL,
                           n.cores = 4) {
    
    
    # Last time of observation:
    if(is.null(time.last.obs)){
        d0 = min(df.ww.cl$date)
        last.time = as.numeric(last.date - d0)   
    }
    else{
        last.time = time.last.obs
    }
    
    
    message(paste0(' --- Forecasting\nLast time = ', 
                   last.time, ' (', d0+last.time ,')'))
    
    # Load baseline parameters
    # prm = load_prm(paste0('parameters/prm-',loc,'/'))
    prm = load_prm_fit(path = paste0('parameters/prm-',loc,'/'),
                       weight = prm.abc$weight)
    
    # Remove any inference done beyond last observed point:
    idx.future = which(prm$b > last.time)
    if(length(idx.future) > 0){
        prm$b <-prm$b[-idx.future]
        prm$v <-prm$v[-idx.future]
        
        # Remove associated columns (:fitted parameters) in `post.abc`:
        # TODO: write a function for that (used for `vv` too)
        nm  = names(post.abc)
        foo = paste('v',idx.future,sep='_')
        idx.rm = which(nm %in% foo)
        if(length(idx.rm)>0) post.abc = post.abc[, -idx.rm]
    }
    
    idx.future = which(prm$bb > last.time)
    if(length(idx.future) > 0){
        prm$bb <-prm$bb[-idx.future]
        prm$vv <-prm$vv[-idx.future]
        
        # Remove associated columns (:fitted parameters) in `post.abc`:
        nm  = names(post.abc)
        foo = paste('vv',idx.future,sep='_')
        idx.rm = which(nm %in% foo)
        if(length(idx.rm)>0) post.abc = post.abc[, -idx.rm]
    }
    
    if(!is.null(time.horizon)) prm$horizon <- time.horizon
    
    # Run simulations
    ss = simul_from_post(post.abc, prm, ci, n.cores)
    
    sim.post = ss %>% 
        mutate(date = d0 + time,
               wwloc = loc)
    
    return(list(
        sim.post = sim.post, 
        post.abc = post.abc,
        ci = ci,
        prm = prm, 
        prm.abc = prm.abc,
        last.date = last.date,
        date.first.obs = d0))
}


