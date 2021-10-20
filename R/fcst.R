

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
                 dat = NULL, 
                 ci = 0.95, 
                 n.cores = 1) {
    
    prm = fitobj$prm
    h = prm$horizon
    prm$horizon <- h + horizon.fcst
    
    ss = wem::simul_from_post(post.abc = fitobj$post.abc, 
                              prm      = prm, 
                              ci       = ci, 
                              n.cores  = n.cores)
    
    if(!is.null(dat)) {
        d0 = min(dat$obs$date)
        ss = ss %>% 
            mutate(date = d0 + time)
    }
    
    res = list(
        sim.post = ss, 
        post.abc = fitobj$post.abc,
        prm = prm, 
        prm.abc = fitobj$prm.abc,
        ci = ci)
    
    if(!is.null(dat)) res = c(res, list(date.first.obs = d0))
    
    return(res)
}


