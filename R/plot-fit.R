###
###   FUNCTION TO GENERATE PLOTS ASSOCIATED WITH THE FITTING
###


#' @title Visualize Fitting Results and Statistics
#'
#' @param fitobj List. Fitted object. Output of function \code{fit()}.
#' @param ci Numerical. Confidence interval.
#'
#'
#' @return A list of \code{ggplot2} objects:
#' \itemize{
#' \item error:  Plot statistical errors from fitting different data sources. 
#' \item statistic.summary:  Plot statistics of the posteriors. 
#' \item posterior.dist:  Plot posterior distributions. 
#' \item fitted.observation:  Plot fitted curve versus observational data. 
#' }
#'
#' @export
#'
plot_fit_result <- function(fitobj, ci=0.95) {
    
    # Plot statistical errors:
    error = plot_abc_err(fit = fitobj$fit, 
                         prm.abc = fitobj$prm.abc)
    
    # Plot statistical summary of the posterior
    statistic.summary = fitobj$post.ss$plot
    
    # Plot posterior distributions:
    posterior.dist = plot_post_distrib_abc(fitobj$post.abc, fitobj$df.prior)
    
    # Check the fit against the data
    fitted.observation = plot_fit_abc_vs_obs(
      prm      = fitobj$prm, 
      post.abc = fitobj$post.abc, 
      obs.long = fitobj$obs.long, 
      hosp.var = fitobj$hosp.var,
      case.var = fitobj$case.var,
      ci       = ci)
    
    return(list(error              = error,
                statistic.summary  = statistic.summary,
                posterior.dist     = posterior.dist,
                fitted.observation = fitted.observation))
}

#' @title Visualize Fitting Errors
#'
#' @param fit List. Fitted objects. Output of function fit() 
#' @param prm.abc List. Required variables for abc fitting.
#'
#' @return
#' @export
#'
plot_abc_err <- function(fit, prm.abc) {
    dferr = fit$err
    
    npost = prm.abc$n * prm.abc$accept
    
    # Thin beyond threshold for lighter plot
    nremain = nrow(dferr) - npost
    if(nremain > 1000){
        idx = round(seq(npost+1, nrow(dferr), length.out = 500))
        dferr = dferr[c(1:npost, idx),]
    }
    
    if('fit.period' %in% names(dferr)){
      dferr = select(dferr, -fit.period)
    }
    
    g = dferr %>%
        mutate(x=1:nrow(dferr)) %>%
        select(-i) %>%
        pivot_longer(cols = -x) %>%
        filter(!grepl('total',name)) %>%
        ggplot(aes(x=x, y=value)) + 
        # geom_area(aes(fill=name)) +
        geom_step(aes(color=name))+
        geom_vline(xintercept = prm.abc$n * prm.abc$accept, linetype='dashed')+
        theme(panel.grid.minor = element_blank(),
              axis.ticks = element_blank())+
        scale_x_log10()+
        scale_y_log10()+
        scale_fill_brewer(palette = 'Pastel1')+
        xlab('rank') + ylab('error')+
        ggtitle('Ranked ABC errors',
                paste0('Posteriors: best ',
                       npost,
                       ' (',round(prm.abc$accept*100,2),'%)'))
    
    
    g2 = dferr %>% 
        slice_head(n=npost) %>% 
        pivot_longer(cols=-i) %>%
        ggplot(aes(x=name, y=value)) +
        geom_boxplot()+
        theme(panel.grid = element_blank()) +
        xlab('') + ylab('')+
        ggtitle('Errors from posteriors only')
    
    
    return(g | g2)
}

#' @title Visualize posterior distribution
#'
#' @param post.abc Posterior results. Output of function fit()
#' @param df.prior Dataframe. Prior parameters for fitting. 
#' @param df.true 
#'
#' @return
#' @export
#'
plot_post_distrib_abc <- function(post.abc, 
                                  df.prior, 
                                  df.true = NULL) {
    
    post.abc.long = post.abc %>%
        pivot_longer(cols = 1:ncol(post.abc)) 
    
    post.abc.ss = post.abc.long %>%
        group_by(name) %>%
        summarise(m = mean(value))
    
    col.post = 'tomato3'
    
    df = post.abc.long
    
    if(!is.null(df.true)) 
        df = left_join(post.abc.long, df.true, by='name')
    
    g = df %>%
        ggplot(aes(x = value)) + 
        geom_histogram(data=df.prior, fill='lightgrey', aes(y=..density..),
                       bins = 30) +
        geom_histogram(bins = 20, aes(y=..density..), fill=col.post, alpha=0.3) + 
        geom_density(size=0.7, col=col.post) +
        geom_vline(data= post.abc.ss, aes(xintercept = m), col=col.post) + 
        facet_wrap(~name, scales = 'free')+
        theme(strip.text = element_text(face='bold'),
              panel.grid = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) + 
        xlab('') + ylab('') + 
        ggtitle('Posterior distributions')
    
    if(!is.null(df.true)) 
        g = g + geom_vline(aes(xintercept = true.value), linetype = 'dashed')
    
    return(g)
}


#' @title Visualize fitting performance versus observational data
#'
#' @param prm List. All parameters for model found in \code{wem-prm.csv}.
#' @param post.abc List. Contains posterior results. Output of function fit()
#' @param obs.long Dataframe. Observational data
#' @param hosp.var String. Type of hospital (e.g., \code{NULL}, \code{'hosp.adm'}, \code{'hosp.occ'})
#' @param case.var String. Type of date for clinical cases (e.g., \code{'report'} and \code{'episode'})
#' @param ci Numerical. Percentage of confidence interval.
#'
#' @return
#' @export
#'
plot_fit_abc_vs_obs <- function(prm, 
                                post.abc, 
                                obs.long,
                                hosp.var,
                                case.var,
                                ci = 0.95) {
    
    message('Plotting ABC fit vs. observations ...')
    
    # Simulate epidemics from 
    # posterior distributions: 
    npf    = names(post.abc)
    n.post = nrow(post.abc)
    tmp = list()
    for(i in 1:n.post){
        if(i%%10==0) message(paste0(i,'/',n.post))
        prmi = prm
        for(j in seq_along(npf)){
            prmi[[npf[j]]] <- post.abc[[ npf[j] ]][i]
        }
        tmp[[i]] = simul(prm = prmi)$ts
        tmp[[i]]$iter = i
    }
    # Dataframe of simulation
    df.post = bind_rows(tmp)
    
    # Add dates
    df.post$date = lubridate::ymd(min(obs.long$date) + df.post$time)
    
    # Determine type of date for clinical cases (reported date or episode date)
    if(case.var=='report'){
      df.post = df.post %>%
        mutate(clin.case = report)
    }
    if(case.var=='episode'){
      df.post = df.post %>%
        mutate(clin.case = report.episode)
    }
    
    if(!is.null(hosp.var)){
      #---- Summary stats of posterior simulations:
      # Determine hospital type (new admissions or occupancy)
      if(hosp.var=='hosp.adm'){
          df.post = df.post %>%
              mutate(hospital = hosp.admission)
      }
      if(hosp.var=='hosp.occ'){
          df.post = df.post %>%
              mutate(hospital = Hall)
      }
        
      df.ss = df.post %>%
          group_by(date) %>%
          summarize(report.m = mean(clin.case), 
                    report.qhi = quantile(clin.case, probs = 0.5 + ci/2),
                    report.qlo = quantile(clin.case, probs = 0.5 - ci/2),
                    ww.m = mean(WWreport), 
                    ww.qhi = quantile(WWreport, probs = 0.5 + ci/2),
                    ww.qlo = quantile(WWreport, probs = 0.5 - ci/2),
                    hosp.m = mean(hospital),
                    hosp.qhi = quantile(hospital, probs = 0.5 + ci/2),
                    hosp.qlo = quantile(hospital, probs = 0.5 - ci/2),
                    .groups = 'keep')
    }
    
    if(is.null(hosp.var)){
        df.ss = df.post %>%
            group_by(date) %>%
            summarize(report.m = mean(clin.case), 
                      report.qhi = quantile(clin.case, probs = 0.5 + ci/2),
                      report.qlo = quantile(clin.case, probs = 0.5 - ci/2),
                      ww.m = mean(WWreport), 
                      ww.qhi = quantile(WWreport, probs = 0.5 + ci/2),
                      ww.qlo = quantile(WWreport, probs = 0.5 - ci/2),
                      .groups = 'keep')
    }
    
    df.ss.long = df.ss %>% 
        pivot_longer(cols = -date) %>%
        mutate(type = ifelse(grepl('report',name),'clin',
                             ifelse(grepl('hosp', name), 'hosp','ww'))) %>%
        mutate(stat = ifelse(grepl('qhi',name),'qhi',
                             ifelse(grepl('qlo',name),'qlo','m'))) %>%
        select(-name) %>%
        pivot_wider(names_from = 'stat')
    
    obs.long$type = stringr::str_remove(obs.long$name,'\\.obs')
    
    
    # --- Plots ---
    
    g.chk = df.ss.long %>%
        ggplot(aes(x=date, color=type, fill=type)) + 
        #
        # observations:
        geom_point(data=obs.long, aes(x=date, y=value), size=1, alpha=0.8) +
        geom_step(data=obs.long, aes(x=date, y=value), size=0.2) +
        #
        # posterior:
        geom_line(aes(y = m), size=1) + 
        geom_ribbon(aes(ymin=qlo, ymax=qhi), alpha=0.2, size=0.2) + 
        facet_wrap(~type, scales = 'free_y', ncol=1) +
        theme(panel.grid.minor = element_blank())+
        scale_x_date(date_breaks = '2 months', date_labels = '%b \'%y')+
        guides(fill='none', color='none')+
        ggtitle(paste('Check fit')) + 
        xlab('') + ylab('')
    
    message(paste0('\nPlot ABC fit vs. observations done ','.\n'))
    
    return(g.chk)
}
