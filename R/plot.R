###
###  PLOTTING FUNCTIONS
###



#' @title Visualize model parameters
#' 
#' @description Provide a synthetic view of the model parameter in a dashboard-type format.
#'
#' @param prm List of all model parameters as returned by the function \code{wem::load_prm()}.
#' @param nrow Integer. Number of rows in the facetted plot (default = 3)
#' @param textsize Intager. Font size for the text in the plot (default = 10)
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @import ggplot2
#'
plot_dashboard_prm <- function(prm, nrow = 3, textsize = 10) {

  .mkdf <- function(prm, regex) {
    tmp = unlist(prm[grepl(regex, names(prm))])
    df = data.frame(var= names(tmp), value = tmp)
    return(df)
  }

  # --- Durations
  g.dur = .mkdf(prm, '^dur')  %>%
    ggplot(aes(x=var, y=value)) + geom_bar(stat='identity') +
    xlab('') + ylab('') + ggtitle('Durations')+
    theme(panel.grid.major.y = element_blank(),
          text = element_text(size = textsize))+
    coord_flip()

  # --- Proportions
  g.prop = .mkdf(prm, 'prop$')  %>%
    ggplot(aes(x=var, y=value)) + geom_bar(stat='identity') +
    xlab('') + ylab('') + ggtitle('Proportions')+
    theme(panel.grid.major.y = element_blank(),
          text = element_text(size = textsize))+
    coord_flip(ylim = c(0,1)) +
    scale_y_continuous(breaks = seq(0,1,by=0.2))

  # --- Time lags
  g.lag = .mkdf(prm, '\\.lag')  %>%
    ggplot(aes(x=var, y=value)) + geom_bar(stat='identity') +
    xlab('') + ylab('days') + ggtitle('Time lags')+
    theme(panel.grid.major.y = element_blank(),
          text = element_text(size = textsize))+
    coord_flip()

  # --- Number of subcompartments
  g.n = .mkdf(prm, '^n[A-Z]+$')  %>%
    ggplot(aes(x=var, y=value)) + geom_bar(stat='identity') +
    theme(panel.grid = element_blank(),
          text = element_text(size = textsize))+
    xlab('') + ylab('number') + ggtitle('Number of subcompartments')

  # --- Infectiousness
  g.inf =  .mkdf(prm, '^inf') %>%
    mutate(v = stringr::str_extract(var, 'inf\\.[A-Z]+'),
           x = as.numeric(stringr::str_extract(var, '\\d+$'))) %>%
    ggplot(aes(x=x, y=value)) +
    geom_point(color = 'tomato3', alpha=0.5)+
    geom_area(fill = 'tomato3', alpha=0.5)+
      theme(text = element_text(size = textsize))+
    xlab('Days since infection') + ylab('') + ggtitle('Infectiousness profile')+
    facet_grid(~v)

  # --- Fecal shedding
  g.shed =  .mkdf(prm, '^shed') %>%
    mutate(v = stringr::str_extract(var, 'shed\\.[A-Z]+'),
           x = stringr::str_extract(var, '\\d*$')) %>%
    mutate(x = as.numeric(ifelse(x=='', 1,x)),
           v = stringr::str_remove(v,'^shed.')) %>%
    ggplot(aes(x=x, y=value)) +
    geom_point(color = 'orange3', alpha=0.5)+
    geom_area(fill = 'orange3', alpha=0.5)+
    xlab('Days in state') + ylab('') +
    scale_y_log10()+
      theme(text = element_text(size = textsize))+
    ggtitle('Fecal shedding profile')+
    facet_wrap(~v, nrow=2)

  # --- Time-dependent interventions
  tmp = .mkdf(prm, '^transm.[tv]') %>%
    mutate(v = stringr::str_extract(var, '^transm.[tv]'),
           x = as.numeric(stringr::str_extract(var, '\\d*$')))
  b = filter(tmp, v=='transm.t')
  v = filter(tmp, v=='transm.v')
  dj = left_join(b,v,by='x')

  g.int = dj %>%
    ggplot(aes(x=value.x, y=value.y)) +
    geom_line(color='red2')+geom_point(color='red2')+
    xlab('time') + ylab('') +
      theme(text = element_text(size = textsize)) +
    ggtitle('Multiplier on transmission')

  # --- Time-dependent hospitalization rate
  tmp = .mkdf(prm, '^hosp.rate.[tv]') %>%
    mutate(v = stringr::str_extract(var, '^hosp.rate.[tv]$'),
           x = as.numeric(stringr::str_extract(var, '\\d*$')))
  b = filter(tmp, v=='hosp.rate.t')
  v = filter(tmp, v=='hosp.rate.v')
  dj = left_join(b,v,by='x')

  g.hosp = dj %>%
    ggplot(aes(x=value.x, y=value.y)) +
    geom_step(color='green3',size=2)+geom_point(color='green3')+
    xlab('time') + ylab('') +
      theme(text = element_text(size = textsize)) +
    ggtitle('Multiplier on hosp. rate')

  # --- Time-dependent asymptomatic fraction
  tmp = .mkdf(prm, '^asymp.prop.[tv]$') %>%
    mutate(v = stringr::str_extract(var, '^asymp.prop.[tv]$'),
           x = as.numeric(stringr::str_extract(var, '\\d*$')))
  b = filter(tmp, v=='asymp.prop.t')
  v = filter(tmp, v=='asymp.prop.v')
  dj = left_join(b,v,by='x')

  g.asymp = dj %>%
    ggplot(aes(x=value.x, y=value.y)) +
    geom_step(color='steelblue2',size=2)+geom_point(color='green3')+
    xlab('time') + ylab('') +
      theme(text = element_text(size = textsize)) +
    ggtitle('Multiplier on Asymp. fraction')

  # --- Dashboard
  patchwork::wrap_plots(
    g.dur,
    g.prop,
    g.lag,
    g.n,
    g.inf,
    g.shed,
    g.int,
    g.hosp,
    g.asymp,
    nrow = nrow)
}

#' @title Plot all compartments of the model
#' 
#' 
#' @param df Dataframe of time series as returned by \code{simul()}.
#' @param ncol Integer. Number of column in the faceted plot (default=6).
#' @param textsize Integer. Size of the text for the faceted plot (default=12).
#'
#' @return A \code{ggplot2} object displaying a facet plot for all time series. 
#' 
#' @export
#'
plot_all <- function(df, ncol = 6, textsize=12) {
  df.long = pivot_longer(df, -time)

  g = df.long %>%
    ggplot(aes(x=time, y=value)) +
    geom_line(size=1) +
      theme(text=element_text(size=textsize)) + 
    facet_wrap(~name, scales = 'free_y', ncol=ncol)
  return(g)
}


#' @title Plot time series of main epidemiological variables.
#'
#' @param df Dataframe of time series as returned by \code{simul()}.
#'
#' @return A \code{ggplot2} object displaying a facet plot.
#' @export
#'
plot_epidemic <- function(df) {
  # df = sim$ts
  df.long = pivot_longer(df, -time)
  
  g = df.long %>%
    filter(name %in% c('S','inc', 'V','hosp.admission', 
                       'report', 'R')) %>%
    ggplot(aes(x=time, y=value)) +
    geom_line(size=1) +
    # theme(text=element_text(size=textsize)) + 
    facet_wrap(~name, scales = 'free_y', ncol=3)
  return(g)
}


#' @title Plot the Effective Reproduction Number
#'
#' @description Plot the effective reproduction number as calculated by
#' the function \code{estimate_Rt()} 
#'
#' @param Rdata Dataframe. Rt values from posterior simulation. Output of function \code{estimate_Rt()}
#' @param min.date String. Minimum date to show Rt values on the plot.
#' @param max.date String. Maximum date to show Rt values on the plot.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
plot_Rt <- function(Rdata,min.date,max.date){

  df = Rdata 

  # minimum and maximum date
  date.min = as.Date(min.date)
  date.max = as.Date(max.date)
  
  
  # characteristics of x-axis 
  xaxis = scale_x_date(date_breaks = '2 months',
                       date_labels = '%b \'%y',
                       limits = c(date.min,date.max))
  

  # --- Plot for the mean estimates
  
  g = df %>%
    ggplot()+
    geom_line(aes(x=date, y=Reff.m), colour ='blue') +
    geom_ribbon(aes(x=date,ymin=Reff.lo, ymax=Reff.hi), alpha=0.2, fill='steelblue2') +
    geom_hline(yintercept = 1, linetype='dashed', color='grey50')+
    xaxis+
    xlab('')+
    ylab('Effective reproduction number')
    
  return(g)
}

#' @title Plot forecasts
#'
#' @param fcst.obj 
#' @param var 
#' @param dat 
#'
#' @return
#' @export
#'
plot_fcst <- function(var, fcst.obj, dat, 
                      date_breaks = '1 month',
                      show.fit = TRUE) {
  
  simpost = fcst.obj$sim.post
  names(simpost)
  
  d.last = max(dat$obs$date)
  
  simpost.before = filter(simpost, date <= d.last)
  simpost.after  = filter(simpost, date > d.last)
  
  var.m  = paste0(var,'.m')
  var.lo = paste0(var,'.lo')
  var.hi = paste0(var,'.hi')
  
  if(var == 'report') var.obs = 'clin.obs'
  if(var == 'ww')     var.obs = 'ww.obs'
  if(var == 'hosp')   var.obs = 'hosp.obs'
  
  col.fcst = 'steelblue2'
  
  g = ggplot(mapping = aes(x=date)) + 
    geom_point(data = dat$obs, aes_string(y=var.obs)) + 
    geom_line(data = simpost.before, aes_string(y=var.m), 
              linetype='dashed', color = col.fcst)+
    geom_line(data = simpost.after, aes_string(y=var.m),
              color = col.fcst) +
    geom_ribbon(data = simpost.after, 
                mapping = aes_string(ymin=var.lo, ymax=var.hi), 
                alpha = 0.2, fill = col.fcst)+
    xlab('') + ylab(var)
  
  g = g + scale_x_date(date_breaks = date_breaks, 
                       date_labels = '%Y-%b-%d')
  # g
  return(g)
}



#' @title Visualize observational data, simulated epidemic and transmission rate.
#' 
#' @description Plot the observed data and the simulated epidemic generated 
#' from the initial parameter values. Also plot the time-varying transmission rate.
#' This function is typically used to check an initial guess for priors 
#' (e.g., transmission rate (\code{beta}), wastewater scaling factor (\code{ww.scale})) 
#' before launching the procedure that fits selected model parameters to data.
#'
#' @param data Data object (list) as returned by \code{wem::build_data()}
#' @param prm List of all parameter values.
#' @param include.cases Logical. Display observations of reported cases?
#' @param include.ww  Logical. Display observed viral concentration in wastewater? 
#' @param include.hosp Logical. display hospital observations?
#' @param log.scale Logical. Use log scale?
#'
#' @return A \code{ggplot2} object.
#' 
#' @export
#'
#' @examples
#' # Load data sets examples from `wem` package:
#' data('cases')
#' data('hosp')
#' data('wwviralconc')
#' 
#' # Build the data object:
#' dat = build_data(cases = cases, 
#'                  hosp = hosp, 
#'                  ww = wwviralconc, 
#'                  hosp.type = 'hosp.adm', 
#'                  case.date.type = 'report')
#'
#' # Load example of mdoel parameters
#' prm = model_prm_example()
#' 
#' # Plot data, simulation and 
#' # time-dependent transmission rate (beta)
#' g = plot_simobs_beta(data          = dat,
#'                      prm           = prm,
#'                      include.cases = TRUE,
#'                      include.ww    = TRUE,
#'                      include.hosp  = FALSE,
#'                      log.scale     = FALSE)
#' plot(g)
#'
plot_simobs_beta <- function(data,
                             prm,
                             include.cases = TRUE,
                             include.ww    = TRUE,
                             include.hosp  = FALSE,
                             log.scale     = FALSE){
  
  # rename/separate data
  obs      = data[['obs']]
  obs.long = data[['obs.long']]
  hosp.var = data[['hosp.var']]
  case.var = data[['case.var']]
  
  
  # time range
  xrng = range(prm[['transm.t']], obs.long$time)
  xaxis = scale_x_continuous(limits = xrng)
  
  # Plot transmission rate (beta) 
  d0 = min(obs$date)
  g.beta.init = plot_multvec(prm, 
                             xname = 'transm.t',
                             yname = 'transm.v',
                             d0=d0) + 
    xaxis + 
    xlab('time') + 
    theme(panel.grid.minor.y = element_blank())
  
  
  # -- Single inital simulation 
  sim.init = simul(prm)
  df.init  = sim.init$ts
  
  # plot single simulation with initial prms
  g.init = plot_obs_simulation(df.init, obs.long, 
                               dates.break   = prm$transm.t,
                               hosp.var      = hosp.var,
                               case.var      = case.var,
                               include.cases = include.cases,
                               include.ww    = include.ww,
                               include.hp    = include.hosp,
                               logscale      = log.scale) + 
    xaxis + 
    xlab('') + 
    theme(plot.margin = margin(b=0))
  
  g.final = patchwork::wrap_plots(g.init, g.beta.init, ncol = 1)
  
  return(g.final)
}



#' @title Plot Observational and Simulated Data
#' 
#' @description For each wastewater sampling location, it plots observational data and corresponding simulation for
#' reported cases, hospital admissions, hospital occupancy, and measured RNA concentration in wastewater
#'
#' @param df.init Dataframe. Time-series simulation of all wwmodel's compartments  
#' @param obs.long Dataframe. Time-series observational data in a long format created by \code{build_obs()}
#' @param dates.break Logical. Default is \code{TRUE}. Display break dates when effective reproduction number changes
#' @param hosp.var String. Type of hospitalization (e.g., \code{NULL}, \code{'hosp.adm'}, \code{'hos.occ'})
#' @param case.var String. Type of date for clinical cases (e.g., \code{'report'} and \code{'episode'})
#' @param include.cases Logical.Default is \code{TRUE}. Display simulated and observational reported cases
#' @param include.ww Logical. Default is \code{TRUE}. Display simulated and observational measured RNA concentration  
#' @param include.hp Logical. Default is \code{TRUE}. Display simulated and observational hospitalization
#' @param logscale Logical. Default is \code{FALSE}. Logarithmic scale for x-axis
#'
#' @return Visualization for simulation and observational data 
#' @export
#'
plot_obs_simulation <- function(df.init,
                                obs.long,
                                dates.break,
                                hosp.var,
                                case.var,
                                include.cases = TRUE,
                                include.ww = TRUE,
                                include.hp = TRUE,  
                                logscale = FALSE) {
  
  df.long = pivot_longer(df.init, -time)
  
  if(case.var == 'report')simplot = filter(df.long,name %in% c('report','WWreport'))
  if(case.var == 'episode')simplot = filter(df.long,name %in% c('report.episode','WWreport'))
  
  if(!is.null(hosp.var)){
    if(hosp.var == 'hosp.adm' && case.var == 'report') simplot = filter(df.long,name %in% c('report','WWreport','hosp.admission'))
    if(hosp.var == 'hosp.occ' && case.var == 'report') simplot = filter(df.long,name %in% c('report','WWreport','Hall')) 
    if(hosp.var == 'hosp.adm' && case.var == 'episode') simplot = filter(df.long,name %in% c('report.episode','WWreport','hosp.admission'))
    if(hosp.var == 'hosp.occ' && case.var == 'episode') simplot = filter(df.long,name %in% c('report.episode','WWreport','Hall')) 
  }
  
  if(!include.cases){
    simplot = filter(simplot, !grepl('^report',name))
    obs.long = filter(obs.long, !grepl('clin.obs',name))
  }
  
  if(!include.ww){
    simplot = filter(simplot, !grepl('WWreport',name))
    obs.long = filter(obs.long, !grepl('ww.obs',name))
  }
  
  if(!include.hp){
    simplot = filter(simplot, 
                     ! name %in% c('hosp.admission','Hall'))
    obs.long = filter(obs.long, !grepl('hosp.obs',name))
  }
  
  simplot  = .plottype(simplot)
  obs.long = .plottype(obs.long)
  
  g = simplot %>%
    ggplot(aes(x=time, y=value, colour=name)) +
    geom_step(data = obs.long, size=0.4)+
    geom_point(data = obs.long, size=1, fill='white', shape=21)+
    geom_line(size = 1)+
    scale_color_brewer(palette= 'Set2')+
    theme(panel.grid.minor = element_blank())+
    ggtitle('Simulation and observational data')
  
  g = g + facet_wrap(~plottype, scales='free_y', ncol=1)
  
  if(!is.null(dates.break)){
    g = g + geom_vline(xintercept = dates.break,
                       color = 'black', alpha=0.5, 
                       size=0.5, linetype = 'longdash')
  }
  
  if(logscale) 
    g = g + scale_y_log10()
  
  return(g)
}



#' @title Plot Transmission Rate (Beta)
#'
#' @param prm 
#' @param xname 
#' @param yname 
#' @param d0 
#' @param xmax 
#'
#' @return
#' @export
#'
plot_multvec <- function(prm, xname, yname, d0, xmax = NULL) {
  # DEBUG
  if(0){
    xname = 'transm.t'
    yname = 'transm.v'
  }
  
  dfplot = data.frame(
    x = prm[[xname]],
    y = prm[[yname]]
  ) %>%
    mutate(date = d0 + x,
           i = row_number()) %>%
    mutate(a = i%%2==0)
  
  ymax = max(dfplot$y)
  
  g = dfplot %>%
    ggplot(aes(x=x,y=y)) +
    geom_segment(aes(x=x,xend=x, y=y, yend=1.2*ymax*a),
                 linetype='dotted', color='grey50')+
    geom_line() +
    geom_point(shape=21,fill='white',stroke=1)+
    geom_label(aes(y= 1.2*ymax * a,
                   label=paste0(format(date,'%b-%d'),
                                ' \n ', x,' | ',y)),
               fill = 'wheat1', alpha=0.3, size=2,
               fontface='bold', color='grey50')+
    # geom_text(aes(x=x, label= format(date,'%b-%d')),y=1.1*ymax,
    #           angle=60, hjust=0, size=3)+
    scale_y_continuous(limits = c(-0.2, 1.5*ymax))+
    ggtitle('Transmission Rate') +
    xlab(xname) + ylab(yname)
  
  if(!is.null(xmax)) g = g + coord_cartesian(xlim=c(0,xmax))
  return(g)
}


#' @title Type of Datasource for plotting
#'
#' @param d 
#'
#' @return
#'
.plottype <- function(d) {
  res = d %>%  
    mutate(plottype = ifelse(name %in% c('hosp','Hall','hosp.obs'), 'Hospitalization',
                             ifelse(grepl('[wW][wW]',name), 'WW', 
                                    'Clinical')))
  return(res)
}

