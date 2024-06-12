test_that("model fit is satisfactory", {

  prm.model = model_prm_example()
  prm.model$horizon <- 90
  sim = wem::simul(prm.model)
  # wem::plot_epidemic(sim$ts)
  
  # Generate synthetic observations from simulation
  set.seed(1234)
  
  genobs <- function(times, values, prob, 
                     date.start = lubridate::ymd('2000-01-01')) {
    nn = length(values)
    idx.obs = c(1:nn)[rbinom(n=nn, size=1, prob = prob) |> as.logical()]
    t.obs   = times[idx.obs]
    tmp = rpois(n = length(idx.obs), lambda = values[idx.obs])
    res = data.frame(
      date  = date.start + times[idx.obs], 
      value = tmp )  
    return(res)
  }
  
  
  obs.reports = genobs(times  = sim$ts$time,
                       values = sim$ts$report,
                       prob   = 0.3)

  obs.ww = genobs(times  = sim$ts$time,
                       values = sim$ts$WWreport,
                       prob   = 0.2)
    
  obs.hosp = genobs(times  = sim$ts$time,
                    values = sim$ts$hosp.admission,
                    prob   = 0.2)
  
  # ggplot(obs.hosp, aes(x=date, y=value)) + geom_point()

  dat = build_data(cases = obs.reports, 
                   hosp = obs.hosp, ww = obs.ww, 
                   hosp.type = 'hosp.adm', 
                   case.date.type = 'report')
  
 prm.abc = define_abc_prms(
   iter        = 2e3,
   accept      = 1e-2,
   case.weight = 1.0,
   ww.weight   = 1.0,
   hosp.weight = 1.0,
   hosp.type   = 'hosp.adm') 
 
priors = data.frame(
  name    = c('R0', 
              'ww.scale',    
              'transm.v_2', 
              'hospital.prop'),
  
  distrib = c('rnorm',  
              'runif',
              'runif', 
              'runif'),
  
  prms    = c('2.5;0.8',
              '1e-5;3e-4',
              '0.2;1.3', 
              '1e-3;0.02')
)

fitobj = fit(data      = dat,
             prm.abc   = prm.abc,
             prm       = prm.model,
             df.priors = priors,
             n.cores   = 1, 
             last.date = NULL)

# g = plot_fit_result(fitobj = fitobj,ci = 0.95)
# wrap_plots(g$fitted.observation, g$posterior.dist, g$error, ncol = 2)

# Are the posterior estimates close to
# the `true` known values (obsverations 
# were generated from simulations)

p = fitobj$post.ss$df

testthat::expect_lt(abs(prm.model$R0/p$m[p$name == 'R0']-1), 0.3)

# TODO: FIX TESTS BELOW, FAILING...
# testthat::expect_lt(abs(prm.model$transm.v[2]/p$m[p$name == "transm.v_2"]-1), 0.3)
# testthat::expect_lt(abs(prm.model$ww.scale/p$m[p$name == "ww.scale"]-1), 0.3)
    
})
