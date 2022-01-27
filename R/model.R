

#' SEIR compartmental model with `nI` I, `nIH` IH, `nA` A, `nH` H and `nZ` Z compartments.
#' @param t Numerical. Simulation time. 
#' @param x List. Initial epidemiological states.   x = (S, Vw, V, Ev1, Ev2, ... nEv,
#' @param parms List. Model parameters.
#' 
seir <- function(t, x, parms)
{
    # Unpack `parms`
    # Not using `with()` because it's faster without:
    nA = parms$nA
    nE = parms$nE
    nEv = parms$nEv
    nI = parms$nI
    nIH = parms$nIH
    nH = parms$nH
    nZ = parms$nZ
    hosp.prop = parms$hosp.prop
    asymp.prop = parms$asymp.prop
    h.vac = parms$h.vac
    alpha.vac = parms$alpha.vac
    delta = parms$delta
    beta = parms$beta
    nepsilon = parms$nepsilon
    nepsilon.vac = parms$nepsilon.vac
    ntau = parms$ntau
    nmu = parms$nmu
    ntheta = parms$ntheta
    neta = parms$neta
    nell = parms$nell
    popSize = parms$popSize
    transm.t = parms$transm.t 
    transm.v = parms$transm.v
    vacc.rate.t = parms$vacc.rate.t 
    vacc.rate.v = parms$vacc.rate.v
    hosp.prop.t = parms$hosp.prop.t 
    hosp.prop.v = parms$hosp.prop.v
    asymp.prop.t = parms$asymp.prop.t 
    asymp.prop.v = parms$asymp.prop.v
    hosp.rate.vacc.t = parms$hosp.rate.vacc.t 
    hosp.rate.vacc.v = parms$hosp.rate.vacc.v
    asymp.prop.vacc.t = parms$asymp.prop.vacc.t 
    asymp.prop.vacc.v = parms$asymp.prop.vacc.v
    inf.A = parms$inf.A 
    inf.I = parms$inf.I
    inf.IH = parms$inf.IH
    rel.inf.a = parms$rel.inf.a
    eff.inf = parms$eff.inf
    r = parms$r 
    d = parms$d 
    tau.immu.R = parms$tau.immu.R
    tau.immu.V = parms$tau.immu.V
    tau.immu.R.t = parms$tau.immu.R.t
    tau.immu.V.t = parms$tau.immu.V.t
    tau.immu.R.v = parms$tau.immu.R.v
    tau.immu.V.v = parms$tau.immu.V.v
    
    # define compartment parameters in vector x
    
    n_E_Ev_I = nEv + nE + nI
    n_E_IH   = n_E_Ev_I + nIH
    n_E_A    = n_E_IH + nA
    n_E_Z    = n_E_A + nH + nZ
    
    S  = x[1]
    Vw = x[2]
    V  = x[3]
    E  = x[4:(nE+3)]
    Ev = x[(nE+4):(nEv+nE+3)]
    I  = x[(nEv+nE+4):(n_E_Ev_I+3)]
    IH = x[(n_E_Ev_I+4):(n_E_IH+3)]
    A  = x[(n_E_IH+4):(n_E_A+3)]
    H  = x[(n_E_A+4):(n_E_A+nH+3)]
    Z  = x[(n_E_A+nH+4):(n_E_Z+3)]
    R  = x[n_E_Z+4] 
    D  = x[n_E_Z+5]
    cuminc     = x[n_E_Z+6]
    cumincsymp = x[n_E_Z+7]
    cumHospAdm = x[n_E_Z+8]
    
    # Run a few checks to catch obvious user entry errors
    docheck = TRUE
    if(docheck & (t < 1)){
        stopifnot(nE == length(E))
        stopifnot(nEv == length(Ev))
        stopifnot(nA == length(A))
        stopifnot(nI == length(I))
        stopifnot(nIH == length(IH))
        stopifnot(nH == length(H))
        
        stopifnot(length(A) == length(inf.A))
        stopifnot(length(I) == length(inf.I))
        stopifnot(length(IH) == length(inf.IH))
        
        stopifnot(grepl('^E', names(E)))
        stopifnot(grepl('^Ev', names(Ev)))
        stopifnot(grepl('^A', names(A)))
        stopifnot(grepl('^I', names(I)))
        stopifnot(grepl('^IH', names(IH)))
        stopifnot(grepl('^H', names(H)))
        stopifnot(grepl('^Z', names(Z)))
    }
    
    # The "d" in front of the variable name "S" 
    # means we are now considering an differential equation
    # "dS" <=> 'dS(t) / dt'
    
    # --- Transmission rate 
    
    # normalize infectiousness profile:
    inf.A.norm  = inf.A / sum(inf.A) * length(inf.A)
    inf.I.norm  = inf.I / sum(inf.I) * length(inf.I)
    
    sum.A  = sum(inf.A.norm * A)
    sum.I  = sum(inf.I.norm * I)
    sum.IH = sum(inf.I.norm[1:nIH] * IH)
    
    if(!is.numeric(transm.v)){
        beta_t = beta
    } 
    else{
        # Apply time-dependent transmission rate
        mult   = broken_line(x=t, b = transm.t, v = transm.v)
        beta_t = beta * mult
    }
    
    if(!is.numeric(vacc.rate.v)){
        rt = r 
    }else{
        # Apply time-dependent vaccine rate
        rt = broken_line(x=t, b=vacc.rate.t, v=vacc.rate.v)
    }
    
    #--- time-dependent immunity
    if(!is.numeric(tau.immu.R.v)){
        tau.immu.R_t = tau.immu.R
    }else{
        tau.immu.R_t = broken_line(x=t,
                                   b=tau.immu.R.t,
                                   v=tau.immu.R.v)
    }
    if(!is.numeric(tau.immu.V.v)){
        tau.immu.V_t = tau.immu.V
    }else{
        tau.immu.V_t = broken_line(x=t,
                                   b=tau.immu.V.t,
                                   v=tau.immu.V.v)
    }
    #---
    
    # calculate incidence
    infrate  = beta_t * S * (rel.inf.a  *sum.A + sum.I + sum.IH) / popSize
    
    ## susceptible
    dS = tau.immu.R_t*R + tau.immu.V_t*V - rt*S - infrate
    
    ## vaccinated but not immuned
    dVw = rt*S - d * Vw
    
    ## vaccinated with full protection 
    vac.infrate = (1-eff.inf) * beta_t * V * (rel.inf.a  *sum.A + sum.I + sum.IH) / popSize
    dV = d*Vw - vac.infrate - tau.immu.V_t*V
    
    ## Exposed
    if(nE == 1) dE = infrate - nepsilon * E
    
    if(nE > 1){
        dE = nepsilon * ( c(0, E[1:(nE-1)]) - E )
        dE[1] = dE[1] + infrate
    }
    
    ## Vaccinated and exposed 
    if(nEv == 1) dEv = vac.infrate - nepsilon.vac * Ev
    
    if(nEv > 1){
        dEv = nepsilon.vac * ( c(0, Ev[1:(nEv-1)]) - Ev )
        dEv[1] = dEv[1] + vac.infrate
    }
    
    if(!is.numeric(hosp.prop.v)){
        h_t = hosp.prop
    }else{
        # Time-dependent hospital proportion
        h_t = broken_line(x=t,
                          b = hosp.prop.t,
                          v = hosp.prop.v)
    }
    
    if(!is.numeric(asymp.prop.v)){
        alpha.t = asymp.prop 
    }else{
        # Time dependent asymptomatic proportion
        alpha.t = broken_line(x=t, 
                                 b = asymp.prop.t, 
                                 v = asymp.prop.v)
    }
    
    if(is.null(hosp.rate.vacc.v)){
        h.vac_t = h.vac
    }else{
        # Vaccinated Time-dependent hospital rate
        h.vac_t = broken_line(x = t, 
                              b = hosp.rate.vacc.t, 
                              v = hosp.rate.vacc.v)
    }
    
    if(is.null(asymp.prop.vacc.v)){
        alpha.vac.t    =  alpha.vac
    }else{
        # Vaccinated Time dependent asymptomatic proportion
        alpha.vac.t = broken_line(x = t, 
                                  b = asymp.prop.vacc.t,
                                  v = asymp.prop.vacc.v)
    }
    
    # print(paste('DEBUG alpha.t =', alpha.t, 't =',t))
    
    ## infectious (symp + sympHosp + asymp)
    if(nI == 1) dI = (1-h_t)*(1-alpha.t) * nepsilon * E[nE] + (1-h.vac_t)*(1-alpha.vac.t) * nepsilon.vac * Ev[nEv]- ntau * I
    if(nIH == 1) dIH = h_t*(1-alpha.t) * nepsilon * E[nE] + h.vac_t*(1-alpha.vac.t) * nepsilon.vac * Ev[nEv]- nmu * IH
    if(nA == 1) dA = alpha.t * nepsilon * E[nE] + alpha.vac.t * nepsilon.vac * Ev[nEv] - ntheta * A
    
    if(nI > 1) {
        dI = ntau * (c(0, I[1:(nI-1)]) - I)
        dI[1] = dI[1] + (1-h_t)*(1-alpha.t) * nepsilon * E[nE] + (1-h.vac_t)*(1-alpha.vac.t) * nepsilon.vac * Ev[nEv]
    }
    if(nIH > 1) {
        dIH = nmu * (c(0, IH[1:(nIH-1)]) - IH)
        dIH[1] = dIH[1] + h_t*(1-alpha.t) * nepsilon * E[nE] + h.vac_t*(1-alpha.vac.t) * nepsilon.vac * Ev[nEv]
    }
    if(nA > 1){
        dA = ntheta * (c(0, A[1:(nA-1)]) - A)   
        dA[1] = dA[1] + alpha.t * nepsilon * E[nE] + alpha.vac.t * nepsilon.vac * Ev[nEv]
    } 
    
    
    ## Hospital admission and occupancy
    
    hospadm =  nmu * IH[nIH]
    
    if(nH == 1) dH = hospadm - nell * H
    
    if(nH > 1){
        dH = nell * ( c(0, H[1:(nH-1)]) - H )
        dH[1] = dH[1] + hospadm
    }
    
    ## Cumulative hospital ADMISSIONS
    dcumHospAdm = hospadm
    
    ## not infectious but shedding RNA copies to WW
    if(nZ == 1) dZ = ntau*I[nI] + ntheta*A[nA] - neta * Z
    
    if(nZ > 1){
        dZ = neta * (c(0,Z[1:nZ-1]) - Z)
        dZ[1] = dZ[1] + ntau*I[nI] + ntheta*A[nA]
    }
    
    ## recovered
    dR = neta * Z[nZ] + (1-delta) *  nell * H[nH] - tau.immu.R_t*R
    
    ## death
    dD = delta * nell * H[nH]
    
    ## Cumulative infection (symp+sympHosp+asymp) incidence
    dcuminc     = infrate + vac.infrate
    dcumincsymp = infrate * (1-alpha.t) + vac.infrate * (1-alpha.vac.t)
    
    ## IMPORTANT: all the derivatives must be 
    ## saved into a big list which is as long as "x"
    res=c(dS, 
          dVw, dV, 
          dE, dEv, 
          dI, dIH, dA, dH, 
          dZ, dR, dD, 
          dcuminc, dcumincsymp, dcumHospAdm)
    
    return(list(res))
}
