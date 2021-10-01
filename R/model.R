

#' SEIR compartmental model with `nI` I, `nIH` IH, `nA` A, `nH` H and `nZ` Z compartments.
#' @param t Numerical. Simulation time. 
#' @param x List. Initial epidemiological states.   x = (S, Vw, V, Ev1, Ev2, ... nEv,
#'                                                       E1, E2, ... nE, I1, I2, ..., nI, 
#'                                                       IH1, IH2, ..., nIH,A1, A2, ..., nA,
#'                                                       H1, H2, ... ,nH,Z1, Z2, ... nZ, 
#'                                                       W, R, D, cuminc, cumHospAdm)
#' @param parms List. Model parameters.
#' 
seir <- function(t, x, parms)
{
    with(
        c(as.list(parms), as.list(x)), 
        {
        
        # define compartment parameters in vector x
        S  = x[1]
        Vw = x[2]
        V  = x[3]
        E  = x[4:(nE+3)]
        Ev = x[(nE+4):(nEv+nE+3)]
        I  = x[(nEv+nE+4):(nEv+nE+nI+3)]
        IH = x[(nEv+nE+nI+4):(nEv+nE+nI+nIH+3)]
        A  = x[(nEv+nE+nI+nIH+4):(nEv+nE+nI+nIH+nA+3)]
        H  = x[(nEv+nE+nI+nIH+nA+4):(nEv+nE+nI+nIH+nA+nH+3)]
        Z  = x[(nEv+nE+nI+nIH+nA+nH+4):(nEv+nE+nI+nIH+nA+nH+nZ+3)]
        R  = x[(nEv+nE+nI+nIH+nA+nH+nZ+4)] 
        D  = x[(nEv+nE+nI+nIH+nA+nH+nZ+5)]
        cuminc     = x[(nEv+nE+nI+nIH+nA+nH+nZ+6)]
        cumincsymp = x[(nEv+nE+nI+nIH+nA+nH+nZ+7)]
        cumHospAdm = x[(nEv+nE+nI+nIH+nA+nH+nZ+8)]
        
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
        
        # Apply time-dependent transmission rate
        mult   = broken_line(x=t, b = transm.t, v = transm.v)
        beta_t = beta * mult
        
        # Apply time-dependent vaccine rate
        mult.rt   = broken_line(x=t, b=vacc.rate.t, v=vacc.rate.v)
        rt = r * mult.rt
        
        # calculate incidence
        infrate  = beta_t * S * (rel.inf.a  *sum.A + sum.I + sum.IH) / popSize
        
        ## susceptible
        dS = tau.immu*(R+V) - rt*S - infrate
        
        ## vaccinated but not immuned
        dVw = rt*S - d * Vw
        
        ## vaccinated with full protection 
        vac.infrate = (1-eff.inf) * beta_t * V * (rel.inf.a  *sum.A + sum.I + sum.IH) / popSize
        dV = d*Vw - vac.infrate - tau.immu*V
        
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
        
        # Time-dependent hospital rate
        mmult = step_line(x=t, b = hosp.rate.t, v = hosp.rate.v)
        h_t = h * mmult
        
        # Vaccinated Time-dependent hospital rate
        mmult.vac = step_line(x = t, 
                              b = hosp.rate.vacc.t, 
                              v = hosp.rate.vacc.v)
        h.vac_t = h.vac * mmult.vac
        
        # Time dependent asymptomatic proportion
        mult.alpha = broken_line(x=t, 
                                 b = asymp.prop.t, 
                                 v = asymp.prop.v)
        alpha.t    = min(1.0, mult.alpha * alpha)
        
        # Vaccinated Time dependent asymptomatic proportion
        mult.alpha.vac = broken_line(x = t, 
                                     b = asymp.prop.vacc.t,
                                     v = asymp.prop.vacc.v)
        alpha.vac.t    = min(1.0, mult.alpha.vac * alpha.vac)
        
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
        dR = neta * Z[nZ] + (1-delta) *  nell * H[nH] - tau.immu*R
        
        ## death
        dD = delta * nell * H[nH]
        
        ## Cumulative infection (symp+sympHosp+asymp) incidence
        dcuminc     = infrate + vac.infrate
        dcumincsymp = infrate * (1-alpha.t) + vac.infrate * (1-alpha.vac.t)
        
        ## IMPORTANT: all the derivatives must be 
        ## saved into a big list which is as long as "x"
        res=c(dS, dVw, dV, dE, dEv, dI, dIH, dA, dH, dZ, dR, dD, 
              dcuminc, dcumincsymp, dcumHospAdm)
        list(res)
    })
}
