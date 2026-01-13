### ======================================================================== ###
### constant in the likelihood function for student-t distribution
### c(df, N)/T defined on page 7
### ======================================================================== ###
source("00_LL_utils_func.R")
### ======================================================================== ###


### ======================================================================== ###
### Log-likelihood functions for Beta-t-EGARCH model
### params=[Omega, Kappa, Phi, Eta_bar]
### Eta_bar = 1/eta [0,1]
### log_vol: h_t = omega+phi*h_(t-1)+kappa*ut
### ======================================================================== ###
LL_vol_BetaTEGARCH <-  function(params, r, vol_cal=FALSE){
    n <- length(r)
    omega <- params[1]; kappa <- params[2]; phi <- params[3]
    df <- 1/params[4]   ### eta is eta_bar=1/d.o.f
    
    # log_vol <- matrix(log(sd(r)), n, 1)
    log_vol <- matrix(omega/(1-phi), n, 1)
    ut <- matrix(0, n-1, 1)
    
    for (t in 1:(n-1)){
        ut[t] <- (df+1)*r[t]^2/((df-2)*exp(2*log_vol[t])+r[t]^2)-1
        log_vol[t+1] <- omega+phi*log_vol[t]+kappa*ut[t]
    }
    
    LL_t <- (df+1)/2*log(1+r^2*exp(-2*log_vol)/(df-2))
    LL_vec <- LL_constant(df, 1)-log_vol-LL_t
    
    if(vol_cal) return(log_vol)
    
    return(-sum(LL_vec))
}
### ======================================================================== ###


