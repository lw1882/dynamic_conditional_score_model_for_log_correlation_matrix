### ======================================================================== ###
source("main_00_init.R")
### ======================================================================== ###
load("output/residual_vol.RData"); N <- ncol(residual_vol)
### ======================================================================== ###
source("02_log_corr_matrix_t_func.R")
### ======================================================================== ###


### ======================================================================== ###
### params = [omega, kappa, phi, eta]; omega is vector; kappa, phi, eta are constants
omega_param <- rep(0.01, N*(N-1)/2)*runif(N*(N-1)/2, 0.8, 1)
kappa_param <- 0.01
phi_param <- 0.95
eta_param <- 0.12   ### 1/df, i.e. eta_bar
params_corr_ini <- c(omega_param, kappa_param, phi_param, eta_param)
### ======================================================================== ###
LL_log_corr_mat(params_corr_ini, residual_vol)   # -175.3277
### ======================================================================== ###
res_corr <- solnp(pars=params_corr_ini, fun=LL_log_corr_mat,
                  LB=c(rep(-10, N*(N-1)/2), 0, 0, 0),
                  UB=c(rep(10, N*(N-1)/2), 1, 1, 0.5),
                  control=list(delta=1e-6,tol=1e-12),
                  eps=residual_vol)
sqrt(diag(solve(res_corr$hessian)))
res_corr$pars/sqrt(diag(solve(res_corr$hessian)))
### ======================================================================== ###


# ### ======================================================================== ###
# ### compute correlation matrix time series
# ### ======================================================================== ###
corr_mat_Lst <- LL_log_corr_mat(res_corr$pars, eps=residual_vol,
                                Theta_Cal=TRUE)
corr_ts <- t(sapply(1:length(corr_mat_Lst), function(t) {
    # elimin_lower%*%vec(corr_mat_Lst[[t]])
    corr <- corr_mat_Lst[[t]]
    corr[lower.tri(corr, diag = F)]
}))
plot(xts(corr_ts, order.by=(index(residual_vol))))
# ### ======================================================================== ###
save(corr_mat_Lst, file="output/corr_mat_Lst.RData")
# ### ======================================================================== ###

