### ======================================================================== ###
source("main_00_init.R")
source("01_vol_beta_t_egarch_func.R")
### ======================================================================== ###


### ======================================================================== ###
assets_returns <- read_csv("data/asset_returns.csv")
assets_returns <- as.xts(assets_returns[,-1], order.by=assets_returns$Index, 
                         format="%Y-%m-%d")
tail(assets_returns)
### ======================================================================== ###
### standardize the assets returns
# assets_returns <- sweep(assets_returns, 2, colMeans(assets_returns), "-")
# assets_returns <- sweep(assets_returns, 2, colSds(assets_returns), "/")
colMeans(assets_returns); colSds(assets_returns)
### ======================================================================== ###
N <- ncol(assets_returns)
### ======================================================================== ###


### ======================================================================== ###
### Parameter estimation [LL] ["BETATEGARCH"]
### assets_returns: each row for one asset [standardized]
### ======================================================================== ###
col_names <- c("Omega", "Kappa", "Phi", "Eta_bar")
params_vol <- matrix(nrow=N, ncol=length(col_names))
colnames(params_vol) <- col_names
rownames(params_vol) <- toupper(colnames(assets_returns))
sderror_vol <- params_vol
for(i in 1:N){
    ### ==================================================================== ###
    params_i <- c(0.0001, 0.01, 0.95, 0.15)
    res <- solnp(pars=params_i, fun=LL_vol_BetaTEGARCH, 
                 LB=c(-10, 0, 0, 0), UB=c(10, 1, 1, 0.5), 
                 r=assets_returns[, i])
    se <- sqrt(diag(solve(res$hessian)))
    
    params_vol[i,] <- res$par   # estimated parameters
    sderror_vol[i,] <- se   # standard errors
}
### ======================================================================== ###
params_vol/(sderror_vol)   # t statistics
1/params_vol[, "Eta_bar"]   # d.o.f
### ======================================================================== ###
assets_log_vol <- sapply(1:N, function(i) {
    LL_vol_BetaTEGARCH(c(t(params_vol[i, ])), assets_returns[, i], vol_cal=TRUE)
}); colnames(assets_log_vol) <- colnames(assets_returns)
assets_vol <- exp(assets_log_vol)
plot(xts(assets_vol, order.by=index(assets_returns)), legend.loc="topleft",
     col=c("darkblue", "darkred", "darkgreen"))
### ======================================================================== ###
residual_vol <- assets_returns/assets_vol
colMeans(residual_vol); colSds(residual_vol)
residual_vol <- sweep(residual_vol, 2, colMeans(residual_vol), "-")
residual_vol <- sweep(residual_vol, 2, colSds(residual_vol), "/")
save(residual_vol, file="output/residual_vol.RData")
### ======================================================================== ###

