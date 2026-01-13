### ======================================================================== ###
### constant in the likelihood function for t copula
### c(df, N)/T defined on page 7
### ======================================================================== ###
source("00_LL_utils_func.R")
### ======================================================================== ###


### ======================================================================== ###
### define elimination/permutation matrix needed
### constant matrices only depend on N 
### ======================================================================== ###
matrix_lst <- define_matrix(N)
elimin_diag <- matrix_lst[[1]]
elimin_S <- matrix_lst[[2]]
P_Permut <- matrix_lst[[3]]
### ======================================================================== ###


### ======================================================================== ###
### ep_t: residuals at time t from volatility step 
### theta: correlation matrix at time t
### ======================================================================== ###
score_corr <- function(ep_t, theta, df) {
    N <- ncol(theta)
    eiv <- eigen(logm(theta))
    eigen_vectors <- eiv$vectors   # each column: eigenvector
    eigen_values <- eiv$values
    # sum(abs(expm(eigen_vectors%*%diag(eigen_values)%*%t(eigen_vectors))-theta))
    
    ### elements for matrix At
    At_P1 <- kronecker(eigen_vectors, eigen_vectors, FUN="*")
    # eigen_vectors %x% eigen_vectors
    
    At_P2 <- matrix(nrow=N, ncol=N)
    for(i in 1:N){
        At_P2[i, ] <- (exp(eigen_values[i])-exp(eigen_values))/(eigen_values[i]-eigen_values)
        At_P2[i, i] <- exp(eigen_values[i])
    }   ### symmetric
    At_P2 <- diag(as.vector(vec(At_P2)))
    At <- At_P1%*%At_P2%*%t(At_P1)
    
    score_P1 <- -solve(elimin_diag%*%At%*%t(elimin_diag))%*%(elimin_diag%*%At%*%t(elimin_S))
    Mt <- t(rbind(score_P1, diag(N*(N-1)/2)))%*%t(P_Permut)%*%At
    
    theta_inv <- solve(theta)
    score_P2 <- vec(c(((df+N)/(df-2+t(ep_t)%*%theta_inv%*%ep_t)))*
                        theta_inv%*%ep_t%*%t(ep_t)%*%theta_inv-theta_inv)
    
    ut <- 1/2*Mt%*%score_P2
    return(ut)
}


### ======================================================================== ###
### re-construct the correlation matrix using Hansen algorithm
### ======================================================================== ###
recon_corr_mat <- function(log_theta_vecl, N){
    log_theta_recon <- diag(rep(0, N))
    log_theta_recon[lower.tri(log_theta_recon)] <- log_theta_vecl
    log_theta_recon <- log_theta_recon + t(log_theta_recon)
    
    for(i in 1:1000){
        theta_recon <- expm(log_theta_recon)
        diag(log_theta_recon) <- diag(log_theta_recon) - log(diag(theta_recon))
        
        if(sum(abs(diag(theta_recon)-rep(1,N)))<1e-10){
            break
        }
    }
    return(theta_recon)
}
### check algorithm
# Q <- cor(residual_vol)
# log_theta_vecl <- elimin_lower%*%vec(logm(Q))
# recon_corr_mat(log_theta_vecl, N); Q
### ======================================================================== ###


### ======================================================================== ###
### log-likelihood for correlation matrix estimation
### eps: residuals from the volatility step
### ======================================================================== ###
LL_log_corr_mat <-  function(params_corr, eps, Theta_Cal=FALSE){
    n <- nrow(eps); N <- ncol(eps); NLogCorr <- N*(N-1)/2 ### number of log-corr 
    ### ==================================================================== ###
    ### omega:vector; kappa:constant; phi:constant
    ### ==================================================================== ###
    omega_corr <- matrix(params_corr[1:NLogCorr], nrow=NLogCorr)
    kappa_corr <- matrix(0, nrow=NLogCorr, ncol=NLogCorr)
    diag(kappa_corr) <- rep(params_corr[length(params_corr)-2], NLogCorr)
    phi_corr <- matrix(0, nrow=NLogCorr, ncol=NLogCorr)
    diag(phi_corr) <- rep(params_corr[length(params_corr)-1], NLogCorr)
    df <- 1/params_corr[length(params_corr)]
    ### ==================================================================== ###
    ### log-parameters: vectorize the lower triangular part of the log-correlation matrix
    log_corr_vecl <- matrix(solve(diag(NLogCorr)-phi_corr)%*%omega_corr,
                            nrow=NLogCorr, ncol=n)
    corr_mat_lst <- list()   ### for t=1:T(n)
    corr_mat_lst[[1]] <- recon_corr_mat(log_corr_vecl[, 1], N)
    ut_mat <- matrix(0, nrow=NLogCorr, ncol=n-1)
    
    LL <- tryCatch({
        for(t in 1:(n-1)) {
            ut_mat[, t] <- score_corr(matrix(eps[t, ]), corr_mat_lst[[t]], df)
            
            log_corr_vecl[, t+1] <- (omega_corr + 
                                         phi_corr%*%log_corr_vecl[, t, drop=FALSE] +
                                         kappa_corr%*%ut_mat[, t, drop=FALSE])
            
            ### reconstruct theta matrix [already have non-diag, fill diag]
            corr_mat_lst[[t+1]] <- recon_corr_mat(log_corr_vecl[, t+1], N)
        }
        
        LL_const <- (LL_constant(df, N)-N*LL_constant(df, 1))*n
        LL_vol <- sum((df+1)*log(1+eps^2/(df-2)))
        # LL_t <- sapply(1:n, function(t){
        #     inv_theta <- solve(corr_mat_lst[[t]])
        #     ep_t <- matrix(eps[t, ], ncol=1)
        #     return(log(det(corr_mat_lst[[t]]))+(df+N)*log(1+1/(df-2)*t(ep_t)%*%inv_theta%*%ep_t))
        # })
        LL_t <- sum(sapply(1:n, function(t){
            ep_t <- matrix(eps[t, ], ncol=1)
            chol_theta <- chol(corr_mat_lst[[t]])
            inv_chol_theta <- solve(chol_theta)   ### t(inv()) -- inv(t())
            return(sum(log(diag(chol_theta)))*2+
                       (df+N)*log(1+1/(df-2)*t(ep_t)%*%(inv_chol_theta%*%t(inv_chol_theta))%*%ep_t))
        }))
        LL_C <- -1/2*(LL_t-LL_vol)
        LL <- -(LL_const+LL_C)
    }, error=function(x){return(10^(20)*runif(1,1,2))})
    
    print(LL)
    
    if(Theta_Cal) return(corr_mat_lst)
    
    return(LL)
}
### ======================================================================== ###