### ======================================================================== ###
### constant in the likelihood function 
### c(df, N)/T defined on page 7
### ======================================================================== ###
LL_constant <- function(df, N) {
    return(log(gamma((df+N)/2))-log(gamma(df/2))-N/2*log((df-2)*pi))
}
### ======================================================================== ###


### ======================================================================== ###
### define related elimination and permutation matrices for symmetric matrix Q 
### ======================================================================== ###
define_matrix <- function(N) {
    elimin_matrix <- elimination.matrix(N)
    C <- matrix(0, N, N)
    C[lower.tri(C, diag = TRUE)] <- 1:nrow(elimin_matrix)
    
    elimin_diag <- elimin_matrix[diag(C),]
    elimin_lower <- elimin_matrix[C[lower.tri(C)], , drop=F]
    
    tmp <- as.vector(t(matrix(1:(N^2), nrow=N, ncol=N)))
    K <- sparseMatrix(i=1:(N^2), j=tmp, x=1, dims=c(N^2, N^2))
    elimin_upper <- as.matrix(elimin_lower %*% K)
    
    elimin_S <- elimin_lower + elimin_upper
    
    D <- matrix(0, N, N)
    D[row(D)==col(D)] <- 1:dim(D)[1]
    D[lower.tri(D)] <- (dim(D)[1]+1):(dim(D)[1]*(dim(D)[1]+1)/2)
    D[upper.tri(D)] = t(D)[upper.tri((D))]
    P_Permut <- diag(dim(D)[1]*(dim(D)[1]+1)/2)[vec(D),]
    
    ### ==================================================================== ###
    ### check for defined matrix ###
    ### ==================================================================== ###
    # Q <- matrix(1:(N^2), N, N)
    # elimin_lower%*%vec(Q)
    # elimin_diag%*%vec(Q)
    # elimin_upper%*%vec(Q)
    
    # Q[upper.tri(Q)] <- t(Q)[upper.tri((Q))]
    # P_Permut%*%matrix(c(diag(Q), t(elimin_lower%*%vec(Q))), ncol=1)-vec(Q)
    
    return(list(elimin_diag, elimin_S, P_Permut))
}
### ======================================================================== ###
