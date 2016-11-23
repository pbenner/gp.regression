

approximate.posterior.irls.ldB2_exact <- function(W, K, n){
    isWneg <- any(W<0)
    if (isWneg){ # switch between Cholesky and LU decomposition mode
        A <- sweep(K,2,as.matrix(W,n,1),"*") + diag(n) 
        # Multiply W against K row by row elementwise, and add an identity matrix
        lu_matrices <- expand(lu(A))# LU decomposition, A = P*L*U
        diagonal_of_U <- diag(lu_matrices$U)
        sign_of_U  <- prod(sign(diagonal_of_U))
        # sign_of_U is 1 or -1 depending on the number of <0 elements in diagonal_of_U
        if(sign_of_U != det(lu_matrices$P)){ #det(P) is 1 or -1
            ldB2 <- Inf # log becomes complex for negative values, encoded by inf
        } else {
            ldB2 <- sum(log(abs(diagonal_of_U)))/2 #FIXME ldB2 is not necesary!
        }
        Q <- solve( lu_matrices$U, solve(lu_matrices$L,lu_matrices$P) )
        #implement somewhere: solveKiW = @(r) bsxfun(@times,W,Q*r)
    }
    else {
        rootW <- sqrt(W)
        L <- chol(diag(n) + rootW %*% t(rootW) * K)
        ldB2 <- sum(log(diag(L))) 
        Q = FALSE # Q is not necesary for solveKiW when cholsky decomposition is used
    }
    result <- list(ldB2=ldB2, Q=Q)
}
