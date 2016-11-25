# Copyright (C) 2013-2015 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# approximation of the posterior distribution
# ------------------------------------------------------------------------------

approximate.posterior.step <- function(f, yp, mean, L, likelihood, link, n)
{
    # d: d/dx log p(y|f)
    d <- gradient(likelihood, link, f, yp, n)
    # W: negative Hessian of log p(y|f)
    W <- -hessian(likelihood, link, f, yp, n)
    # here we use the matrix inversion lemma (cf. Rasmussen 2006, Eq. 3.27)
    # for (K^-1 + W)^-1 = (K^-1 + W^{1/2} I W^{1/2})^-1,
    # if W has negative elements we need to store the signs in I, i.e. I = sign(W)
    I <- sign(W)
    V <- sqrt(abs(W))
    B <- I + (V %*% t(L)) %*% (L %*% V)
    b <- W %*% (f - mean) + d
    if (all(diag(I) == 1)) {
        # Hessian is positive definite
        K <- solve(chol(B))
        a <- b - V %*% K %*% (t(K) %*% ((V %*% t(L)) %*% (L %*% b)))
    }
    else {
        a <- b - V %*% solve(B) %*% ((V %*% t(L)) %*% (L %*% b))
    }
    f <- mean + t(L) %*% L %*% a

    return (f)
}

approximate.posterior <- function(gp, epsilon=0.00001, verbose=FALSE, method="newton", ...)
{
    print("approximating posterior")
    xp         <- gp$xp
    yp         <- gp$yp
    L          <- gp$prior.sigma.L
    link       <- gp$link
    mean       <- link$link(gp$prior.mean)
    likelihood <- gp$likelihood
    if (is.infinite(mean) || is.nan(mean)) {
        stop("Gaussian process has invalid prior mean!")
    }
    # number of positions where measurements are available
    n <- dim(xp)[[1]]
    # f, fold
    if (!is.null(gp$approximation.d) && nrow(gp$approximation.d) == n) {
        f     <- L %*% (t(L) %*% gp$approximation.d)
        f.old <- f
    }
    else {
        f     <- as.matrix(rep(0, n))
        f.old <- f
    }
    i <- 0
    if (method == "newton"){
        repeat {
            i <- i + 1
            # run Newton steps until convergence
            f <- approximate.posterior.step(f, yp, mean, L, likelihood, link, n)
            if (verbose) {
                print(sprintf("Newton step... %d (error: %f)", i, norm(f - f.old)))
            }
            if (norm(f - f.old) < epsilon) {
                break
            }
            f.old <- f
            }
    } else if (method == "irsl"){
        f <- approximate.posterior.irls() 
    }
    # evaluate the derivative at the current position
    d <- gradient(likelihood, link, f, yp, n)
    W <- -hessian(likelihood, link, f, yp, n)
    I <- sign(W)
    V <- sqrt(abs(W))
    B <- I + (V %*% t(L)) %*% (L %*% V)
    if (all(diag(I) == 1)) {
        # Hessian is positive definite
        gp$approximation.B          <- B
        gp$approximation.B.inv.chol <- solve(chol(B))
        gp$approximation.B.inv      <- NULL
        gp$approximation.d          <- d
        gp$approximation.V          <- V
    }
    else {
        gp$approximation.B          <- B
        gp$approximation.B.inv.chol <- NULL
        gp$approximation.B.inv      <- solve(B)
        gp$approximation.d          <- d
        gp$approximation.V          <- V
    }
    gp
}

#' Numerically stable mode finding. This function and related subroutines are
#' essentialy translations of GPML 4.0 (BSD). See GPML for nomenclature.
#' 
#' @param gp model Gaussian process
#' @param alpha
#' @param mean mean of gp
#' @param K covariance matrix
    #alpha matches with gpml.
    #f, dlp and W like in line 100 (I think they are updated the same).
    # to do: test search_line
    #        write solveKiW,
    #        test solveKiW
    #        write contents of K.fun,
    #        test K.fun
    #        complete code
    #        test code
    #        tidy the code, perhaps use a structure to put all the relevant parameters in.

approximate.posterior.irls <- function(gp, mean, n){
    alpha <- matrix(0, n) #make sure that alpha is a column vector
    maxit <-  20 #settings
    W_vectorMin <- 0.0
    tol <- 1e-6
    K <- gp$kernelf(gp$xp) # I get K = 10 here. it is supposed to be 100 according to gpml.
    f <- K %*% alpha + mean
    d <- gradient(gp$likelihood, gp$link, f, gp$yp, n)
    W_vector <- as.matrix((diag(-hessian(gp$likelihood, gp$link, f, gp$yp, n))))
    Psi_new <- approximate.posterior.irls.psi(gp, alpha, mean, K)
    Psi_old <- Inf
    it = 0
    while(Psi_old - Psi_new > tol && it < 20){# change this to repeat -> break
#     while(1){
        W_vector <- pmax(W_vector,W_vectorMin)
        b <- W_vector * (f - mean) + d
        r = K %*% b# up to here, all the parameters behave the same as gpml
        if(any(W_vector<0)){#whats inside here is making something different from sloveKiW.
            A <- sweep(K,2,as.matrix(W,n,1),"*") + diag(n) 
            # Multiply W_vector against K row by row elementwise, and add an identity matrix
            Q <- solve(A)
            dalpha <- b - sweep(W_vector, 2, Q %*% r, "*") - alpha
        }
        else {
            rootW <- sqrt(W_vector)
            B <- diag(n) + rootW %*% t(rootW) * K # (c.f. Rasmussen 2006, Eq. 3.26)
            L <- chol(B)
            #temp <- solve(L * t(L), sweep(r, 2, rootW, "*"))#this sweep is strange, r and rootW should be both vectors.
            temp <- solve(L, solve(t(L), sweep(r, 1, rootW, "*")))#I get different temp value here!
            # bellow is consistent with gpml!!!
            #dalpha <- b - sweep(temp, 2,  rootW, "*") - alpha
            dalpha <- b - sweep(temp, 1,  rootW, "*") - alpha #sweep(...) gives different result
            #In above, I get dalpha with slightly different value compared to gpml.
                    }
        #update parameters after search
        Psi_old <- Psi_new
        optimisation_step <- approximate.posterior.irls.search_line(gp, alpha, dalpha, mean, K)
        Psi_new = optimisation_step$objective
        alpha <- alpha + dalpha * optimisation_step$minimum
        f <- K %*% alpha + mean
        d <- gradient(gp$likelihood, gp$link, f, gp$yp, n)
        W_vector <- as.matrix(diag(-hessian(gp$likelihood, gp$link, f, gp$yp, n)))
        it = it + 1
    }
    return(alpha)
}

#' Computes psi. The irls search is performed by minimising psi, however
#' psi_line is used as the objective function in actuality. This is because psi takes
#' vectors as inputis rather than a scalar, which makes it difficult to use in a
#' typical optimisation routine.
#' 
#' @param gp model Gaussian process
#' @param alpha
#' @param mean
#' @param K covariance matrix

approximate.posterior.irls.psi  <- function(gp, alpha, mean, K){#changing alpha and mean works fine.
  # K, yp and hyper parameters also tested.
  # f is calculated so don't need to test, same in matlab.
    alpha <- matrix(alpha,length(alpha)) #make sure that alpha is a column vector
    z <- K %*% alpha
    f <- z + mean
    lp <- logp(gp$likelihood, gp$yp, f)
    psi <- t(alpha) %*% z /2 - sum(lp)
    return(psi)
}

#' Computes psi_lin, the objective function of the search_line routine.
#' Given alpha and dalpha as constant vectors, search_line optimises 
#' s, a scalar, in order to minimise psi.
#' 
#' @param alpha vector or a matrix with one column
#' @param dalpha vector or a matrix with one column
#' @param s scalar
#' @param mean 
#' @param K covariance matrix
#' @param gp model Gaussian process

approximate.posterior.irls.psi_line <- function(s, alpha, dalpha, mean, K, gp){
    psi <- approximate.posterior.irls.psi(gp, alpha + s * dalpha, mean, K)
    return(psi)
}


#' Runs optimisation routine using psi_line as the objective function.
#' 
#' @param alpha
#' @param dalpha
#' @param s
#' @param mean
#' @param K covariance matrix
#' @param gp model Gaussian process

approximate.posterior.irls.search_line <- function(gp, alpha, dalpha, mean, K, s_interval=c(0,2)){
    smin_line <- 0 
    smax_line <- 2           # min/max line search steps size range
    nmax_line <- 10          # maximum number of line search steps
    thr_line <- 1e-4           
    result <- optimize(approximate.posterior.irls.psi_line,
             s_interval, alpha, dalpha, mean, K, gp) # DON'T optimise for alpha, its a vector. Optimise for s.
}

approximate.posterior.summary <- function(gp, k1, k2, k3, ...)
{
    if (!is.null(gp$approximation.B.inv.chol)) {
        # Hessian is positive definite
        # (c.f. Rasmussen 2006, Algorithm 3.2)
        mean     <- drop(gp$link$link(gp$prior.mean) + k2 %*% gp$approximation.d)
        v        <- t(gp$approximation.B.inv.chol) %*% (gp$approximation.V %*% k1)
        variance <- diag(k3 - t(v) %*% v)#cf. Rasmussen 2006, Eq. 3.29
    }
    else {
        mean     <- drop(gp$link$link(gp$prior.mean) + k2 %*% gp$approximation.d)
        v        <- gp$approximation.V %*% k1
        variance <- diag(k3 - t(v) %*% gp$approximation.B.inv %*% v)
    }
    list(mean = mean, variance = variance)
}
