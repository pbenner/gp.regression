# Copyright (C) 2016 Gen Kamita
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
    alpha <- matrix(0, n)
    maxit <-  100 #setting
    W_vectorMin <- 0.0
    tol <- 1e-6
    K <- gp$kernelf(gp$xp) 
    Psi_new <- approximate.posterior.irls.psi(gp, alpha, mean, K)
    f <- K %*% alpha + mean
    it <- 0
    repeat{
        d <- gradient(gp$likelihood, gp$link, f, gp$yp, n)
        W_vector <- as.matrix(-hessian(gp$likelihood, gp$link, f, gp$yp, n, form = 'vector'))
        W_vector <- pmax(W_vector,W_vectorMin)
        b <- W_vector * (f - mean) + d
        r = K %*% b
        if(any(W_vector<0)){
            A <- sweep(K, 2, as.matrix(W_vector, n, 1), "*") + diag(n) 
            # Multiply W_vector against K row by row elementwise, and add an identity matrix
            Q <- solve(A)
            dalpha <- b - sweep(W_vector, 2, Q %*% r, "*") - alpha
        }
        else {
            rootW <- sqrt(W_vector)
            B <- diag(n) + rootW %*% t(rootW) * K # (c.f. Rasmussen 2006, Eq. 3.26)
            L <- chol(B)
            temp <- solve(L, solve(t(L), sweep(r, 1, rootW, "*")))
            dalpha <- b - sweep(temp, 1,  rootW, "*") - alpha 
                    }
        #update parameters after search
        optimisation_step <- approximate.posterior.irls.search_line(gp, alpha, dalpha, mean, K)
        Psi_old <- Psi_new
        Psi_new = optimisation_step$objective
        alpha <- alpha + dalpha * optimisation_step$minimum
        f <- K %*% alpha + mean
        it = it + 1
        if (Psi_old - Psi_new < tol || it > maxit) break
    }
    return(f)
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
             s_interval, alpha, dalpha, mean, K, gp) 
}
