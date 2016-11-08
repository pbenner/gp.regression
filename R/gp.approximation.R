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
        repeat { #REFACTOR: make newton algorythm a separate method.
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

approximate.posterior.irls <- function(alpha, mean, K, link){
# Numerically stable mode finding. Influenced by GPML 4.0 (BSD)
# not suer if I need mean here. can I remove K?
#parameters thac can be optional. Put somewhere else later.
    maxit = 20
    Wmin = 0.0
    tol = 1e-6
    smin_line = 0 
    smax_line = 2           # min/max line search steps size range
    nmax_line = 10          # maximum number of line search steps
    thr_line = 1e-4           

    K <- gp$kernelf(gp$xp)
    
}

approximate.posterior.psi  <- function(gp, alpha, K, m){#changing alpha and m works fine.
  # K, yp and hyper parameters also tested.
  # f is calculated so don't need to test, same in matlab.
    alpha = matrix(alpha,length(alpha)) #make sure that alpha is a column vector
    z <- K %*% alpha
    #print(z)
    f <- z + m
    #print('f')
    #print(f)
    #print('f')
    lp = logp(gp$likelihood, gp$yp, f)#I need gp$y here, GPML puts input y into likfun object.
    psi = t(alpha) %*% z /2 - sum(lp)# works with vector and matrix inputs.
    return(psi)
}

approximate.posterior.psi_line <- function(alpha, dalpha, s, K, m, gp){
    psi <- approximate.posterior.psi(gp, alpha + s * dalpha, K, m)
    return(psi)
}#debug this function.

approximate.posterior.search_line <- function(interval=c(0,2), gp, s, alpha, dalpha, m, K){
    #‘f’ will be called as ‘f(x, ...)’ for a numeric value of x.
    alpha <- optimize(approximate.posterior.psi_line,
             interval, dalpha, s, m, K, gp)
    #alpha matches with gpml.
    #f, dlp and W like in line 100 (I think they are updated the same).
    # to do: test search_line
    #        write solveKiW,
    #        test solveKiW
    #        write contents of K.fun,
    #        test K.fun
    #        complete code
    #        test code
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
