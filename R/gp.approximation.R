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

approximate.posterior <- function(gp, epsilon=0.00001, verbose=FALSE, ...)
{
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

approximate.posterior.summary <- function(gp, k1, k2, k3, ...)
{
    if (!is.null(gp$approximation.B.inv.chol)) {
        # Hessian is positive definite
        # (c.f. Rasmussen 2006, Algorithm 3.2)
        mean     <- drop(gp$link$link(gp$prior.mean) + k2 %*% gp$approximation.d)
        v        <- t(gp$approximation.B.inv.chol) %*% (gp$approximation.V %*% k1)
        variance <- diag(k3 - t(v) %*% v)
    }
    else {
        mean     <- drop(gp$link$link(gp$prior.mean) + k2 %*% gp$approximation.d)
        v        <- gp$approximation.V %*% k1
        variance <- diag(k3 - t(v) %*% gp$approximation.B.inv %*% v)
    }
    list(mean = mean, variance = variance)
}
