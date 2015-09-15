# Copyright (C) 2013 Philipp Benner
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

#' Compute the Kullback-Leibler divergence between two models
#' 
#' @param model0 first distribution
#' @param model1 second distribution
#' @param ... arguments to be passed to methods
#' @export

kl.divergence <- function(model0, model1, ...)
{
    UseMethod("kl.divergence")
}

#' Compute the Kullback-Leibler divergence between two Gaussian
#' processes
#' 
#' @param model0 first Gaussian process
#' @param model1 second Gaussian process
#' @param x positions at which the Gaussian processes are evaluated
#' @param ... unused
#' @method kl.divergence gp
#' @export

kl.divergence.gp <- function(model0, model1, x, ...)
{
    gp0 <- model0
    gp1 <- model1

    list[mean0, variance0, covariance0] <- summarize(gp0, x, return.covariance=TRUE)
    list[mean1, variance1, covariance1] <- summarize(gp1, x, return.covariance=TRUE)

    stopifnot(!is.null(covariance0))
    stopifnot(!is.null(covariance1))

    n      <- length(mean0)

    L0inv  <- solve(chol(covariance0))
    L1inv  <- solve(chol(covariance1))
  
    tmp1   <- log(det((covariance1 %*% L0inv) %*% t(L0inv)))
    tmp2   <- sum(diag(L1inv %*% (t(L1inv) %*% covariance0)))
    tmp3   <- drop(t(mean0 - mean1) %*% (L1inv %*% (t(L1inv) %*% (mean0 - mean1))))

    return (1/2*(tmp1 + tmp2 + tmp3 - n))
}

#' Compute the entropy of a model
#' 
#' @param model for instance a Gaussian process
#' @param ... arguments to be passed to methods
#' @export

entropy <- function(model, ...)
{
    UseMethod("entropy")
}

#' Compute the entropy of a Gaussian process
#' 
#' @param model Gaussian process
#' @param x positions at which the Gaussian process is evaluated
#' @param ... unused
#' @method entropy gp
#' @export

entropy.gp <- function(model, x, ...)
{
    gp <- model
    list[mean, variance, covariance] <- summarize(gp, x, return.covariance=TRUE)
    n <- length(mean)

    return (1/2*(n*log(2*pi*exp(1)) + log(det(covariance))))
}

#' Compute the gain in information between two models
#' 
#' @param model for instance a Gaussian process
#' @param ... arguments to be passed to methods
#' @export

information.gain <- function(model, ...)
{
    UseMethod("information.gain")
}

#' Compute the gain in information between two Gaussian processes
#' 
#' @param model first Gaussian process
#' @param gp2 second Gaussian process
#' @param x positions at which the Gaussian process is evaluated
#' @param ... unused
#' @method information.gain gp
#' @export

information.gain.gp <- function(model, gp2, x1, x2, ep=0.0, ...)
{
    gp1 <- model
    list[mean1, variance1, cov1] <- summarize(gp1, x1, return.covariance=TRUE)
    list[mean2, variance2, cov2] <- summarize(gp2, x2, return.covariance=TRUE)

    diag(cov1) <- diag(cov1) + ep
    diag(cov2) <- diag(cov2) + ep

    Linv <- solve(chol(cov2))

    return (1/2*log(det(cov1 %*% Linv %*% t(Linv))))
}

# ------------------------------------------------------------------------------
if (FALSE) {

    xp <- c(1, 2, 3)
    yp <- c(0.7, 0.7, 0.7)
    ep <- c(0.01, 0.01, 0.01)

    gp0 <- new.gp(0.5, kernel.squared.exponential(1, 1))
    gp1 <- new.gp(0.5, kernel.squared.exponential(1, 2))
    gp2 <- posterior(gp1, xp, yp, ep)

    kl.divergence(gp0, gp0, xp)

    kl.divergence(gp0, gp1, xp)

    kl.divergence(gp1, gp2, xp)

    entropy(gp1, xp)

    entropy(gp2, xp)

}
