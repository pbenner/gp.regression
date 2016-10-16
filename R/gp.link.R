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

#' Create a new link function object
#' 
#' @param type name of the link function
#' @param ... unused
#' @export

new.link <- function(type = "probit", ...)
{
    if (type == "probit") {
        new.link.probit(...)
    }
    else if (type == "logistic") {
        new.link.logistic(...)
    }
    else if (type == "null") { # A null link function doesn't do anything. 
    #For likelihood models that doesn't require link functions.
        new.link.null(...)
    }
    else {
        stop("Unknown type.")
    }
}

new.link.probit <- function(...) {
    # qnorm: quantile function  (qnorm = pnorm^-1)
    # pnorm: cumulative density (pnorm = Int dnorm)
    # dnorm: density function   (dnorm = d/dx pnorm)
    result   <- list(link                = qnorm, # link function
                     response            = pnorm, # inverse link function (response)
                     response.derivative = dnorm)
    class(result) <- c("link.probit", "link")
    result
}

new.link.logistic <- function(...) {
    dr <- function(x, n=1) {
             if (n == 1) exp(x)/(1+exp(x))
        else if (n == 2) exp(x)/(1+exp(x)) - exp(2*x)/(1+exp(x))^2
        else stop("Invalid derivative")
    }
    result   <- list(link                = function(x) log(exp(x) - 1),
                     response            = function(x) log(exp(x) + 1),
                     response.derivative = dr)
    class(result) <- c("link.logistic", "link")
    result
}

new.link.null <- function(...) {
    # A null link function doesn't do anything. 
    # It can be used in likelihood models that doesn't require link functions for Laplace approximation.
    result   <- list(link                = function(x) x, # link function
                     response            = NULL,
                     response.derivative = NULL,
    class(result) <- c("link.null", "link")
    result
}


#' Summarize the posterior of a Gaussian process equipped with a probit link function
#' 
#' @param model probit link object
#' @param p.mean mean of the posterior Laplace approximation
#' @param p.variance variance of the posterior Laplace approximation
#' @param ... unused
#' @method summarize link.probit
#' @export

summarize.link.probit <- function(model, p.mean, p.variance, ...)
{
    n        <- length(p.mean)
    mean     <- rep(0, n)
    variance <- rep(0, n)

    for (i in 1:n) {
        # prepare the covariance matrix
        sigma <- matrix(c(1+p.variance[i], p.variance[i], p.variance[i], 1+p.variance[i]), 2, 2)
        # predictive mean
        # (cf. Rasmussen 2006, Eq. 3.82)
        mean[i]     <- pnorm(p.mean[i], mean=0, sd=sqrt(1+p.variance[i]))
        # and predictive variance
        # (cf. math/probit.nb)
        variance[i] <- pmvnorm(upper=c(p.mean[i], p.mean[i]), mean=c(0, 0), sigma=sigma) - mean[i]^2
    }
    return (list(mean = mean, variance = variance))
}

#' Summarize the posterior of a Gaussian process equipped with a logistic link function
#' 
#' There is no analytic solution for the expectation and variance, hence return only
#' the MAP estimate.
#' 
#' @param model logistic link object
#' @param p.mean mean of the posterior Laplace approximation
#' @param p.variance variance of the posterior Laplace approximation
#' @param ... unused
#' @method summarize link.logistic
#' @export

summarize.link.logistic <- function(model, p.mean, p.variance, ...)
{
    mean = model$response(p.mean)
    confidence.1 = model$response(p.mean - 2.0*sqrt(p.variance))
    confidence.2 = model$response(p.mean + 2.0*sqrt(p.variance))
    variance = (confidence.1-confidence.2)^2/4^2

    return (list(mean = mean, variance = variance))
}
# ------------------------------------------------------------------------------
if (FALSE) {

    # binomial observations
    # --------------------------------------------------------------------------
    xp <- c(1,2,3,4)
    yp <- matrix(0, 4, 2)
    yp[1,] <- c(2, 14)
    yp[2,] <- c(4, 12)
    yp[3,] <- c(7, 10)
    yp[4,] <- c(15, 8)

    gp <- new.gp(0.5, kernel.squared.exponential(1, 0.25),
                 likelihood=NULL,
                 link=new.link("probit"))
    gp <- posterior(gp, xp, yp)
    summarize(gp, 1:100/20)

    plot(gp, 1:100/20)

    # gamma distributed observations
    # --------------------------------------------------------------------------
    n  <- 1000
    xp <- 10*runif(n)
    yp <- rgamma(n, 1, 2)

    gp <- new.gp(1.0, kernel.squared.exponential(1.0, 5.0),
                 likelihood=new.likelihood("gamma", 1.0),
                 link=new.link("logistic"))
    # add some tiny noise to the diagonal for numerical stability
    gp <- posterior(gp, xp, yp, ep=0.01, verbose=TRUE)
    summarize(gp, 0:10/5)

    plot(gp, 1:100/10)

}
