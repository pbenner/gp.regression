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

#' @name gp.regression
#' @docType package
#' @title Gaussian Process Regression
#' @author Philipp Benner \email{Philipp.Benner@@mis.mpg.de}
#' @import ggplot2
#' @import scales
#' @importFrom Matrix chol
#' @importFrom mvtnorm pmvnorm
#' @importFrom gridExtra grid.arrange
#' @useDynLib gp.regression
NULL

library("Matrix")

#' Generate a new Gaussian process
#' 
#' @param mean prior mean
#' @param kernelf kernel function
#' @param dim dimension of the x-values
#' @param likelihood a likelihood model of not Gaussian
#' @param link a link function used to transform the target domain
#' @export

new.gp <- function(mean, kernelf, dim=1, likelihood=new.likelihood("normal", 1.0), link=NULL)
{
    stopifnot(length(mean) == 1)

    if (any(class(likelihood) == "likelihood.normal") && !is.null(link)) {
        stop("The Gaussian likelihood should not be combined with a link function!")
    }
    if (any(class(link) == "link.probit") && !is.null(likelihood)) {
        stop("The probit link function should not be compined with a particular likelihood model.")
    }
    if (any(class(likelihood) == "likelihood.gamma") && is.null(link)) {
        stop("A link function is required for the gamma likelihood model!")
    }
    if (any(class(likelihood) == "likelihood.student_t") && is.null(link)) {
       link = new.link.null()
    }

    gp <- list(xp               = NULL,       # data
               yp               = NULL,
               dim              = dim,        # dimension of the x-values
               kernelf          = kernelf,    # kernel function
               likelihood       = likelihood, # likelihood model (requires a link function)
               link             = link,       # link function
               prior.mean       = mean,       # mean
               prior.sigma.L    = NULL,       # cholesky decomposition of covariance matrix
               prior.sigma.Linv = NULL        # inverse cholesky decomposition
               )
    class(gp) <- "gp"

    return (gp)
}

#' Show the dimension of the GP
#' 
#' @param x Gaussian process
#' @method dim gp
#' @export

dim.gp <- function(x) {
    x$dim
}

#' Summarize a distribution
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

summarize <- function(model, ...)
{
    UseMethod("summarize")
}

#' Compute posterior expectation and variance of a GP
#' 
#' @param model object
#' @param x where to evaluate the posterior
#' @param return.covariance if true, return the covariance matrix as well
#' @param ... unused
#' @method summarize gp
#' @export

summarize.gp <- function(model, x, return.covariance = FALSE, ...)
{
    gp <- model

    if (is.null(x)) {
        return (NULL)
    }
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    # check dimension of the x-values
    stopifnot(dim(x)[2] == gp$dim)
    # xp and yp can be NULL, then simplify compute the prior
    if (is.null(gp$xp) || is.null(gp$yp)) {
        mean       <- rep(gp$prior.mean, dim(x)[1])
        covariance <- gp$kernelf(x)
        variance   <- diag(covariance)
    }
    else {
        k1 <- gp$kernelf(gp$xp, x) # K(X , X*)
        k2 <- t(k1)                # K(X*, X )
        k3 <- gp$kernelf(x)        # K(X*, X*)

        # check if a link is available
        if (!is.null(gp$link)) {
            list[mean, variance] <- approximate.posterior.summary(gp, k1, k2, k3, ...)
            # propagate mean and variance through the link function
            # (c.f. Rasmussen 2006, Eq. 3.25)
            covariance           <- NULL
            if (!class(gp$link)[[1]] == "link.null"){
                list[mean, variance] <- summarize(gp$link, mean, variance, ...)
            }
        }
        else {
            # posterior mean
            mean       <- drop(gp$prior.mean + (k2 %*% gp$prior.sigma.Linv) %*% (t(gp$prior.sigma.Linv) %*% (gp$yp - gp$prior.mean)))
            # posterior covariance
            covariance <- k3 - (k2 %*% gp$prior.sigma.Linv) %*% (t(gp$prior.sigma.Linv) %*% k1)
            variance   <- diag(covariance)
        }
    }
    if (return.covariance) {
        return (list(mean = mean, variance = variance, covariance = covariance))
    }
    else {
        return (list(mean = mean, variance = variance))
    }
}

#' Compute posterior of a model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

posterior <- function(model, ...)
{
    UseMethod("posterior")
}

#' Compute posterior of a Gaussian process
#' 
#' @param model Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param ep uncertainty of measurements (optional)
#' @param ... unused
#' @method posterior gp
#' @export

posterior.gp <- function(model, xp, yp, ep=NULL, ...)
{
    gp <- model

    # check arguments
    if (!is.null(xp) && !is.matrix(xp)) {
        xp <- as.matrix(xp, 1)
    }
    # yp is usually a vector unless a link function is used
    if (!is.null(yp) && !is.matrix(yp)) {
        yp <- as.matrix(yp)
    }
    if (!is.null(ep) && !is.vector(ep)) {
        ep <- as.vector(ep)
    }
    # check dimension of the x-values
    stopifnot(dim(xp)[2] == gp$dim)

    # extend covariance matrix
    if (is.null(gp$xp)) {
        # first time the posterior is computed
        A12 <- NULL
        A22 <- gp$kernelf(xp)
    }
    else {
        A12 <- gp$kernelf(gp$xp, xp)
        A22 <- gp$kernelf(xp)
    }
    # update gp measurements
    gp$xp <- rbind (gp$xp, xp)
    gp$yp <- rbind (gp$yp, yp)
    # add Gaussian noise to measurements?
    if (any(class(gp$likelihood) == "likelihood.normal")) {
        if (is.null(ep)) {
            # use the variance specified in the likelihood model
            A22 <- A22 + diag(gp$likelihood$variance, dim(xp)[1])
        }
        else {
            # ignore the likelihood model and use ep instead
            A22 <- A22 + diag(ep, dim(xp)[1])
        }
    }
    else {
        if (!is.null(ep) && any(ep != 0.0)) {
            warning("Adding Gaussian noise to observations although no Gaussian likelihood is used!")
            # numerical stability
            A22 <- A22 + diag(ep, dim(xp)[1])
        }
    }
    # update cholesky decomposition and its inverse
    tmp <- cholesky.inverse.update(gp$prior.sigma.L, gp$prior.sigma.Linv, A12, A22)#error occures here when t-likelihood is used with ep=NULL.
    gp$prior.sigma.L    <- tmp$K
    gp$prior.sigma.Linv <- tmp$Kinv

    # check if a link is available
    if (!is.null(gp$link)) {
        # if so, then call a specialized method to
        # approximate the posterior of the gaussian process
        gp <- approximate.posterior(gp, ...)
    }

    gp
}

#' Draw samples from a probability distribution
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

draw.sample <- function(model, ...)
{
    UseMethod("draw.sample")
}

#' Draw samples from a Gaussian process
#' 
#' @param model Gaussian process
#' @param x locations where to evaluate the Gaussian process
#' @param ep uncertainty of measurements (optional)
#' @param ... arguments to be passed to the algorithms
#' @method draw.sample gp
#' @export

draw.sample.gp <- function(model, x, ep=NULL, ...)
{
    gp <- model

    x <- as.matrix(x)
    y <- rep(0.0, nrow(x))

    for (i in 1:length(x)) {
        list[mean, variance] <- summarize(gp, x[i,], ...)
        y[i] <- rnorm(1, mean = mean, sd = sqrt(variance))
        gp <- posterior(gp, x[i,], y[i], ep=ep, ...)
    }
    y
}

#' Compute marginal likelihood of a model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

marginal.likelihood <- function(model, ...)
{
    UseMethod("marginal.likelihood")
}

#' Compute marginal likelihood of a Gaussian process
#' 
#' @param model Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param ep uncertainty of measurements (optional)
#' @param ... unused
#' @method marginal.likelihood gp
#' @export

marginal.likelihood.gp <- function(model, xp, yp, ep=NULL, ...)
{
    gp1 <- model
    gp2 <- posterior(gp1, xp, yp, ep, ...)

    ml <- function(gp) {
        if (is.null(gp$xp)) {
            return (0.0)
        }
        # check if a link is available
        if (!is.null(gp$link)) {
            k <- gp$kernelf(gp$xp)
            list[mean, variance] <- approximate.posterior.summary(gp, k, k, k, ...)
            list[mean, variance] <- summarize(gp$link, mean, variance, ...)
            result <- -1/2*drop((mean %*% gp$prior.sigma.Linv) %*% (t(gp$prior.sigma.Linv) %*% mean))
            result <- result + logp(gp$likelihood, gp$yp, mean) - 1/2*log(det(gp$approximation.B))
        }
        else {
            # log marginal likelihood
            result <- -1/2*drop((t(gp$yp - gp$prior.mean) %*% gp$prior.sigma.Linv) %*% (t(gp$prior.sigma.Linv) %*% (gp$yp - gp$prior.mean)))
            result <- result - sum(log(diag(gp$prior.sigma.L))) - dim(gp$xp)[1]/2*log(2*pi)
        }
        result
    }
    ml(gp2) - ml(gp1)
}

#' Compute gradient of the marginal likelihood with respect to model parameters
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

gradient.marginal.likelihood <- function(model, ...)
{
    UseMethod("gradient.marginal.likelihood")
}

#' Compute gradient of the marginal likelihood of a Gaussian process with respect
#' to the kernel parameters
#' 
#' @param model Gaussian process
#' @param n number of parameters of the model
#' @param xp positions of measurements
#' @param yp measured values
#' @param ep uncertainty of measurements (optional)
#' @param ... unused
#' @method gradient.marginal.likelihood gp
#' @export

gradient.marginal.likelihood.gp <- function(model, n, xp, yp, ep=NULL, ...)
{
    gp <- posterior(model, xp, yp, ep, ...)

    # initialize gradient to zero
    gradient <- rep(0.0, n)
    
    alpha <- gp$prior.sigma.Linv %*% (t(gp$prior.sigma.Linv) %*% (gp$yp - gp$prior.mean))
    tmp   <- alpha%*%t(alpha) - gp$prior.sigma.Linv%*%t(gp$prior.sigma.Linv)
    
    # loop over parameters
    for (i in 1:n) {
        k0d <- gp$kernelf(xp, gradient=TRUE, i=i)
        gradient[i] <- 1/2*sum(diag(tmp%*%k0d))
    }
    gradient
}

#' Maximize marginal likelihood with respect to the model parameters
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

maximize.marginal.likelihood <- function(model, ...)
{
    UseMethod("maximize.marginal.likelihood")
}

#' Compute gradient of the marginal likelihood of a Gaussian process with respect
#' to the kernel parameters
#' 
#' @param moden Gaussian process
#' @param get.kernel function (specified as a variable or string) that returns the (programatic) kernel function for the given parameters.
#' @param init initial parameters
#' @param xp positions of measurements
#' @param yp measured values
#' @param ep uncertainty of measurements (optional)
#' @param eta adaptive step size parameter for the RProp algorithm
#' @param epsilon stop criterium
#' @param step.init initial step size
#' @param verbose print norm(gradient)^2 and parameter values at every step
#' @param ... unused
#' @method maximize.marginal.likelihood gp
#' @export

maximize.marginal.likelihood.gp <- function(model, get.kernel, init, xp, yp, ep=NULL,
                                            eta=0.01, epsilon=0.01, step.init = 0.1, verbose=FALSE, ...)
{
    gp <- model
    # conditional marginal likelihood maximization is not supported
    stopifnot(is.null(gp$xp))
    # GP with link function is not supported
    stopifnot(is.null(gp$link))
    # initialize vector of parameters
    pt <- as.vector(init)
    # number of parameters
    n  <- length(pt)
    # initialize adaptive step size
    step.size <- rep(step.init, n)
    # and the old gradient just to one
    gradient.old <- rep(1.0, n)

    repeat {
        # get kernel for the current parameters
        gp$kernelf <- do.call(get.kernel, as.list(pt))
        # compute the gradient at the current position pt
        gradient.new <- gradient.marginal.likelihood(gp, n, xp, yp, ep=ep, ...)

        if (verbose) {
            print(sprintf("1/n*norm(gradient)^2 = %f, parameters: %f",
                           sum(gradient.new^2)/n, pt))
            print(sprintf("marginal likelihood: %f",
                          marginal.likelihood(gp, xp, yp, ep)))
        }
        # check if stop condition is met
        if (sum(gradient.new^2)/n < epsilon) {
            return (gp)
        }
        # update step size
        for (i in 1:n) {
            if (gradient.old[i]*gradient.new[i] > 0.0) {
                step.size[i] = (1.0+eta)*step.size[i]
            }
            else {
                step.size[i] = (1.0-eta)*step.size[i]
            }
        }
        # update parameters
        pt <- pt + sign(gradient.new)*step.size
        # save current gradient
        gradient.old <- gradient.new
    }
}

# ------------------------------------------------------------------------------
if (FALSE) {

    xp <- c(1, 2, 3)
    yp <- c(0.7, 0.7, 0.7)
    ep <- c(0.01, 0.01, 0.01)

    gp1 <- new.gp(0.5, kernel.squared.exponential(1, 1),
                  likelihood = new.likelihood("normal", 0.01))
    gp1 <- posterior(gp1, xp, yp)

    plot(gp1, 0:40/10)

    # update in two steps
    gp2 <- new.gp(0.5, kernel.squared.exponential(1, 1))
    gp2 <- posterior(gp2, xp[  1], yp[  1], ep[  1])
    gp2 <- posterior(gp2, xp[2:3], yp[2:3], ep[2:3])

    # evaluate marginal likelihood
    gp3 <- new.gp(0.5, kernel.squared.exponential(1, 1))
    marginal.likelihood(gp3, xp, yp, ep)
    # [1] -2.292447

    # evaluate conditional marginal likelihood
    gp4 <- new.gp(0.5, kernel.squared.exponential(1, 1))
    gp4 <- posterior(gp4, xp[1], yp[1], ep[1])
    marginal.likelihood(gp3, xp[2:3], yp[2:3], ep[2:3])
    # [1] -1.648935

}
