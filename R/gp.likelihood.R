# Copyright (C) 2015 Philipp Benner
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

#' Create a new likelihood model
#' 
#' @param type name of the likelihood model
#' @param ... unused
#' @export

new.likelihood <- function(type, ...)
{
    if (type == "normal" || type == "gaussian") {
        new.likelihood.normal(...)
    }
    else if (type == "gamma") {
        new.likelihood.gamma(...)
    }
    else if (type == "student_t" || type == "t") {
        new.likelihood.student_t(...)
    }
    else {
        stop("Unknown type.")
    }
}

new.likelihood.normal <- function(variance, ...) {
    result <- list(variance = variance)
    class(result) <- c("likelihood.normal", "likelihood")
    result
}

new.likelihood.gamma <- function(alpha, ...) {
    result <- list(alpha = alpha)
    class(result) <- c("likelihood.gamma", "likelihood")
    result
}

#' Evaluate the log probability or log density of the likelihood model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods
#' @export

logp <- function(model, ...)
{
    UseMethod("logp")
}

#' Evaluate the log probability of a probit model
#' 
#' @param model probabilistic model
#' @param y where to evaluate the density
#' @param mean the mean of the distribution
#' @param ... arguments to be passed to methods
#' @method logp NULL
#' @export

logp.NULL <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean)
    }
    stopifnot(length(mean) == nrow(y))

    result <- 0.0

    for (i in 1:nrow(y)) {
        result <- result + y[[i,1]]*log(mean[i]) + y[[i,2]]*log(1.0-mean[i])
    }
    result
}

#' Evaluate the log density of the gamma distribution
#' 
#' @param model probabilistic model
#' @param y where to evaluate the density
#' @param mean the mean of the distribution (i.e. mean = shape*scale)
#' @param ... arguments to be passed to methods
#' @method logp likelihood.gamma
#' @export

logp.likelihood.gamma <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean)
    }
    if (!is.vector(y)) {
        y <- as.vector(y)
    }
    stopifnot(length(mean) == 1 || length(mean) == length(y))

    shape = model$alpha
    scale = mean/shape

    sum(dgamma(y, shape=shape, scale=scale, log=TRUE))
}

#' Compute the gradient of a likelihood model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods

gradient <- function(model, ...)
{
    UseMethod("gradient")
}

gradient.NULL <- function(likelihood, link, f, yp, n) {
    # d: d/dx log p(y|f)
    d <- as.matrix(rep(0, n))
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # counts at position x[[i]]
        c1  <- yp[[i,1]]
        c2  <- yp[[i,2]]
        # value of the response derivative evaluated at fx
        Nfx <- link$response.derivative(fx)
        # response evaluated at fx
        Pfx <- link$response(fx)
        # gradient
        d[[i]] <- c1*Nfx/(0+Pfx) - c2*Nfx/(1-Pfx)
    }
    return (d)
}

gradient.likelihood.gamma <- function(likelihood, link, f, yp, n) {
    # d: d/dx log p(y|f)
    d     <- as.matrix(rep(0, n))
    # parameter of the gamma likelihood
    alpha <- likelihood$alpha
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        yx  <- yp[[i]]
        # value of the response derivative evaluated at fx
        Nfx <- link$response.derivative(fx)
        # response evaluated at fx
        Pfx <- link$response(fx)
        # gradient
        d[[i]] <- -alpha*Nfx/Pfx*(1-yx/Pfx)
    }
    return (d)
}

#' Compute the Hessian of a likelihood model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods

hessian <- function(model, ...)
{
    UseMethod("hessian")
}

hessian.NULL <- function(likelihood, link, f, yp, n) {
    # W: Hessian of log p(y|f)
    W <- diag(n)
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # counts at position x[[i]]
        c1  <- yp[[i,1]]
        c2  <- yp[[i,2]]
        # value of the response derivative evaluated at fx
        Nfx <- link$response.derivative(fx)
        # response evaluated at fx
        Pfx <- link$response(fx)
        # Hessian
        W[[i,i]] <- c1*(-fx*Nfx/(0+Pfx) - Nfx^2/(0+Pfx)^2) -
                    c2*(-fx*Nfx/(1-Pfx) + Nfx^2/(1-Pfx)^2)
    }
    return (W)
}

hessian.likelihood.gamma <- function(likelihood, link, f, yp, n) {
    # W: Hessian of log p(y|f)
    W     <- diag(n)
    # parameter of the gamma likelihood
    alpha <- likelihood$alpha
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        yx <- yp[[i]]
        # value of the first and second response derivative evaluated at fx
        Nfx1 <- link$response.derivative(fx, 1)
        Nfx2 <- link$response.derivative(fx, 2)
        # response evaluated at fx
        Pfx <- link$response(fx)
        # Hessian
        W[[i,i]] <- alpha*Nfx1^2/Pfx^2*(1 - 2*yx/Pfx) -
                    alpha*Nfx2  /Pfx  *(1 -   yx/Pfx)
    }
    return (W)
}
