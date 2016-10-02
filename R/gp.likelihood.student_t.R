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

#' Create a new likelihood model
#' 
#' @param type name of the likelihood model
#' @param ... unused
#' @export

new.likelihood.gamma <- function(alpha, ...) {
    result <- list(alpha = alpha)
    class(result) <- c("likelihood.gamma", "likelihood")#classname likelihood.gamma
    result
}

logp.likelihood.gamma <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean)
    }
    if (!is.vector(y)) {
        y <- as.vector(y)
    }
    stopifnot(length(mean) == 1 || length(mean) == length(y))#length of mean has to be 1 or the same as length of y, but why?

    shape = model$alpha
    scale = mean/shape# scale (1/beta) is determined by alpha and mean, therefore not necessary as an input paramter.

    sum(dgamma(y, shape=shape, scale=scale, log=TRUE))#This is just plane old gamma distribution, but why the sum? approximation of integration? maybe see SUMMATION of wikipedia page of Gamma distribution.
}

gradient.likelihood.gamma <- function(likelihood, link, f, yp, n) {#What is f and yp??? yp must be y datapoints as far as the other codes go. f must be the posterior step.
    # d: d/dx log p(y|f)
    d     <- as.matrix(rep(0, n))
    # parameter of the gamma likelihood
    alpha <- likelihood$alpha
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]] #and yes, observation at x[[i]] is yp[[i]]
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
