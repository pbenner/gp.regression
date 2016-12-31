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
#' @param ... udfsed
#' @export

# implementation of non-standardized student's t likelihood.
new.likelihood.weibull <- function(shape, ..) {
    result <- list(shape = shape) 
    class(result) <- c("likelihood.weibull", "likelihood")
    result
}

logp.likelihood.weibull <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean)#figure out how to set mean
    }
    if (!is.vector(y)) {
        y <- as.vector(y)
    }
    shape <- model$shape #I want to keep the shape parameter constant,
    scale <- mean / ( gamma( 1 + 1 / shape) ) #so here I calculate the scale parameter from mean.
    result <- dweibull( y, shape = shape, scale = scale, log = TRUE)
    return (as.matrix(result))
}

gradient.likelihood.weibull <- function(likelihood, link, f, yp, n) {
    # d: d/dx log p(y|f)
    d <- as.matrix(rep(0, n))
    # parameter of the weibull likelihood
    df <- likelihood$ka
    sigma <- likelihood$sigma
    #calculate mu here.
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        yx  <- yp[[i]]
        r <- yx - fx
        rsqwr <- r*r
        a <- rsqwr+df*sigma^2
        # gradient
        d[[i]] <- (df+1)*r/a
    }
    return (d)
}

#' Compute the Hessian of a likelihood model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods

hessian.likelihood.weibull <- function(likelihood, link, f, yp, n, form = "matrix") {
    # W: Hessian of log p(y|f)
    W <- vector(mode = "numeric", length = n)
    # parameter of the weibull likelihood
    ka <- likelihood$ka
    sigma <- likelihood$sigma
    sn2 = sigma^2
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        # Hessian
        W[[i,i]] <- (df+1)*(rsqwr-df*sn2)/a^2;#check df is correctly defined, likely to need +1.
    }
    if (form == "vector") return(W)
    else return(diag(n))
}
