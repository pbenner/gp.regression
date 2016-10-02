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

new.likelihood.student_t <- function(alpha, ...) {
    result <- list(alpha = alpha)
    class(result) <- c("likelihood.student_t", "likelihood")
    result
}

logp.likelihood.student_t <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean)
    }
    if (!is.vector(y)) {
        y <- as.vector(y)
    }

   #implement here. 
}

gradient.likelihood.student_t <- function(likelihood, link, f, yp, n) {
    # d: d/dx log p(y|f)
    d     <- as.matrix(rep(0, n))
    # parameter of the student_t likelihood
    alpha <- likelihood$alpha#change here
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
        d[[i]] <- -alpha*Nfx/Pfx*(1-yx/Pfx)#implement here
    }
    return (d)
}

#' Compute the Hessian of a likelihood model
#' 
#' @param model probabilistic model
#' @param ... arguments to be passed to methods

hessian.likelihood.student_t <- function(likelihood, link, f, yp, n) {
    # W: Hessian of log p(y|f)
    W     <- diag(n)
    # parameter of the student_t likelihood
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
                    alpha*Nfx2  /Pfx  *(1 -   yx/Pfx)#implement here
    }
    return (W)
}
