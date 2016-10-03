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
new.likelihood.student_t <- function(df, sigma, ..) {
    #df: degree of freedom, sigma: scale parameter.
    result <- list(df = df, sigma = sigma) 
    class(result) <- c("likelihood.student_t", "likelihood")
    result
}

logp.likelihood.student_t <- function(model, y, mean, ...)
{
    if (!is.vector(mean)) {
        mean <- as.vector(mean) #mean defined against what?
    }
    if (!is.vector(y)) {
        y <- as.vector(y)
    }
    df <- model$df
    sigma <- model$sigma
    result <- sigma * df(x, df, 0) #How do a specify x here? mean?
    return (result)
}

gradient.likelihood.student_t <- function(likelihood, link, f, yp, n) {
    # d: d/dx log p(y|f)
    d     <- as.matrix(rep(0, n))
    # parameter of the student_t likelihood
    df <- likelihood$df
    sigma <- model$sigma
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        yx  <- yp[[i]]
        # gradient
        r = student_t.r(y, mu) #copied and pasted 3 lines from gpml
        d[[i]] <- (nu+1)*r./a
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
    df <- likelihood$df
    for (i in 1:n) {
        # current f value at x[[i]]
        fx  <- f[[i]]
        # observation at position x[[i]]
        yx <- yp[[i]]
        r <- student_t.r(y, mu)
        #copied and pasted 3 lines from gpml
        rsqwr <- student_t.rsqwr(r)
        # Hessian
        W[[i,i]] <- (nu+1)*(rsqwr-nu*sn2)./a.^2;
    }
    return (W)
}

student_t.r <- function(y, mu){
    r = y - mu
}

student_t.rsqwr <- function(r){
    rsqwr = r*r
}
