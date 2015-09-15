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

call.kernel <- function(name, x, y, ...)
{
  storage.mode(x)   <- "double"
  if (!is.null(y)) {
    storage.mode(y)   <- "double"
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(y) && !is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (is.null(y)) {
    .Call(name, x, x, ..., PACKAGE="gp.regression")
  }
  else {
    .Call(name, x, y, ..., PACKAGE="gp.regression")
  }
}

#' Create a squared exponential kernel
#' 
#' @param l length scale
#' @param var noise variance
#' @export

kernel.squared.exponential <- function(l, var)
{
    storage.mode(l)   <- "double"
    storage.mode(var) <- "double"
    f <- function(x, y=NULL, gradient=FALSE, i = 0) {
        if (gradient) {
            storage.mode(i) <- "integer"
            call.kernel("squared_exponential_gradient", x, y, l, var, i)
        }
        else {
            call.kernel("squared_exponential_kernel", x, y, l, var)
        }
    }
    return (f)
}

#' Create a gamma exponential kernel
#' 
#' @param l length scale
#' @param var noise variance
#' @param gamma exponent
#' @export

kernel.gamma.exponential <- function(l, var, gamma)
{
    storage.mode(l)     <- "double"
    storage.mode(var)   <- "double"
    storage.mode(gamma) <- "double"
    f <- function(x, y=NULL, gradient=FALSE, i = 0) {
        if (gradient) {
            storage.mode(i) <- "integer"
            call.kernel("gamma_exponential_gradient", x, y, l, var, gamma, i)
        }
        else {
            call.kernel("gamma_exponential_kernel", x, y, l, var, gamma)
        }
    }
    return (f)
}

#' Create an Ornstein-Uhlenbeck kernel
#' 
#' @param l length scale
#' @export

kernel.ornstein.uhlenbeck <- function(l)
{
    storage.mode(l)   <- "double"
    f <- function(x, y=NULL, gradient=FALSE, i = 0) {
        if (gradient) {
            storage.mode(i) <- "integer"
            call.kernel("ornstein_uhlenbeck_gradient", x, y, l,  i)
        }
        else {
            call.kernel("ornstein_uhlenbeck_kernel", x, y, l)
        }
    }
    return (f)
}

#' Create a Matern kernel
#' 
#' @param l length scale
#' @param var noise variance
#' @export

kernel.matern <- function(l, nu)
{
    storage.mode(l)  <- "double"
    storage.mode(nu) <- "double"
    f <- function(x, y=NULL, gradient=FALSE, i = 0) {
        if (gradient) {
            storage.mode(i) <- "integer"
            call.kernel("matern_gradient", x, y, l, nu, i)
        }
        else {
            call.kernel("matern_kernel", x, y, l, nu)
        }
    }
    return (f)
}
