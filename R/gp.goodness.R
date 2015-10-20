
#' Compute the R^2 of a model
#' 
#' @param model regression model
#' @param ... arguments to be passed to methods
#' @export

r.squared <- function(model, ...)
{
    UseMethod("r.squared")
}

#' Compute the R^2 of a Gaussian process
#'
#' @param gp Gaussian process
#' @param xp optional x-values
#' @param yp optional y-values
#' @method r.squared gp
#' @export

r.squared.gp <- function(model, xp = NULL, yp = NULL)
{
    gp <- model

    if (is.null(xp) && is.null(yp)) {
        xp <- gp$xp
        yp <- gp$yp
    }
    stopifnot(!is.null(xp) && !is.null(yp))

    m <- mean(yp)
    f <- summarize(gp, xp)$mean
    ss.tot <- sum((yp - m)^2)
    ss.res <- sum((yp - f)^2)
    1 - ss.res/ss.tot
}

#' Compute the R^2 of a heteroscedastic Gaussian process
#'
#' @param gp Gaussian process
#' @param xp optional x-values
#' @param yp optional y-values
#' @method r.squared gp.heteroscedastic
#' @export

r.squared.gp.heteroscedastic <- function(model, xp = NULL, yp = NULL)
{
    r.squared(gp$gp, xp=xp, yp=yp)
}
