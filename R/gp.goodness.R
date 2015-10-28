
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
#' @param adjusted compute adjusted R^2
#' @method r.squared gp
#' @export

r.squared.gp <- function(model, xp = NULL, yp = NULL, adjusted = FALSE)
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
    r2 <- 1 - ss.res/ss.tot

    if (adjusted) {
        n  <- nrow(gp$xp)
        k  <- dim (gp)
        stopifnot(n-k-1 > 0)
        r2 <- 1.0 - (1.0-r2)*(n-1.0)/(n-k-1)
    }
    r2
}

#' Compute the R^2 of a heteroscedastic Gaussian process
#'
#' @param gp Gaussian process
#' @param xp optional x-values
#' @param yp optional y-values
#' @param ... options passed to the r.squared method
#' @method r.squared gp.heteroscedastic
#' @export

r.squared.gp.heteroscedastic <- function(model, xp = NULL, yp = NULL, ...)
{
    r.squared(gp$gp, xp=xp, yp=yp, ...)
}
