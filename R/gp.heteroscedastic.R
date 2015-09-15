
#' Generate a new heteroscedastic Gaussian process
#' 
#' @param gp Gaussian process on the observations
#' @param gp.h Gaussian process on empirical variances
#' @param transform transformation of empirical variances
#' @param transform.inv inverse transformation of empirical variances
#' @export

new.gp.heteroscedastic <- function(gp, gp.h, transform = identity, transform.inv = identity)
{
    stopifnot(length(mean) == 1)

    gp <- list(gp   = gp,
               gp.h = gp.h,
               transform     = transform,
               transform.inv = transform.inv)

    class(gp) <- "gp.heteroscedastic"

    return (gp)
}

#' Compute posterior of a Gaussian process
#' 
#' @param model Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param ep uncertainty of measurements (optional)
#' @param ep.init initial noise estimate
#' @param epsilon terminate iteration if error is below epsilon
#' @param step step size for updating the mean
#' @param verbose print mean squared error in each iteration
#' @param ... unused
#' @method posterior gp.heteroscedastic
#' @export

posterior.gp.heteroscedastic <- function(model, xp, yp, ep = NULL,
                                         ep.init = 1.0,
                                         epsilon = 0.01,
                                         step = 0.2,
                                         verbose = FALSE, ...)
{
    # check arguments
    if (!is.null(xp) && !is.matrix(xp)) {
        xp <- as.matrix(xp, 1)
    }
    if (!is.null(yp) && !is.vector(yp)) {
        yp <- as.vector(yp)
    }
    if (!is.null(ep) && !is.vector(ep)) {
        ep <- as.vector(ep)
    }
    if (!is.null(ep) && !is.vector(ep)) {
        ep.init <- as.vector(ep.init)
    }
    # empirical variance for each data point
    empirical.variance <- function(mean, data) {
        sapply(mean-data, function(x) sum(x^2))
    }
    # new gp.1 evaluated at measurement locations xp
    gp.1 <- model$gp
    gp.1 <- posterior(gp.1, xp, yp, ep=ep.init, verbose=verbose...)
    # keep a copy of the mean
    mean <- summarize(gp.1, xp, ...)$mean
    repeat {
        # recompute empirical variance
        ep.empirical <- model$transform(empirical.variance(mean, yp))
        # reset gp.h
        gp.h <- model$gp.h
        gp.h <- posterior(gp.h, xp, ep.empirical, ep=ep, verbose=verbose, ...)
        # recompute gp.2
        gp.2 <- model$gp
        gp.2 <- posterior(gp.2, xp, yp, ep=model$transform.inv(summarize(gp.h, xp, ...)$mean), verbose=verbose, ...)
        # compare gp.2 with gp.1
        gp.1.mean <- summarize(gp.1, xp, ...)$mean
        gp.2.mean <- summarize(gp.2, xp, ...)$mean
        error <- max(abs(gp.1.mean-gp.2.mean))
        if (verbose) {
            print(sprintf("Training heteroscedastic model... (error: %f)", error))
        }
        if (error < epsilon) {
            # return new model
            model$gp   <- gp.1
            model$gp.h <- gp.h
            return (model)
        }
        mean <- (1.0-step)*mean + step*gp.2.mean
        gp.1 <- gp.2
    }
}

#' Compute posterior expectation and variance of a heteroscedastic GP
#' 
#' @param model GP object
#' @param x where to evaluate the posterior
#' @param ... unused
#' @method summarize gp.heteroscedastic
#' @export

summarize.gp.heteroscedastic <- function(model, x, ...)
{
    s <- summarize(model$gp, x)
    # add noise estimate to gp
    list(mean     = s$mean,
         variance = s$variance + model$transform.inv(summarize(model$gp.h, x, ...)$mean)
         )
}

#' Compute marginal likelihood of a heteroscedastic Gaussian process
#' 
#' @param model heteroscedastic Gaussian process
#' @param xp positions of measurements
#' @param yp measured values
#' @param ... unused
#' @method marginal.likelihood gp.heteroscedastic
#' @export

marginal.likelihood.gp.heteroscedastic <- function(model, xp, yp, ...)
{
    ep <- summarize(model$gp.h, xp, ...)$mean
    marginal.likelihood(model$gp, xp, yp, ep=ep, ...)
}

# ------------------------------------------------------------------------------
if (FALSE) {

    data("mcycle", package = "MASS")

    x <- seq(from=0, to=max(mcycle$times), length.out=100)

    # --------------------------------------------------------------------------
    gp <- new.gp.heteroscedastic(
        new.gp(1.0, kernel.squared.exponential(5, 100)),
        new.gp(0.0, kernel.squared.exponential(5, 0.1)))
    gp <- posterior(gp, mcycle$times, mcycle$accel, 0.1,
                    step = 0.1,
                    epsilon = 0.0001,
                    verbose=T)
    # --------------------------------------------------------------------------
    gp <- new.gp.heteroscedastic(
        new.gp( 0.0, kernel.squared.exponential(4, 100)),
        new.gp(10.0, kernel.squared.exponential(4,  10),
               likelihood=new.likelihood("gamma", 1),
               link=new.link("logistic")),
        transform     = sqrt,
        transform.inv = function(x) x^2)
    gp <- posterior(gp, mcycle$times, mcycle$accel, 0.00001,
                    step = 0.1,
                    epsilon = 0.000001,
                    verbose=T)
    # --------------------------------------------------------------------------

    x11(type="cairo")

    plot(gp, x)

    plot(gp$gp.h, x)

}
