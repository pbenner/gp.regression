
library(gp.regression)

data("mcycle", package = "MASS")

x <- seq(from=0, to=max(mcycle$times), length.out=100)

# homoscedastic gaussian process
# ------------------------------------------------------------------------------

gp <- new.gp(0.0, kernel.squared.exponential(5, 100))
gp <- posterior(gp, mcycle$times, mcycle$accel, 20.0)

plot(gp, x)

# heteroscedastic gaussian process
# ------------------------------------------------------------------------------

x11(type="cairo")

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

plot(gp, x)

plot(gp$gp.h, x)
