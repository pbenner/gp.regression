
library(gp.regression)

data("mcycle", package = "MASS")

# homoscedastic gaussian process
# ------------------------------------------------------------------------------

gp <- new.gp(0.0, kernel.exponential(5, 100))
gp <- posterior(gp, mcycle$times, mcycle$accel, 20.0)

p <- plot(gp, seq(from=0, to=max(mcycle$times), length.out=100))
p

# heteroscedastic gaussian process
# ------------------------------------------------------------------------------

x11(type="cairo")

gp <- new.gp.heteroscedastic(
    new.gp(1.0, kernel.exponential(5, 100)),
    new.gp(0.0, kernel.exponential(5, 0.1)))

gp <- posterior(gp, mcycle$times, mcycle$accel, 0.1,
                step = 0.1,
                epsilon = 0.0001,
                verbose=T)

p <- plot(gp, seq(from=0, to=max(mcycle$times), length.out=100))
p

# maximize marginal likelihood
# ------------------------------------------------------------------------------

gp <- new.gp(0.0, kernel.exponential(5, 100))

gp <- maximize.marginal.likelihood(gp, kernel.exponential, c(10,1),
                                   mcycle$times, mcycle$accel, 20.0,
                                   step.init = 0.0001, eta = 0.01,
                                   verbose=TRUE)

gp <- posterior(gp, mcycle$times, mcycle$accel, 20.0)

p <- plot(gp, seq(from=0, to=max(mcycle$times)))
p
