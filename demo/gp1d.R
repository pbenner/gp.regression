library(MASS)

xp <- beav1$time
yp <- beav1$temp
ep <- 0.1

gp <- new.gp(36.8, kernel.exponential(200, 0.1))
gp <- posterior(gp, xp, yp, ep)

plot(gp, 1:300*10)

# compute the marginal likelihood of the data
# ------------------------------------------------------------------------------

gp <- new.gp(36.8, kernel.exponential(200, 0.1))

marginal.likelihood(gp, xp, yp, ep)

# draw samples from the prior distribution
# ------------------------------------------------------------------------------

gp <- new.gp(36.8, kernel.exponential(200, 0.1))

plot (x, draw.sample(gp, x, ep=0.000001), type="l", lty=1, ylim=c(36,38))
lines(x, draw.sample(gp, x, ep=0.000001), type="l", lty=2)
lines(x, draw.sample(gp, x, ep=0.000001), type="l", lty=3)
