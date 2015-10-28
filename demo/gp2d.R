
np <- 400
xp <- cbind(x1 = runif(np), x2 = runif(np))
yp <- sin(pi*xp[,1]) + cos(2*pi*xp[,2]) + rnorm(np, 0, 1)

gp <- new.gp(0.5, kernel.squared.exponential(0.5, 1), dim=2)
gp <- posterior(gp, xp, yp, 1)

# plot both dimensions
x  <- as.matrix(expand.grid(x = 1:100/100, y = 1:100/100))
plot(gp, x, plot.scatter=TRUE, plot.variance=FALSE)

# reduce plot to one dimension
x  <- as.matrix(expand.grid(x = 1:100/100, y = 0.0))
plot(gp, x, plot.scatter=TRUE, plot.variance=FALSE, covariates=1)
