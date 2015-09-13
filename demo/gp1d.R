gp <- new.gp(0.5, kernel.exponential(1, 1))

xp <- c(1, 2, 3)
yp <- c(0.7, 0.7, 0.7)
# measurement noise
ep <- c(0.01, 0.01, 0.01)

marginal.likelihood(gp, xp, yp, ep)

gp <- posterior(gp, xp, yp, ep)

plot(gp, 1:100/20)
