library(MASS)

xp <- beav1$time
yp <- beav1$temp
ep <- 0.1

gp <- new.gp(36.8, kernel.exponential(200, 0.1))
gp <- posterior(gp, xp, yp, ep)

plot(gp, 1:300*10)

marginal.likelihood(gp, xp, yp, ep)
