x  <- as.matrix(expand.grid(x = 1:20/4, y = 1:20/4))
gp <- new.gp(x, 0.5, kernel.exponential(2, 1))

xp <- matrix(0,2, nrow=2)
xp[1,] <- c(2,2)
xp[1,] <- c(4,4)

yp <- c(0.2, 0.8)
# measurement noise
ep <- c(0.001, 0.001)

gp <- posterior(gp, xp, yp, ep)

plot(gp)
