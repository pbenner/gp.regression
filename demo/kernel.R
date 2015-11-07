
gp1 <- new.gp(0.0, kernel.squared.exponential(5, 1))
gp2 <- new.gp(0.0, kernel.ornstein.uhlenbeck(1, 1))
gp3 <- new.gp(0.0, kernel.matern(10, 1, 3/10))

x <- seq(0, 100, length.out=1000)

f1 <- draw.sample(gp1, x, ep=0.000001)
f2 <- draw.sample(gp2, x, ep=0.000001)
f3 <- draw.sample(gp3, x, ep=0.000001)

plot (x, f1, type="l", lty=1, ylim=c(-4,4), ylab="f", col=1)
lines(x, f2+1.5, type="l", col=2)
lines(x, f4-1.5, type="l", col=3)
legend("bottomright",
       c("Squared exponential", "Ornstein-Uhlenbeck", "Matern"), lty=c(1,1,1), col=c(1,2,3))

# combine kernels
# ------------------------------------------------------------------------------

gp <- new.gp(0.0,
             kernel.combined(
                 kernel.squared.exponential(2, 50),
                 kernel.linear(1, 5, 1)))

x <- seq(0, 20, length.out=1000)

f1 <- draw.sample(gp, x, ep=0.000001)
f2 <- draw.sample(gp, x, ep=0.000001)
f3 <- draw.sample(gp, x, ep=0.000001)

plot (x, f2, type="l", lty=1, ylab="f", col=1)
plot (x, f1, type="l", lty=1, ylab="f", col=1, ylim=c(min(f1,f2,f3),max(f1,f2,f3)))
lines(x, f2, type="l", col=2)
lines(x, f3, type="l", col=3)
