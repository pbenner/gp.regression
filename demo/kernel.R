
gp1 <- new.gp(0.0, kernel.exponential(5, 1))
gp2 <- new.gp(0.0, kernel.ornstein.uhlenbeck(1))
gp3 <- new.gp(0.0, kernel.matern(1/2, 1))

x <- seq(0, 100, length.out=1000)

f1 <- draw.sample(gp1, x, ep=0.000001)
f2 <- draw.sample(gp2, x, ep=0.000001)
f3 <- draw.sample(gp3, x, ep=0.000001)

plot (x, f1, type="l", lty=1, ylim=c(-4,4), ylab="f")
lines(x, f2, type="l", lty=2)
