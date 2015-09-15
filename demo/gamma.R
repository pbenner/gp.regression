
gp <- new.gp(1.0, kernel.squared.exponential(1.0, 5.0),
             likelihood=new.likelihood("gamma", 1.0),
             link=new.link("logistic"))

n  <- 1000
xp <- 10*runif(n)
yp <- rgamma(n, 1, 2)

# add some tiny noise to the diagonal for numerical stability
gp <- posterior(gp, xp, yp, ep=0.01, verbose=TRUE)

plot(gp, 1:100/10)

