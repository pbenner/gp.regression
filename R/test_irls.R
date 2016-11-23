source('~/Documents/voleon_problem1/script2.R')
gp <- new.gp(0, kernel.squared.exponential(1, 1),
             likelihood=new.likelihood("t", 3, 3))
gp$xp = c(1,1)
gp$yp = c(1,1)
print(approximate.posterior.irls(gp = gp, mean = 0,n = 2))