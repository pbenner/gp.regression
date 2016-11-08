source('~/Documents/voleon_problem1/script2.R')
gp <- new.gp(0, kernel.squared.exponential(1, 1),
             likelihood=new.likelihood("t", 3,3))
gp$xp = c(1,1)
gp$yp = c(1,1)
K = matrix(1,2,2)
K[[1,2]] = 2
K[[2,1]] = 3
K[[2,2]] = 4
print(approximate.posterior.psi(gp, c(1,1), K, c(1,1)))