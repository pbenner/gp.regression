source('~/Documents/Scripts/r/gp.regression/R/gp.approximation.R')
library(Matrix)
W = c(-1,2)
K = matrix(c(1,3,2,4),2,2)
print(approximate.posterior.irls.ldB2_exact(W,K,2) )#matches with matlab
W = c(1,2)
print(approximate.posterior.irls.ldB2_exact(W,K,2) )#matches with matlab
