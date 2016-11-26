source('~/Documents/voleon_problem1/script2.R')
result = matrix(NaN, 2^6,1)
index = 1



for (kernel_hyp2 in c(10^2,1)) {#this guy here needs to be squared because, unlike gpml, gp.regression requires "sigma^2".
for (kernel_hyp1 in c(1,10)) {
      for (lik_hip1 in c(3,13)) {
        for (lik_hip2 in c(3,13)) {
          for (x in c(1,10)) {#at teh 17th iteration, I get different x
            for (y in c(1,10)) {#y = 10 doesn't agree with gpml.
              for (m in c(1,10)) {#when m = 1 and y = 10, it doesn't agree with gpml
                gp <- new.gp(mean, kernel.squared.exponential(kernel_hyp1, kernel_hyp2),
                             likelihood=new.likelihood("t", lik_hip1, lik_hip2))
                gp$xp = c(x)
                gp$yp = c(y)
                result[[index]] = approximate.posterior.irls(gp = gp, mean = m,n = 1)
                index = index + 1
              }
              
            }
            
          }
          
        }
        
      }
    }
}
print(result[1:2])
