source('~/Documents/voleon_problem1/script2.R')
gp <- new.gp(0, kernel.squared.exponential(1, 1),
#check how this hyperparameter is supposed to come into play.
likelihood=new.likelihood("t", 3,3))
index = 1
result = matrix(NaN, 2^6,1)
K_result = list()
for (alpha in c(1,10)){#check how the loop is organised.
    for (dalpha in c(1,10)){
        for (x in c(1,10)){
            for (y in c(1,10)){#something is wrong here
                for (s in c(1,10)){
                    for (m in c(1,10)){
                        gp$xp = x
                        gp$yp = y
                        K = gp$kernelf(x)#the resulting covariance matrix is different
                        K_result[[index]] = K
                        result[[index]] = approximate.posterior.irls.psi_line(alpha,dalpha,s,m,K,gp)
                        #getting similar results but in a mixedup sequence
                        index = index + 1
                    }
                }
            }
        }
    }
}
print(result)
plot(result)
