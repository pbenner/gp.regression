source('~/Documents/voleon_problem1/script2.R')

kernel_hyp2 = 10
kernel_hyp1 = 1
lik_hip1 = 3
lik_hip2 = 3
x = 1
y = 10
m = 1
gp <- new.gp(m, kernel.squared.exponential(kernel_hyp1, kernel_hyp2),
                           likelihood=new.likelihood("t", lik_hip1, lik_hip2))
gp$xp = c(x)
gp$yp = c(y)

print(gp$kernelf(gp$xp))