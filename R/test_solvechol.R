test_solvechol <- function(W_vector, K, r, n){
rootW <- sqrt(W_vector)
B <- diag(n) + rootW %*% t(rootW) * K # (c.f. Rasmussen 2006, Eq. 3.26)
L <- chol(B)
#temp <- solve(L * t(L), sweep(r, 2, rootW, "*"))#this sweep is strange, r and rootW should be both vectors.
temp <- solve(L, solve(t(L), sweep(r, 1, rootW, "*"))
}