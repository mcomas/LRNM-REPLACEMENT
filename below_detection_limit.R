x = c(0, 0, 0, 0, 0.2, 0.8)
L = c(0.1, 0.2, 0.1, 0.3, 0, 0)
L
d = length(x) - 1
Dz = sum(x == 0)
Dnz = sum(x != 0)
c1 = sbp_basis(t(t(c(1, 1, 1, 1, -1, -1))))
log(x+L) %*% c1

L1 = mean(log(L[1:4])) - mean(log(x[5:6]))
L1
x
B = ilr_basis(6)
coordinates(x+L, B)
Bnz = ilr_basis(2)

B = matrix(0, ncol = d, nrow = d+1)
B[,1] = sbp_basis(t(t(c(1,1,1,1,-1,-1))), silent = TRUE)
B[1:Dz, 2:Dz] = ilr_basis(Dz)
B[-(1:Dz),-(1:Dz)] = ilr_basis(Dnz)
B

M = 0
S = 1

b = coordinates(x+L, B)
beta = (b[1] - M) / sqrt(S)
b1_m1 = M - dnorm(beta) / pnorm(beta) * sqrt(S)
b1_m2 = S * (1 - beta  * dnorm(beta) / pnorm(beta) - (dnorm(beta)/ pnorm(beta))^2 )
