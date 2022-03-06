library(coda.base)
D = 5
d = D-1
M0 = rep(1, d)
S0 = diag(d)
B0 = ilr_basis(D)

x = c(1, 0, 2, 1, 1)
L = c(0, 0.04, 0, 0, 0)

iZ = x == 0

Dz = sum(iZ)
Dnz = sum(!iZ)

B = matrix(0, ncol = d, nrow = D)
part = cbind(2 * iZ - 1)
B[,1] = sbp_basis(part, silent = TRUE)
if(Dz > 1){
  B[iZ, 2:Dz] = ilr_basis(Dz)
}
B[!iZ,-(1:Dz)] = ilr_basis(Dnz)

bl = log(x+L) %*% B

Bt = t(B) %*% B0

M = Bt %*% M0
S = Bt %*% S0 %*% t(Bt)

i1 = 1
i2 = 2:d
invS2 = MASS::ginv(S[i2,i2])
M1 = as.vector(M[i1] + S[i1,i2] %*% invS2 %*% (bl[i2]-M[i2]))
S1 = S[i1,i1] - S[i1,i2] %*% invS2 %*% S[i2,i1]


beta = (bl[i1] - M1) / sqrt(S1)
b1_m1 = M1 - dnorm(beta) / pnorm(beta) * sqrt(S1)
b1_m2 = S1 * (1 - beta  * dnorm(beta) / pnorm(beta) - (dnorm(beta)/ pnorm(beta))^2 )

hr = c(b1_m1, bl[i2])
xr = composition(hr, B)
xr * x[1] / xr[1]



### one-dim
part = c(-1,1)
B = sbp_basis(cbind(part), silent = TRUE)
bl = log(x+L) %*% B
bl

Bt = t(B) %*% B0

M = Bt %*% M0
S = Bt %*% S0 %*% t(Bt)


beta = (bl[1] - M) / sqrt(S)
b1_m1 = M - dnorm(beta) / pnorm(beta) * sqrt(S)
b1_m2 = S * (1 - beta  * dnorm(beta) / pnorm(beta) - (dnorm(beta)/ pnorm(beta))^2 )

xr = composition(b1_m1, B)
xr * x[1] / xr[1]

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

##



