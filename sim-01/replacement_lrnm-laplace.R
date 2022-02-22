library(coda.base)
library(coda.count)
library(zCompositions)

if(!exists("GEN")) GEN = "count_uniform-size_00010-data_mvtnorm-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

X0 = cmultRepl(X, method = 'CZM')
H0 = coordinates(X0)

d = ncol(H0)
M = colMeans(H0)
S = cov(H0)

B0 = matrix(0, nrow = 1+length(M), ncol = length(M))
Bd = lapply(1:ncol(X0), ilr_basis)

Xr = X
Xr[] = 0

for(i in 1:nrow(X)){
  x = X[i,]
  
  B = B0
  iZ = (x == 0)
  sZ = sum(iZ)
  if(sZ > 0){
    sNZ = sum(!iZ)
    B[,1] = sbp_basis(matrix(2*iZ - 1, ncol = 1), silent = TRUE)
    if(sZ > 1){
      B[iZ,2:sZ] = Bd[[sZ]]
    }
    B[!iZ,-(1:sZ)] = Bd[[sNZ]]

    
    if(sNZ > 1){
      Bt = t(B) %*% Bd[[1+d]]
      Mt = Bt %*% M
      St = Bt %*% S %*% t(Bt)
      
      I1 = 1:sZ
      I2 = -I1
      
      h2 = t(Bd[[sNZ]]) %*% log(x[!iZ])
      invSt2 = MASS::ginv(St[I2,I2])
      Mc = Mt[I1] + St[I1,I2] %*% invSt2 %*% (h2-Mt[I2])
      Sc = St[I1,I1] - St[I1,I2] %*% invSt2 %*% St[I2,I1]
      
      h1 = c_posterior_approximation_vec(x, Mc, Sc, B[,I1,drop = FALSE])[,-I1]
      x_r = composition(c(h1,h2), B)
    }else{
      x_r = composition(c_posterior_approximation_vec(x, M, S, ilr_basis(d+1))[,d+1])
    }
  }else{
    x_r = x / sum(x)
  }
  Xr[i,] = x_r
}

P.rpl = Xr

save(P.rpl, file = sprintf("sim-01/data/replacement_lrnm-laplace-%s.RData", GEN))
