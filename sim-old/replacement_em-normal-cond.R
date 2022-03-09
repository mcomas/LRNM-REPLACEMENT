library(coda.base)
library(coda.count)
library(zCompositions)
library(mvtnorm)
library(tmvtnorm)
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_parliament-seed_00002"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

d = ncol(X) - 1
M = coordinates(colSums(X))
S = diag(d)

B0 = matrix(0, nrow = 1+length(M), ncol = length(M))
Bd = lapply(1:ncol(X), ilr_basis)


iZ = X == 0
sZ = rowSums(iZ)
sNZ = ncol(X) - sZ

wZ = which(sZ > 0 & sNZ > 1)
wNZ1 = which(sNZ == 1)

lB = lapply(wZ, function(i){
  B = B0
  bs = sapply(which(iZ[i,]), function(iz){
    b = rep(0, d)
    b[iz] = 1
    b[!iZ[i,]] = -1
    b
  })
  B[,1:sZ[i]] = sbp_basis(bs, silent = TRUE)
  B[!iZ[i,],-(1:sZ[i])] = ilr_basis(sNZ[i])
  return(B)
})

lBt = lapply(lB, function(B) t(B) %*% Bd[[1+d]])
lBt.inv = lapply(lBt, function(Bt) MASS::ginv(Bt))

########## Start iteration here
IT = 0
CONT = TRUE
loglikN_prev = NA
BL = apply(X/rowSums(X), 2, function(x) min(x[x > 0]))
while(CONT){
  IT = IT + 1
  
  lMt = lapply(lBt, function(Bt) Bt %*% M)
  lSt = lapply(lBt, function(Bt) Bt %*% S %*% t(Bt))
  
  lh2 = lapply(wZ, function(i) as.vector(t(Bd[[sNZ[i]]]) %*% log(X[i,][!iZ[i,]])))
  
  I1 = lapply(wZ, function(i) 1:sZ[i])
  I2 = lapply(I1, function(i1) -i1)
  
  lInvSt2 = mapply(function(st, i2){
    MASS::ginv(st[i2,i2])
  }, lSt, I2, SIMPLIFY = FALSE)
  lMc = mapply(function(mt, st, invst2, h2, i1, i2){
    as.vector(mt[i1] + st[i1,i2] %*% invst2 %*% (h2-mt[i2]))
  }, lMt, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
  lSc = mapply(function(st, invst2, h2, i1, i2){
    st[i1,i1] - st[i1,i2] %*% invst2 %*% st[i2,i1]
  }, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
  
  lMoments.wZ = mapply(function(i, mc, sc, B){
    print(i)
    xl = X[i,]
    xl[iZ[i,]] = BL[iZ[i,]]
    bl = log(xl) %*% B
    M = mtmvnorm(mc, sc, lower = rep(-Inf, length(mc)), upper = bl[1:sZ[i]])
    M$tvar = M$tvar + M$tmean %*% t(M$tmean)
    M
  }, wZ, lMc, lSc, lB, SIMPLIFY = FALSE)
  
  M1.wZ = mapply(function(moments, h2, Bt.inv){
    Bt.inv %*% c(moments$tmean, h2)
  }, lMoments.wZ, lh2, lBt.inv) |> t()
  
  M2.wZ = mapply(function(moments, h2, Bt.inv){
    S = t(Bt.inv) %*% rbind(
      cbind(moments$tvar, moments$tmean %*% t(h2)),
      cbind(h2 %*% t(moments$tmean), h2 %*% t(h2))) %*% Bt.inv
    0.5 * S + 0.5 * t(S)
  }, lMoments.wZ, lh2, lBt.inv)
  
  if(length(wNZ1) > 0){
    lNapprox.wNZ1 = lapply(wNZ1, function(i){
      c_posterior_approximation_vec(X[i,], M, MASS::ginv(S), ilr_basis(d+1))
    })
    
    lMoments.wNZ1 = mapply(function(i, napprox){
      c_moments_lrnm_hermite(X[i,], napprox[,d+1], napprox[1:d,1:d,drop=FALSE], M, S, ilr_basis(d+1), mu_centering = rep(0,d), order = 10)
    }, wNZ1, lNapprox.wNZ1, SIMPLIFY = FALSE)
    
    M1.wNZ1 = mapply(function(moments){
      moments[,d+1]
    }, lMoments.wNZ1) |> t()
    
    M2.wNZ1 = mapply(function(moments, h2, Bt){
      moments[,1:d]
    }, lMoments.wNZ1)
  }
  
  M1 = matrix(0, nrow = nrow(X), ncol = d)
  M1[-wZ,] = coordinates(X[-wZ,])
  M1[+wZ,] = M1.wZ
  
  M2 = array(0, dim = c(d, d, nrow(X)))
  M2[,,-wZ] = apply(coordinates(X[-wZ,]), 1, function(h) h %*% t(h))
  M2[,,+wZ] = M2.wZ
  
  if(length(wNZ1) > 0){
    M1[wNZ1,] = M1.wNZ1
    M2[,,+wNZ1] = M2.wNZ1
  }
  
  M_new = colMeans(M1)
  S_new = apply(M2, 1:2, mean) - M_new %*% t(M_new)
  tol = max((S - S_new)^2)
  CONT = tol > 1e-04
  
  M = M_new
  S = S_new
  
  loglikN_new = sum(dmvnorm(M1, M, S, log=TRUE))
  print(sprintf("%0.3f. Increment: %0.3f", loglikN_new,  loglikN_new-loglikN_prev))
  loglikN_prev = loglikN_new
}

P.rpl = composition(M1)

save(P.rpl, file = sprintf("sim-01/data/replacement_em-normal-cond-%s.RData", GEN))
