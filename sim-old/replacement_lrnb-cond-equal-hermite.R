library(coda.base)
library(coda.count)
library(zCompositions)
library(mvtnorm)
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_parliament-seed_00004"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

d = ncol(X) - 1
M = coordinates(colSums(X))
S = diag(d)

Bd = lapply(1:ncol(X), ilr_basis)


iZ = X == 0
sZ = rowSums(iZ)
sNZ = ncol(X) - sZ

wZ = which(sZ > 0 & sNZ > 1)
wNZ1 = which(sNZ == 1)

lB = lapply(wZ, function(i){
  B = matrix(0, nrow = 1+length(M), ncol = sNZ[i])
  B[,1] = sbp_basis(matrix(2*iZ[i,] - 1, ncol = 1), silent = TRUE)
  B[!iZ[i,],-1] = Bd[[sNZ[i]]]
  return(B)
})

lBt = lapply(lB, function(B) t(B) %*% Bd[[1+d]])
lBt.inv = lapply(lBt, function(B) MASS::ginv(B))

########## Start iteration here
IT = 0
CONT = TRUE
loglikN_prev = NA
while(CONT){
  IT = IT + 1
  
  #### iteration GLOBAL
  
  lMt = lapply(lBt, function(Bt) as.vector(Bt %*% M))
  lSt = lapply(lBt, function(Bt) Bt %*% S %*% t(Bt))
  
  lh2 = lapply(wZ, function(i) as.vector(t(Bd[[sNZ[i]]]) %*% log(X[i,][!iZ[i,]])))
  
  I1 = 1
  I2 = -I1
  
  lInvSt2 = mapply(function(st, i2){
    MASS::ginv(st[i2,i2])
  }, lSt, I2, SIMPLIFY = FALSE)
  lMc = mapply(function(mt, st, invst2, h2, i1, i2){
    as.vector(mt[i1] + st[i1,i2] %*% invst2 %*% (h2-mt[i2]))
  }, lMt, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
  lSc = mapply(function(st, invst2, h2, i1, i2){
    st[i1,i1] - st[i1,i2] %*% invst2 %*% st[i2,i1]
  }, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
  
  lNapprox.wZ = mapply(function(i, mc, sc, h2, B){
    c_lrnm_cond_posterior_approximation_vec(X[i,], mc, MASS::ginv(sc), h2, B)
  }, wZ, lMc, lSc, lh2, lB, SIMPLIFY = FALSE)
  
  lMoments.wZ = mapply(function(i, napprox, mt, st, h2, B){
    c_moments_lrnm_cond_hermite(X[i,], napprox[,2], napprox[1,1,drop=FALSE], mt, st, h2, B, mu_centering = 0, order = 10)
  }, wZ, lNapprox.wZ, lMt, lSt, lh2, lB, SIMPLIFY = FALSE)
  
  M1.wZ = mapply(function(moments, h2, Bt){
    t(Bt) %*% c(moments[,2], h2)
  }, lMoments.wZ, lh2, lBt) |> t()
  
  M2.wZ = mapply(function(moments, h2, Bt){
    t(Bt) %*% rbind(
      cbind(moments[,1], moments[,2] %*% t(h2)),
      cbind(h2 %*% t(moments[,2]), h2 %*% t(h2))) %*% Bt
  }, lMoments.wZ, lh2, lBt)
  
  if(length(wNZ1) > 0){
    #stop("Observations with all except one zero")
    lB.wNZ1 = lapply(1:length(wNZ1), function(i){
      B = matrix(0, nrow = 1+length(M), ncol = sNZ[wNZ1[i]])
      B[,1] = sbp_basis(matrix(2*iZ[wNZ1[i],] - 1, ncol = 1), silent = TRUE)
      B
    })
    lBt.wNZ1 = lapply(lB.wNZ1, function(B){
      t(B) %*% Bd[[1+d]]
    })
    lMt.wNZ1 = lapply(lBt.wNZ1, function(Bt) as.vector(Bt %*% M))
    lSt.wNZ1 = lapply(lBt.wNZ1, function(Bt) Bt %*% S %*% t(Bt))

    
    
    lNapprox.wNZ1 = lapply(1:length(wNZ1), function(i){
      mt = lMt.wNZ1[[i]]
      st =  lSt.wNZ1[[i]]
      c_posterior_approximation_vec(X[wNZ1[i],], mt, 1/st, lB.wNZ1[[i]])
    })
    lMoments.wNZ1 = lapply(1:length(wNZ1), function(i){
      c_moments_lrnm_hermite(X[wNZ1[i],], lNapprox.wNZ1[[i]][,2], lNapprox.wNZ1[[i]][1,1,drop=FALSE], lMt.wNZ1[[i]], lSt.wNZ1[[i]], lB.wNZ1[[i]], order = 10, mu_centering = 0)
      
    })
    
    
    # 
    # lMoments.wNZ1 = mapply(function(i, napprox){
    #   c_moments_lrnm_hermite(X[i,], napprox[,d+1], napprox[1:d,1:d,drop=FALSE], M, S, ilr_basis(d+1), mu_centering = rep(0,d), order = 10)
    # }, wNZ1, lNapprox.wNZ1, SIMPLIFY = FALSE)
    # 
    M1.wNZ1 = mapply(function(moments){
      moments[,2] %*% lBt.wNZ1[[1]]
    }, lMoments.wNZ1) |> t()
    # 
    M2.wNZ1 = mapply(function(moments, h2, Bt){
      t(lBt.wNZ1[[1]]) %*% moments[,1,drop=F] %*% lBt.wNZ1[[1]]
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
  
  H = M1.wZ
}

P.rpl = composition(M1)

save(P.rpl, file = sprintf("sim-01/data/replacement_lrnb-cond-equal-hermite-%s.RData", GEN))
