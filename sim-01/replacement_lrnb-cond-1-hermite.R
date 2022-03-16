library(coda.base)
library(coda.count)
library(zCompositions)
library(mvtnorm)
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-seed_00002"

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
wZany = c(wZ,wNZ1)

lB_all = lapply(1:nrow(X), function(i){
  if(!i %in% wZany){
    return(NULL)
  }
  nz = sZ[i]
  B = matrix(0, nrow = d+1, ncol = d)
  if(nz>1){
    B[iZ[i,],1:(nz-1)] = ilr_basis(nz)
  }
  B[,nz] = sbp_basis(matrix(2*iZ[i,] - 1, ncol = 1), silent = TRUE)
  B[!iZ[i,],-(1:nz)] = Bd[[sNZ[i]]]
  return(B)
})

for(K in 0){
  wZ = wZ[sZ[wZ] > K]
  wNZ1 = wNZ1[sZ[wNZ1] > K]
  
  ### Big iteration
  if(length(wZ) > 0){
    if(K==0){
      lB = lapply(wZ, function(i){
        lB_all[[i]][,(sZ[i]-K):d]
      })
    }else{
      lB = lB_all[wZ]
    }
    lBt = lapply(lB, function(B) t(B) %*% Bd[[1+d]])
    lBt.inv = lapply(lBt, function(B) MASS::ginv(B))
  }
  if(length(wNZ1) > 0){
    if(K==0){
      lBNZ1 = lapply(wNZ1, function(i){
        lB_all[[i]][,(sZ[i]-K):d,drop=FALSE]
      })
    }else{
      lBNZ1 = lB_all[wNZ1]
    }
    lBNZ1t = lapply(lBNZ1, function(B) t(B) %*% Bd[[1+d]])
    lBNZ1t.inv = lapply(lBNZ1t, function(B) MASS::ginv(B))
  }
  ########## Start iteration here
  IT = 0
  CONT = TRUE
  loglikN_prev = NA
  while(CONT){
    IT = IT + 1
    
    #### iteration GLOBAL
    if(length(wZ) > 0){
      lMt = lapply(lBt, function(Bt) as.vector(Bt %*% M))
      lSt = lapply(lBt, function(Bt) Bt %*% S %*% t(Bt))
      
      if(K == 0){
        lh2 = lapply(wZ, function(i) as.vector(t(Bd[[sNZ[i]]]) %*% log(X[i,][!iZ[i,]])))
      }else{
        lh2 = lapply(1:length(wZ), function(i) as.vector(lBt[[i]] %*% M1[wZ[i],])[-1])
      }
      
      if(K == 0){
        I1 = 1
        I2 = -I1
      }else{
        I1 = lapply(lh2, function(h2) 1:(d-length(h2)))
        I2 = lapply(I1, function(i1) -i1)
      }
      
      lInvSt2 = mapply(function(st, i2){
        MASS::ginv(st[i2,i2])
      }, lSt, I2, SIMPLIFY = FALSE)
      lMc = mapply(function(mt, st, invst2, h2, i1, i2){
        as.vector(mt[i1] + st[i1,i2] %*% invst2 %*% (h2-mt[i2]))
      }, lMt, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
      lSc = mapply(function(st, invst2, h2, i1, i2){
        st[i1,i1] - st[i1,i2] %*% invst2 %*% st[i2,i1]
      }, lSt, lInvSt2, lh2, I1, I2, SIMPLIFY = FALSE)
      
      if(K == 0){
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
      }else{
        M1.wZ = mapply(function(mc, h2, Bt){
          t(Bt) %*% c(mc, h2)
        }, lMc, lh2, lBt) |> t()
        
        M2.wZ = mapply(function(mc, sc, h2, Bt){
          t(Bt) %*% rbind(
            cbind(sc + mc %*% t(mc), mc %*% t(h2)),
            cbind(h2 %*% t(mc), h2 %*% t(h2))) %*% Bt
        }, lMc, lSc, lh2, lBt)
      }
    }
    
    if(length(wNZ1) > 0){

      lMt = lapply(lBNZ1t, function(Bt) as.vector(Bt %*% M))
      lSt = lapply(lBNZ1t, function(Bt) Bt %*% S %*% t(Bt))
      
      if(K == 0){
        lNapprox.wNZ1 = mapply(function(i, mt, st, B){
          c_posterior_approximation_vec(X[i,], mt, MASS::ginv(st), B)
        }, wNZ1, lMt, lSt, lBNZ1, SIMPLIFY = FALSE)
        
        lMoments.wNZ1 = mapply(function(i, napprox, mt, st, B){
          c_moments_lrnm_hermite(X[i,], napprox[,2], napprox[1,1,drop=FALSE], mt, st, B, mu_centering = 0, order = 10)
        }, wNZ1, lNapprox.wNZ1, lMt, lSt, lBNZ1, SIMPLIFY = FALSE)
        
        M1.wNZ1 = mapply(function(moments, Bt){
          t(Bt) %*% moments[,2]
        }, lMoments.wNZ1, lBNZ1t) |> t()
        
        M2.wNZ1 = mapply(function(moments, Bt){
          t(Bt) %*% moments[,1] %*% Bt
        }, lMoments.wNZ1, lBNZ1t)
      }else{
        lh2 = lapply(1:length(wNZ1), function(i) as.vector(lBNZ1t[[i]] %*% M1[wNZ1[i],])[-1])
        
        I1 = lapply(lh2, function(h2) 1:(d-length(h2)))
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
        
        M1.wNZ1 = mapply(function(mc, h2, Bt){
          t(Bt) %*% c(mc, h2)
        }, lMc, lh2, lBNZ1t) |> t()
        
        M2.wNZ1 = mapply(function(mc, sc, h2, Bt){
          t(Bt) %*% rbind(
            cbind(sc + mc %*% t(mc), mc %*% t(h2)),
            cbind(h2 %*% t(mc), h2 %*% t(h2))) %*% Bt
        }, lMc, lSc, lh2, lBNZ1t)
      }
    }
    
    
    if(K == 0){
      M1 = matrix(0, nrow = nrow(X), ncol = d)
      M1[-wZ,] = coordinates(X[-wZ,])
    }
    if(K == 0){
      M2 = array(0, dim = c(d, d, nrow(X)))
      M2[,,-wZ] = apply(coordinates(X[-wZ,]), 1, function(h) h %*% t(h))
    }
    if(length(wZ) > 0){
      M1[+wZ,] = M1.wZ
      M2[,,+wZ] = M2.wZ
    }
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
    
  }
}

P.rpl = composition(M1)

save(P.rpl, file = sprintf("sim-01/data/replacement_lrnb-cond-1-hermite-%s.RData", GEN))
