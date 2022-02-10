library(coda.base)
library(coda.count)
library(zCompositions)
library(mclust)

if(!exists("GEN")) GEN = "count_uniform-size_00010-data_parliament-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

X0 = cmultRepl(X, method = 'CZM')
B0 = pc_basis(X0)
H0 = coordinates(X0, B0)

d = ncol(H0)
mixt = Mclust(H0, modelNames = "EII")

B0 = matrix(0, nrow = 1+ncol(H0), ncol = ncol(H0))
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
      I1 = 1:sZ
      I2 = -I1
      
      h2 = t(Bd[[sNZ]]) %*% log(x[!iZ])
      
      # post = prop.table(round(mixt$z[i,], 4))
      # ipost = post > 0
      # l_h1 = lapply((1:mixt$G)[ipost], function(g){
      #   Mt = Bt %*% mixt$parameters$mean[,g]
      #   St = t(Bt) %*% mixt$parameters$variance$sigma[,,g] %*% Bt
      #   
      #   invSt2 = MASS::ginv(St[I2,I2])
      #   Mc = Mt[I1] + St[I1,I2] %*% invSt2 %*% (h2-Mt[I2])
      #   Sc = St[I1,I1] - St[I1,I2] %*% invSt2 %*% St[I2,I1]
      #   
      #   c_posterior_approximation_vec(x, Mc, Sc, B[,I1,drop = FALSE])[,-I1]
      # })
      # h1 = rowSums(matrix(mapply(`*`, post[ipost], l_h1), ncol = sum(ipost)))
      l_post = lapply(1:mixt$G, function(g){
        Mt = Bt %*% mixt$parameters$mean[,g]
        St = t(Bt) %*% mixt$parameters$variance$sigma[,,g] %*% Bt
        
        invSt2 = MASS::ginv(St[I2,I2])
        Mc = Mt[I1] + St[I1,I2] %*% invSt2 %*% (h2-Mt[I2])
        Sc = St[I1,I1] - St[I1,I2] %*% invSt2 %*% St[I2,I1]
        
        N_pars = c_posterior_approximation_vec(x, Mc, Sc, B[,I1,drop = FALSE])
        
        h = N_pars[,-I1]
        log_p_h = dmvnorm(rbind(h), Mc, Sc, log = TRUE)
        log_px_h = dmultinom(x, prob = composition(c(h, h2), B), log = TRUE)
        log_ph_x = 0.5 * log(1/det(N_pars[,I1,drop=FALSE]))
        list(p = mixt$parameters$pro[g] * exp(log_px_h + log_p_h  - log_ph_x),
             m = h)
      })
      h1 = mapply(`*`, 
                  prop.table(sapply(l_post, `[[`, 1)),
                  lapply(l_post, `[[`, 2)) |>
        matrix(nrow = sZ) |>
        rowSums()
      
      x_r = composition(c(h1,h2), B)
    }else{
      warning("Posterior calculation needs to be modified.")
      post = prop.table(round(mixt$z[i,], 4))
      ipost = post > 0
      l_h1 = lapply((1:mixt$G)[ipost], function(g){
        c_posterior_approximation_vec(x, mixt$parameters$mean[,g], mixt$parameters$variance$sigma[,,g], B0)[,d+1]
      })
      h1 = rowSums(mapply(`*`, post[ipost], l_h1))
      x_r = composition(h1)
    }
  }else{
    x_r = x / sum(x)
  }
  Xr[i,] = x_r
}

P.rpl = Xr

save(P.rpl, file = sprintf("sim-01/data/replacement_mixture-lrn-laplace-%s.RData", GEN))
