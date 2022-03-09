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

B0 = matrix(0, nrow = 1+length(M), ncol = length(M))
Bd = lapply(1:ncol(X), ilr_basis)


iZ = X == 0
sZ = rowSums(iZ)
sNZ = ncol(X) - sZ

wZ = which(sZ > 0 & sNZ > 1)
wNZ1 = which(sNZ == 1)

P = X
P[X==0] = 0.5

# lB = lapply(wZ, function(i){
#   B = B0
#   part_nz = cbind(as.integer(-!iZ[i,]))
#   for(j in 1:sZ[i]){
#     B[,j] = part_nz
#     B[which(iZ[i,])[j],j] = 1L
#     B[,j] = sbp_basis(B[,j,drop=F], silent = TRUE)
#   }
#   B[!iZ[i,],-(1:sZ[i])] = Bd[[sNZ[i]]]
#   return(B)
# })
lB = lapply(wZ, function(i){
  B = B0
  B[,1] = sbp_basis(matrix(2*iZ[i,] - 1, ncol = 1), silent = TRUE)
  if(sZ[i] > 1){
    B[iZ[i,],2:sZ[i]] = Bd[[sZ[i]]]
  }
  B[!iZ[i,],-(1:sZ[i])] = Bd[[sNZ[i]]]
  return(B)
})

H = mapply(function(i, B){
  coordinates(P[i,], B)
}, wZ, lB) |> t()

lBt = lapply(lB, function(B) t(B) %*% Bd[[1+d]])
lBt.inv = lapply(lBt, function(B) MASS::ginv(B))

I1 = rep(1, length(wZ))


#### iteration GLOBAL
I2 = -I1

lMt = lapply(lBt, function(Bt) as.vector(Bt %*% M))
lSt = lapply(lBt, function(Bt) Bt %*% S %*% t(Bt))

Ncond_all = lapply(1:length(wZ), function(i){
  i1 = 1:sZ[wZ[i]]
  i2 = -i1
  Mt = lMt[[i]]
  St = lSt[[i]]
  h = H[i,]
  St2.inv = MASS::ginv(St[i2,i2])
  Mc = as.vector(Mt[i1] + St[i1,i2,drop=F] %*% St2.inv %*% (h[i2]-Mt[i2]))
  Sc = St[i1,i1] - St[i1,i2] %*% St2.inv %*% St[i2,i1]
  cbind(Sc, Mc)
})

for(i in 1:length(Ncond_all)){
  H[i, 1:sZ[wZ[i]]] = Ncond_all[[i]][,sZ[wZ[i]]+1]
}

Ncond = sapply(1:length(wZ), function(i){
  i1 = 1
  i2 = -i1
  Mt = lMt[[i]]
  St = lSt[[i]]
  h = H[i,]
  St2.inv = MASS::ginv(St[i2,i2])
  Mc = as.vector(Mt[i1] + St[i1,i2,drop=F] %*% St2.inv %*% (h[i2]-Mt[i2]))
  Sc = St[i1,i1] - St[i1,i2] %*% St2.inv %*% St[i2,i1]
  c(Sc[1,1], Mc[1])
}) |> t()

Napprox = sapply(1:length(wZ), function(i){
  i1 = 1
  i2 = -i1
  x = X[wZ[i],]
  h = H[i,]
  mu = l_lrnm_cond_join_maximum(x, Ncond[i,2], 1/Ncond[i,1,drop=F], h[i2], lB[[i]])
  sigma = -1/l_lrnm_cond_join_d2(mu, x, Ncond[i,2], 1/Ncond[i,1,drop=F], h[i2], lB[[i]])
  c(sigma, mu)
}) |> t()

Moments = sapply(1:length(wZ), function(i){
  i1 = 1
  i2 = -i1
  x = X[wZ[i],]
  h = H[i,]
  Mt = lMt[[i]]
  St = lSt[[i]]
  c_moments_lrnm_cond_hermite(x, Napprox[i,2], Napprox[i,1,drop=F], Mt, St, h[i2], lB[[i]], 10, 0)
}) |> t()

# Ncond = sapply(1:length(wZ), function(i){
#   i1 = I1[i]
#   i2 = I2[i]
#   sel = c(1, (sZ[wZ[i]]+1):ncol(St))
#   Mt = lMt[[i]][sel]
#   St = lSt[[i]][sel,sel]
#   h = H[i,sel]
#   St2.inv = MASS::ginv(St[i2,i2])
#   Mc = as.vector(Mt[i1] + St[i1,i2,drop=F] %*% St2.inv %*% (h[i2]-Mt[i2]))
#   Sc = St[i1,i1] - St[i1,i2] %*% St2.inv %*% St[i2,i1]
#   c(Sc[1,1], Mc[1])
# }) |> t()

# Ncond = sapply(1:length(wZ), function(i){
#   i1 = I1[i]
#   i2 = I2[i]
#   Mt = lMt[[i]]
#   St = lSt[[i]]
#   h = H[i,]
#   St2.inv = MASS::ginv(St[i2,i2])
#   Mc = as.vector(Mt[i1] + St[i1,i2,drop=F] %*% St2.inv %*% (h[i2]-Mt[i2]))
#   Sc = St[i1,i1] - St[i1,i2] %*% St2.inv %*% St[i2,i1]
#   c(Sc[1,1], Mc[1])
# }) |> t()

# Napprox = sapply(1:length(wZ), function(i){
#   x = X[wZ[i],]
#   sel = c(1, (sZ[wZ[i]]+1):ncol(St))
#   mu = l_lrnm_cond_join_maximum(x, Ncond[i,2], Ncond[i,1,drop=F], lB[[i]][,1,drop=FALSE])
#   sigma = -1/l_lrnm_join_d2(mu, x, Ncond[i,2], Ncond[i,1,drop=F], lB[[i]][,1,drop=FALSE])
#   c(sigma, mu)
# }) |> t()
# 
# 
# Moments = sapply(1:length(wZ), function(i){
#   x = X[wZ[i],]
#   h = H[i,]
#   Mt = lMt[[i]]
#   St = lSt[[i]]
#   c_moments_lrnm_hermite(h, I1[i], x, Napprox[i,2], Napprox[i,1,drop=F], Mt, St, lB[[i]], 50, 0)
# }) |> t()

# Napprox = sapply(1:length(wZ), function(i){
#   x = X[wZ[i],]
#   h = H[i,]
#   # h[I1[i]] = Ncond[i,2]
#   h[I1[i]] = l_lrnm_cond1_join_maximum(h, I1[i]-1, x, Ncond[i,2], 1/Ncond[i,1,drop=F], lB[[i]])
#   mu = h[I1[i]]
#   sigma = -1/l_lrnm_cond1_join_d2(h, I1[i]-1, x, Ncond[i,2], 1/Ncond[i,1,drop=F], lB[[i]])
#   c(sigma, mu)
# }) |> t()
# 
# 
# Moments = sapply(1:length(wZ), function(i){
#   x = X[wZ[i],]
#   h = H[i,]
#   Mt = lMt[[i]]
#   St = lSt[[i]]
#   c_moments_lrnm_cond1_hermite(h, I1[i], x, Napprox[i,2], Napprox[i,1,drop=F], Mt, St, lB[[i]], 50, 0)
# }) |> t()

M1.wZ = sapply(1:length(wZ), function(i){
  h = H[i,]
  h[1] = Moments[i,2]
  as.vector(lBt.inv[[i]] %*% h)
}) |> t()

M2.wZ = sapply(1:length(wZ), function(i){
  h = H[i,]
  h[1] = Moments[i,2]
  m2 = h %*% t(h)
  m2[1:sZ[wZ[i]],1:sZ[wZ[i]]] = Ncond_all[[i]][,1:sZ[wZ[i]]]
  m2[1,1] = m2[1,1] + Moments[i,1]
  lBt.inv[[i]] %*% m2 %*% t(lBt.inv[[i]])
})

# M1.wZ = sapply(1:length(wZ), function(i){
#   h = H[i,]
#   h[I1[i]] = Moments[i,2]
#   as.vector(lBt.inv[[i]] %*% h)
# }) |> t()
# 
# M2.wZ = sapply(1:length(wZ), function(i){
#   m2 = M1.wZ[i,] %*% t(M1.wZ[i,])
#   m2[I1[i],I1[i]] = Moments[i,1]
#   lBt.inv[[i]] %*% as.matrix(Matrix::nearPD(m2)$mat) %*% t(lBt.inv[[i]])
# })

M1 = matrix(0, nrow = nrow(X), ncol = d)
M1[-wZ,] = coordinates(X[-wZ,])
M1[+wZ,] = M1.wZ

M2 = array(0, dim = c(d, d, nrow(X)))
M2[,,-wZ] = apply(coordinates(X[-wZ,]), 1, function(h) h %*% t(h))
M2[,,+wZ] = M2.wZ

M_new = colMeans(M1)
S_new = apply(M2, 1:2, mean) - M_new %*% t(M_new)
M2c = lapply(1:length(wZ), function(i){
  M2[,,i] - M1[i,] %*% t(M1[i,])
})
#S_new = Matrix::nearPD(apply(M2, 1:2, mean) - M_new %*% t(M_new))$mat |> as.matrix()

tol = max((S - S_new)^2)
CONT = tol > 1e-04
M = M_new
S = S_new

H = M1.wZ
# I1 = (I1 + 1) %% sZ[wZ]
# I1[I1==0] = sZ[wZ][I1==0]
M

