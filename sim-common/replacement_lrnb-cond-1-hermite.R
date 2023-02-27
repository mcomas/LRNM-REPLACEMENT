library(coda.base)
library(coda.count)
if(!exists("SIM")) GEN = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

fit_ = fit_conditional_lrnm(X, hermite.order = 50)
p_mu = composition(fit_$mu %*% MASS::ginv(fit_$B1sub), ilr_basis(ncol(X)))
P.rpl = fit_$P
for(i in 1:nrow(X)){
  sel = X[i,] == 0
  nz = sum(sel)
  if(nz > 1){
    P.rpl[i,sel] = P.rpl[i,sel] * prop.table(p_mu[sel]) * nz
  }
}

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_lrnb-cond-1-hermite-%s.RData", SIM, GEN))
