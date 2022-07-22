library(coda.count)

if(!exists("SIM")) GEN = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

library(randtoolbox)
Z = sobol(ncol(X) * 100, dim = ncol(X) - 1, normal = TRUE)
P.rpl = fit_lrnm(X, method = 'montecarlo', probs = TRUE, Z = Z, eps = 0.05, max_iter = 100)$P

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_lrnm-montecarlo-%s.RData", SIM, GEN))
