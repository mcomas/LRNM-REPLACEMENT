library(coda.count)


if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))


library(randtoolbox)
Z = sobol(ncol(X) * 100, dim = ncol(X) - 1, normal = TRUE)
P.rpl = fit_lrnm(X, method = 'montecarlo', probs = TRUE, Z = Z, eps = 0.05, max_iter = 100)$P



save(P.rpl, file = sprintf("sim-01/data/replacement_lrnm-montecarlo-%s.RData", GEN))
