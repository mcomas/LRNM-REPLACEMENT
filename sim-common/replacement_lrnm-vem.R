library(coda.count)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "count_dim-size_00030-data_lrnormal-prop80-dim_3-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

fit = fit_lrnm(X, method = 'vem', max_iter = 1000, eps = 0.01)
P.rpl = fit$P

TIME = proc.time() - t0
save(P.rpl, TIME, fit, file = sprintf("%s/data/replacement_lrnm-vem-%s.RData", SIM, GEN))
