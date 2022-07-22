library(coda.base)
library(coda.count)
library(randtoolbox)
if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00050-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

fit_ = fit_mv_conditional_lrnm(X, eps = 0.05, max_iter = 100)
P.rpl = fit_$P

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_lrnm-cond-montecarlo-%s.RData", SIM, GEN))
