library(coda.count)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00050-data_lrnormal-prop80-dim_20-seed_00004"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

fit = fit_lrnm(X, eps = 0.05, mc.nsim = 500, max.iter = 500)
P.rpl = fit$P

TIME = proc.time() - t0
save(P.rpl, TIME, fit, file = sprintf("%s/data/replacement_lrnm-montecarlo-%s.RData", SIM, GEN))
