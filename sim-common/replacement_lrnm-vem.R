library(coda.count)
library(coda.base)
if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00050-data_lrnormal-prop80-dim_5-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

P.rpl = fit_vem_lrnm(X)

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_lrnm-vem-%s.RData", SIM, GEN))
