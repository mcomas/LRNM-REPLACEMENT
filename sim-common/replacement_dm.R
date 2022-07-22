library(coda.count)

if(!exists("SIM")) GEN = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

P.rpl = t(t(X) + as.vector(fit_dm(X)))
P.rpl = P.rpl / rowSums(P.rpl)

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_dm-%s.RData", SIM, GEN))
