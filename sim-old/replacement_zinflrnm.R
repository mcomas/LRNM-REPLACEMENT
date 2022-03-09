library(ZIPPCAlnm)

if(!exists("GEN")) GEN = "count_uniform-size_00030-data_mvtnorm-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

P.rpl = ZIPPCAlnm(X)$Q

save(P.rpl, file = sprintf("sim-01/data/replacement_zinflrnm-%s.RData", GEN))
