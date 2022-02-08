library(zCompositions)

if(!exists("GEN")) GEN = "count_uniform-size_00100-data_mvtnorm-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

P.rpl = cmultRepl(X)

save(P.rpl, file = sprintf("sim-01/data/replacement_gbm-%s.RData", GEN))