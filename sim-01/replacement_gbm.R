library(zCompositions)

if(!exists("GEN")) GEN = "count_uniform-size_100-data_mvtnorm-seed_0001"

###############
load(sprintf("sim-01/%s.RData", GEN))

P.rpl = cmultRepl(X)

save(P.rpl, file = sprintf("sim-01/replacement_gbm-%s.RData", GEN))