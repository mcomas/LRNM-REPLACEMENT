library(zCompositions)


if(!exists("GEN")) GEN = "count_uniform-size_00030-data_lrnormal-dim_3-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

P.rpl = cmultRepl(X, method = 'GBM')

save(P.rpl, file = sprintf("sim-01/data/replacement_gbm-%s.RData", GEN))
