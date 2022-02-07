library(zCompositions)

if(!exists("GEN")) GEN = "count_uniform-size_00100-data_mvtnorm-seed_00001"

###############
load(sprintf("sim-01/%s.RData", GEN))

P.rpl = cmultRepl(X, method = 'CZM')

save(P.rpl, file = sprintf("sim-01/replacement_czm-%s.RData", GEN))