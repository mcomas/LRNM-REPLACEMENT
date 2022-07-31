library(zCompositions)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "count_uniform-size_00050-data_lrnormal-prop80-dim_3-seed_00001"


###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

if(sum(X==0) > 0){
  P.rpl = cmultRepl(X, method = 'GBM')
}else{
  P.rpl = X / rowSums(X)
}

TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_gbm-%s.RData", SIM, GEN))
