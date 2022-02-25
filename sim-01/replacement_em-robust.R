library(zCompositions)

if(!exists("GEN")) GEN = "count_uniform-size_00030-data_mvtnorm-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

X = X / rowSums(X)
dl = apply(X, 2, function(x) min(x[x!=0]))
P.rpl = lrEM(X, label = 0, dl = dl, rob = TRUE, ini.cov = 'multRepl')

save(P.rpl, file = sprintf("sim-01/data/replacement_em-robust-%s.RData", GEN))