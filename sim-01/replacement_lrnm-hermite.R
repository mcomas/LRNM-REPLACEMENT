library(coda.count)


if(!exists("GEN")) GEN = "count_uniform-size_00030-data_parliament-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))

# if(ncol(X) < 6){
  P.rpl = fit_lrnm(X, method = 'hermite', probs = TRUE)$P
# }else{
#   library(randtoolbox)
#   Z = sobol(500, dim = ncol(X) - 1, normal = TRUE)
#   P.rpl = fit_lrnm(X, method = 'montecarlo', probs = TRUE, Z = Z)$P
# }


save(P.rpl, file = sprintf("sim-01/data/replacement_lrnm-hermite-%s.RData", GEN))
