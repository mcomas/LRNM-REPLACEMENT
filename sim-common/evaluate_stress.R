library(coda.base)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "replacement_lrnm-montecarlo-count_uniform-size_00050-data_lrnormal-prop80-dim_5-seed_00001"

FDATA = sub("replacement_.+-data_(.+)", "\\1", GEN)
FCOUNT = sub("replacement_.+-count_(.+)", "\\1", GEN)

###############
load(sprintf('%s/data/count_%s.RData', SIM, FCOUNT))
load(sprintf('%s/data/data_%s.RData', SIM, FDATA))
load(sprintf("%s/data/%s.RData", SIM, GEN))

stress = function(H_orig, H_p){
  sqrt(sum((as.matrix(dist(H_p)) -  as.matrix(dist(H_orig)))^2) / sum(as.matrix(dist(H_orig))^2))
}

iZ = rowSums(X == 0) != 0
results = data.frame(
  gen = GEN,
  metric = 'STRESS',
  value = NA_real_)

if(sum(iZ) >= 5){
  results$value = stress(coordinates(P[iZ,]), coordinates(P.rpl[iZ,]))
}

saveRDS(results, file = sprintf("%s/data/evaluate_stress-%s.rds", SIM, GEN))

# evaluate = function(H_gs, H_p){
#   c('paired.dist' = mean(apply(H_gs-H_p, 1, function(x) sqrt(sum(x^2)))),
#     'cov.frobenius' = norm(cov(H_gs) - cov(H_p), type = 'F'),
#     'stress' = sqrt(sum((as.matrix(dist(H_gs)) -  as.matrix(dist(H_p)))^2) / sum(as.matrix(dist(H_gs))^2)))
# }