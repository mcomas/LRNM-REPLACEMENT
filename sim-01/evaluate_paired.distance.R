library(coda.base)
library(zCompositions)

if(!exists("GEN")) GEN = "replacement_mixture-lrn-laplace-count_uniform-size_00010-data_parliament-seed_00001"
FDATA = sub("replacement_.+-data_(.+)", "\\1", GEN)
FCOUNT = sub("replacement_.+-count_(.+)", "\\1", GEN)

###############
load(sprintf('sim-01/data/count_%s.RData', FCOUNT))
load(sprintf('sim-01/data/data_%s.RData', FDATA))
load(sprintf("sim-01/data/%s.RData", GEN))

paired.distance = function(H_orig, H_p){
  mean(apply(H_orig-H_p, 1, function(x) sqrt(sum(x^2))))
}

iNZ = rowSums(X == 0) != 0

results = data.frame(
  gen = GEN,
  metric = 'Distance',
  value = paired.distance(coordinates(P[iNZ,]), coordinates(P.rpl[iNZ,])))

saveRDS(results, file = sprintf("sim-01/data/evaluate_paired.distance-%s.RData", GEN))

# evaluate = function(H_gs, H_p){
#   c('paired.dist' = mean(apply(H_gs-H_p, 1, function(x) sqrt(sum(x^2)))),
#     'cov.frobenius' = norm(cov(H_gs) - cov(H_p), type = 'F'),
#     'stress' = sqrt(sum((as.matrix(dist(H_gs)) -  as.matrix(dist(H_p)))^2) / sum(as.matrix(dist(H_gs))^2)))
# }