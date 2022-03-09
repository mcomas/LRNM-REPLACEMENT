library(coda.base)
library(zCompositions)

if(!exists("GEN")) GEN = "replacement_gbm-count_uniform-size_00010-data_mvtnorm-seed_00001"
FDATA = sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\4", GEN)


###############
load(sprintf('sim-01/data/%s.RData', FDATA))
load(sprintf("sim-01/data/%s.RData", GEN))

frobenius = function(H_orig, H_p){
  norm(cov(H_orig) - cov(H_p), type = 'F')
}

results = data.frame(
  gen = GEN,
  metric = 'Frobenius',
  value = frobenius(coordinates(P), coordinates(P.rpl)))

saveRDS(results, file = sprintf("sim-01/data/evaluate_frobenius-%s.RData", GEN))

# evaluate = function(H_gs, H_p){
#   c('paired.dist' = mean(apply(H_gs-H_p, 1, function(x) sqrt(sum(x^2)))),
#     'cov.frobenius' = norm(cov(H_gs) - cov(H_p), type = 'F'),
#     'stress' = sqrt(sum((as.matrix(dist(H_gs)) -  as.matrix(dist(H_p)))^2) / sum(as.matrix(dist(H_gs))^2)))
# }