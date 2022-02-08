library(coda.base)
library(zCompositions)

if(!exists("GEN")) GEN = "replacement_gbm-count_uniform-size_00100-data_mvtnorm-seed_00001"
FDATA = sub("replacement_(.+)-count_(.+)-size_([0-9]+)-(.+)", "\\4", GEN)


###############
load(sprintf('sim-01/data/%s.RData', FDATA))
load(sprintf("sim-01/data/%s.RData", GEN))

stress = function(H_orig, H_p){
  sqrt(sum((as.matrix(dist(H_p)) -  as.matrix(dist(H_orig)))^2) / sum(as.matrix(dist(H_orig))^2))
}

results = data.frame(
  gen = GEN,
  metric = 'STRESS',
  value = stress(coordinates(P), coordinates(P.rpl)))

saveRDS(results, file = sprintf("sim-01/data/evaluate_stress-%s.RData", GEN))

# evaluate = function(H_gs, H_p){
#   c('paired.dist' = mean(apply(H_gs-H_p, 1, function(x) sqrt(sum(x^2)))),
#     'cov.frobenius' = norm(cov(H_gs) - cov(H_p), type = 'F'),
#     'stress' = sqrt(sum((as.matrix(dist(H_gs)) -  as.matrix(dist(H_p)))^2) / sum(as.matrix(dist(H_gs))^2)))
# }