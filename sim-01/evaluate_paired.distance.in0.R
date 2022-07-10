library(coda.base)

if(!exists("GEN")) GEN = "replacement_zinflrnm-count_uniform-size_00050-data_parliament-seed_00001"
FDATA = sub("replacement_.+-data_(.+)", "\\1", GEN)
FCOUNT = sub("replacement_.+-count_(.+)", "\\1", GEN)

###############
load(sprintf('sim-01/data/count_%s.RData', FCOUNT))
load(sprintf('sim-01/data/data_%s.RData', FDATA))
load(sprintf("sim-01/data/%s.RData", GEN))

paired.distance = function(H_orig, H_p){
  mean(apply(H_orig-H_p, 1, function(x) sqrt(sum(x^2))))
}

iNZ = which(rowSums(X == 0) != 0)
sapply(iNZ, function(i){
  pZ = X[i,] == 0
  nZ = sum(pZ)
  B = matrix(0, ncol = nZ, nrow = ncol(X))
  B[,1] = sbp_basis(cbind(2 * (pZ) - 1), silent = TRUE)
  if(nZ > 1){
    B[pZ,-1] = ilr_basis(nZ)
  }
  H = coordinates(P[i,], B)
  H.rpl = coordinates(P.rpl[i,], B)
  sqrt(sum((H-H.rpl)^2))
}) |> mean() -> paired.distance.in0


results = data.frame(
  gen = GEN,
  metric = 'Distance0',
  value = paired.distance.in0)

saveRDS(results, file = sprintf("sim-01/data/evaluate_paired.distance.in0-%s.RData", GEN))

# evaluate = function(H_gs, H_p){
#   c('paired.dist' = mean(apply(H_gs-H_p, 1, function(x) sqrt(sum(x^2)))),
#     'cov.frobenius' = norm(cov(H_gs) - cov(H_p), type = 'F'),
#     'stress' = sqrt(sum((as.matrix(dist(H_gs)) -  as.matrix(dist(H_p)))^2) / sum(as.matrix(dist(H_gs))^2)))
# }