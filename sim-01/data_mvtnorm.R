library(mvtnorm)
library(coda.base)

if(!exists("GEN")) GEN = 'seed_0001'

GEN_PATTERN = "seed_([0-9]+)"
SEED = as.integer(sub(GEN_PATTERN, "\\1", GEN))

###########################
set.seed(SEED)

H = mvtnorm::rmvnorm(100, mean = rep(0, 5))
P = as.matrix(composition(H))

save(P, file = sprintf("sim-01/data_mvtnorm-%s.RData", GEN))
