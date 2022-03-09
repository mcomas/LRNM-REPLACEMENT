library(MixSim)
library(coda.base)

if(!exists("GEN")) GEN = 'seed_0001'

GEN_PATTERN = "seed_([0-9]+)"
SEED = as.integer(sub(GEN_PATTERN, "\\1", GEN))

###########################
set.seed(SEED)

n = 100
D = 4
H = matrix(rnorm(D*n), ncol = D)
H = runif(n, 0.75,1) * H/sqrt(rowSums(H^2))

dir = rnorm(D)
dir = dir / sqrt(sum(dir^2))
H = t(t(H) + 5 * dir)

P = composition(H)

pars = list(dir = dir)

save(P, pars, file = sprintf("sim-01/data/data_circle-%s.RData", GEN))
