library(MixSim)
library(coda.base)

if(!exists("GEN")) GEN = 'seed_0001'

GEN_PATTERN = "seed_([0-9]+)"
SEED = as.integer(sub(GEN_PATTERN, "\\1", GEN))

###########################
set.seed(SEED)

rotation = random_rotation_matrix_incl_flip <- function(n){
  QR <- qr(matrix(rnorm(n^2), ncol=n))          # A = QR
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))  # ensure diag(R) > 0
  return(M)
}

pars = MixSim(MaxOmega = 0.05, K = 4, p = 10)
H = do.call(simdataset, c(list(n = 50), pars[c('Pi', 'Mu', 'S')]))$X
P = composition(H)

save(P, pars, file = sprintf("sim-01/data/data_mixture-%s.RData", GEN))
