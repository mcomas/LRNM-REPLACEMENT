library(mvtnorm)
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
M = rotation(5)

MU = rnorm(5)
SIGMA = t(M) %*% diag(rlnorm(5)) %*% M

set.seed(SEED)
H = mvtnorm::rmvnorm(100, mean = MU, sigma = SIGMA)
P = as.matrix(composition(H))

save(P, MU, SIGMA, file = sprintf("sim-01/data/data_mvtnorm-%s.RData", GEN))
