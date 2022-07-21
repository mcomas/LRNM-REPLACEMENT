library(mvtnorm)
library(coda.base)

if(!exists("GEN")) GEN = 'dim_3-seed_00001'

GEN_PATTERN = "dim_([0-9]+)-seed_([0-9]+)"

DIM = as.integer(sub(GEN_PATTERN, "\\1", GEN))
SEED = as.integer(sub(GEN_PATTERN, "\\2", GEN))


###########################
set.seed(SEED)

rotation = random_rotation_matrix_incl_flip <- function(n){
  QR <- qr(matrix(rnorm(n^2), ncol=n))          # A = QR
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))  # ensure diag(R) > 0
  return(M)
}
d = DIM
M = rotation(d)

MU = rnorm(d)
SIGMA = M %*% diag(rlnorm(d, sdlog = 2)) %*% t(M)

H = rmvnorm(100, mean = MU, sigma = SIGMA) 

P = as.matrix(composition(H, ilr_basis(d+1)))

save(P, MU, SIGMA, file = sprintf("sim-01/data/data_lrnormal-%s.RData", GEN))
