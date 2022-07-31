library(sn)
library(coda.base)
if(!exists("SIM")) SIM = 'sim-02a'
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

B = ilr_basis(d+1)

REPEAT = TRUE
while(REPEAT){
  XI = rnorm(d)
  PXI = composition(XI, B)
  REPEAT = min(PXI) > 50/(200*d)
}
REPEAT = TRUE
while(REPEAT){
  EIGVALS = rlnorm(d)
  REPEAT = max(EIGVALS) / sum(EIGVALS) < 0.8
}
OMEGA = M %*% diag(EIGVALS) %*% t(M)
OMEGA = 0.5 * OMEGA + 0.5 * t(OMEGA)
ALPHA = rnorm(d) * 1000 #runif(1, 10, 1000)

H = rmsn(100, xi = XI, Omega = OMEGA, alpha = ALPHA)

P = as.matrix(composition(H, ilr_basis(d+1)))

save(P, XI, OMEGA, ALPHA, file = sprintf("%s/data/data_lrskew-%s.RData", SIM, GEN))
