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

H.iris = coordinates(iris[,1:4])
M = rotation(3)
P = composition(as.matrix(H.iris) %*% M)

save(P, M, file = sprintf("sim-01/data_iris-%s.RData", GEN))
