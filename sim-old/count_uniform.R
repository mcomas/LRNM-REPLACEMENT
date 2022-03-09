library(coda.count)

if(!exists("GEN")) GEN = "size_00100-data_mvtnorm-seed_00001"

GEN_PATTERN = "size_([0-9]+)-data_(.+)-seed_([0-9]+)"
SIZE = as.integer(sub(GEN_PATTERN, "\\1", GEN))
SEED = as.integer(sub(GEN_PATTERN, "\\3", GEN))
FDATA = sub("size_([0-9]+)-(.+)", "\\2", GEN)

############

load(sprintf('sim-01/data/%s.RData', FDATA))

set.seed(SEED)
NO_ZEROS = TRUE
while(NO_ZEROS){
  X = rmultinomial(size = SIZE, p = P)
  NO_ZEROS = sum(X == 0) < 5 | min(colSums(X != 0)) == 0
}

save(X, file = sprintf("sim-01/data/count_uniform-%s.RData", GEN))
