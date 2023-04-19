library(coda.count)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "size_00030-data_lrnormal-prop80-dim_3-seed_00001"

GEN_PATTERN = "size_([0-9]+)-data_(.+)-dim_(.+)-seed_([0-9]+)"

SIZE = as.integer(sub(GEN_PATTERN, "\\1", GEN))
DIM = as.integer(sub(GEN_PATTERN, "\\3", GEN))
SEED = as.integer(sub(GEN_PATTERN, "\\4", GEN))
FDATA = sub("size_([0-9]+)-(.+)", "\\2", GEN)

############

load(sprintf('%s/data/%s.RData', SIM, FDATA))

set.seed(SEED)
NO_ZEROS = TRUE
while(NO_ZEROS){
  X = rmultinomial(size = DIM * SIZE, p = P)
  # Force that at least one observation has non-zero in each part.
  NO_ZEROS = min(colSums(X != 0)) <= 1 
}

save(X, file = sprintf("%s/data/count_dim-%s.RData", SIM, GEN))
