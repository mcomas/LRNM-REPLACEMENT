library(coda.count)

if(!exists("SIM")) SIM = 'sim-01a'
if(!exists("GEN")) GEN = "size_00030-data_lrnormal-dim_3-seed_00001"

GEN_PATTERN = "size_([0-9]+)-data_(.+)-dim_(.+)-seed_([0-9]+)"

SIZE = as.integer(sub(GEN_PATTERN, "\\1", GEN))
SEED = as.integer(sub(GEN_PATTERN, "\\4", GEN))

FDATA = sub("size_([0-9]+)-(.+)", "\\2", GEN)

############

load(sprintf('%s/data/%s.RData', SIM, FDATA))

set.seed(SEED)
NO_ZEROS = TRUE
# while(NO_ZEROS){
  X = rmultinomial(size = SIZE, p = P)
  NO_ZEROS = sum(X == 0) < 5 | min(colSums(X != 0)) == 0
# }

save(X, file = sprintf("%s/data/count_uniform-%s.RData", SIM, GEN))
