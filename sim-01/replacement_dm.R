library(coda.count)


if(!exists("GEN")) GEN = "count_uniform-size_00030-data_parliament-seed_00001"

###############
load(sprintf("sim-01/data/%s.RData", GEN))


P.rpl = t(t(X) + as.vector(fit_dm(X)))
P.rpl = P.rpl / rowSums(P.rpl)

save(P.rpl, file = sprintf("sim-01/data/replacement_dm-%s.RData", GEN))
