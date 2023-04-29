library(fido)

if(!exists("GEN")) GEN = "count_uniform-size_00050-data_lrnormal-prop80-dim_10-seed_00001"

###############
load(sprintf("%s/data/%s.RData", SIM, GEN))
t0 = proc.time()

pib = to_clr(pibble(Y = t(X), X = rbind(rep(1, nrow(X)))))


P.rpl = prop.table(exp(t(apply(pib$Eta,1:2, mean))),1)


TIME = proc.time() - t0
save(P.rpl, TIME, file = sprintf("%s/data/replacement_lrnm-fido-%s.RData", SIM, GEN))
