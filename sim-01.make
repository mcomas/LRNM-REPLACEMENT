L_data = mvtnorm iris mixture parliament
L_seed = $(shell seq 1 2)
L_count =  uniform
L_size = $(shell seq 10 10 50)
L_replacement = gbm czm
L_evaluate = stress frobenius paired.distance

L_CODA = $(foreach data,$(L_data),$(foreach seed,$(L_seed),$(shell printf 'data_%s-seed_%05d' $(data) $(seed))))
CODA = $(foreach coda,$(L_CODA),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_CODA_COUNT = $(foreach count,$(L_count),$(foreach size,$(L_size),$(foreach coda,$(L_CODA),$(shell printf 'count_%s-size_%05d-%s' $(count) $(size) $(coda)))))
CODA_COUNT = $(foreach coda,$(L_CODA_COUNT),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_REPLACEMENT = $(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(shell printf 'replacement_%s-%s' $(replacement) $(coda_count))))
REPLACEMENT = $(foreach replacement,$(L_REPLACEMENT),$(shell printf 'sim-01/data/%s.RData' $(replacement)))

L_EVALUATE = $(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(shell printf 'evaluate_%s-%s' $(evaluate) $(replacement))))
EVALUATE = $(foreach evaluate,$(L_EVALUATE),$(shell printf 'sim-01/data/%s.RData' $(evaluate)))

all : $(CODA) $(CODA_COUNT) $(REPLACEMENT) $(EVALUATE)

sim-01/data/data_mvtnorm-%.RData : sim-01/data_mvtnorm.R 
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/data_iris-%.RData : sim-01/data_iris.R 
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/data_mixture-%.RData : sim-01/data_mixture.R 
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/data_parliament-%.RData : sim-01/data_parliament.R 
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/count_uniform-%.RData : sim-01/count_uniform.R $(CODA)
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/replacement_gbm-%.RData : sim-01/replacement_gbm.R $(CODA_COUNT)
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/replacement_czm-%.RData : sim-01/replacement_czm.R $(CODA_COUNT)
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/evaluate_stress-%.RData : sim-01/evaluate_stress.R $(REPLACEMENT)
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/evaluate_frobenius-%.RData : sim-01/evaluate_frobenius.R $(REPLACEMENT)
	Rscript -e 'GEN = "$*"; source("$<")'

sim-01/data/evaluate_paired.distance-%.RData : sim-01/evaluate_paired.distance.R $(REPLACEMENT)
	Rscript -e 'GEN = "$*"; source("$<")'

clean :
	rm sim-01/data/*.RData
