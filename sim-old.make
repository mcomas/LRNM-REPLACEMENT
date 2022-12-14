L_data = circle mixture parliament 
# mvtnorm iris
L_seed = $(shell seq 1 5)
L_count =  uniform
L_size = $(shell seq 10 20 110)
L_replacement = lrnm-cond-hermite2 em-normal czm lrnm-cond-hermite  lrnb-cond-equal-hermite lrnb-cond-1by1-hermite lrnm zinflrnm dm
# mixture-lrn-laplace one-mixture-lrn-laplace em-lrn-laplace em-normal 
# lrnm-laplace gbm dir-mixture-lrn-laplace
L_evaluate = paired.distance paired.distance.in0
# frobenius stress1


####
####
L_seed_s = $(foreach seed,$(L_seed), $(shell printf '%05d' $(seed)))

L_CODA = $(foreach data,$(L_data),$(foreach seed,$(L_seed_s),$(shell printf 'data_%s-seed_%s' $(data) $(seed))))
CODA = $(foreach coda,$(L_CODA),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_size_s = $(foreach size,$(L_size), $(shell printf '%05d' $(size)))
L_CODA_COUNT = $(foreach count,$(L_count),$(foreach size,$(L_size_s),$(foreach coda,$(L_CODA),$(shell printf 'count_%s-size_%s-%s' $(count) $(size) $(coda)))))
CODA_COUNT = $(foreach coda,$(L_CODA_COUNT),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_REPLACEMENT = $(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(shell printf 'replacement_%s-%s' $(replacement) $(coda_count))))
REPLACEMENT = $(foreach replacement,$(L_REPLACEMENT),$(shell printf 'sim-01/data/%s.RData' $(replacement)))

L_EVALUATE = $(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(shell printf 'evaluate_%s-%s' $(evaluate) $(replacement))))
EVALUATE = $(foreach evaluate,$(L_EVALUATE),$(shell printf 'sim-01/data/%s.RData' $(evaluate)))


all : $(CODA) $(CODA_COUNT) $(REPLACEMENT) $(EVALUATE) sim-01/datasets-summary.RData
data : $(CODA) $(CODA_COUNT)

sim-01/datasets-summary.RData : sim-01/datasets-summary.R $(CODA_COUNT)
	Rscript $<

define DATA_RULE
sim-01/data/data_$(data)-seed_$(seed).RData : sim-01/data_$(data).R
	Rscript -e 'GEN = "seed_$(seed)"; source("$$<")'
endef
$(foreach data,$(L_data),$(foreach seed,$(L_seed_s),$(eval $(DATA_RULE))))

define COUNT_RULE
sim-01/data/count_$(count)-size_$(size)-$(coda).RData : sim-01/count_$(count).R sim-01/data/$(coda).RData
	Rscript -e 'GEN = "size_$(size)-$(coda)"; source("$$<")'
endef
$(foreach count,$(L_count),$(foreach size,$(L_size_s),$(foreach coda,$(L_CODA),$(eval $(COUNT_RULE)))))

define REPLACEMENT_RULE
sim-01/data/replacement_$(replacement)-$(coda_count).RData : sim-01/replacement_$(replacement).R sim-01/data/$(coda_count).RData
	Rscript -e 'GEN = "$(coda_count)"; source("$$<")' > sim-01/data/replacement_$(replacement)-$(coda_count).log
endef
$(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(eval $(REPLACEMENT_RULE))))

define EVALUATE_RULE
sim-01/data/evaluate_$(evaluate)-$(replacement).RData : sim-01/evaluate_$(evaluate).R sim-01/data/$(replacement).RData
	Rscript -e 'GEN = "$(replacement)"; source("$$<")'
endef
$(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(eval $(EVALUATE_RULE))))

clean :
	rm -f sim-01/data/*.RData
	rm -f sim-01/datasets-summary.RData

