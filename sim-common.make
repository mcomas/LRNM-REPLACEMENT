SIM = sim-01a
DATA = lrnormal-prop80
COUNT = dim
######
L_data = $(DATA)
L_seed = $(shell seq 1 5)
L_dim = 3 5 10 15 20
# 15 20
L_count =  $(COUNT)
L_count_size = $(shell seq 50 50 200)
L_replacement = dm lrnm-montecarlo lrnb-cond-1-hermite-new lrnm-cond-montecarlo gbm
L_evaluate = stress time


####
####
L_seed_s = $(foreach seed,$(L_seed), $(shell printf '%05d' $(seed)))

L_CODA = $(foreach data,$(L_data),$(foreach dim,$(L_dim),$(foreach seed,$(L_seed_s),$(shell printf 'data_%s-dim_%d-seed_%s' $(data) $(dim) $(seed)))))
CODA = $(foreach coda,$(L_CODA),$(shell printf '$(SIM)/data/%s.RData' $(coda)))

L_count_size_s = $(foreach size,$(L_count_size), $(shell printf '%05d' $(size)))
L_CODA_COUNT = $(foreach count,$(L_count),$(foreach size,$(L_count_size_s),$(foreach coda,$(L_CODA),$(shell printf 'count_%s-size_%s-%s' $(count) $(size) $(coda)))))
CODA_COUNT = $(foreach coda,$(L_CODA_COUNT),$(shell printf '$(SIM)/data/%s.RData' $(coda)))

L_REPLACEMENT = $(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(shell printf 'replacement_%s-%s' $(replacement) $(coda_count))))
REPLACEMENT = $(foreach replacement,$(L_REPLACEMENT),$(shell printf '$(SIM)/data/%s.RData' $(replacement)))

L_EVALUATE = $(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(shell printf 'evaluate_%s-%s' $(evaluate) $(replacement))))
EVALUATE = $(foreach evaluate,$(L_EVALUATE),$(shell printf '$(SIM)/data/%s.rds' $(evaluate)))

FIGURES = overleaf/$(SIM)-fig01.pdf overleaf/$(SIM)-fig02.pdf 

all : $(CODA) $(CODA_COUNT) $(REPLACEMENT) $(EVALUATE) $(FIGURES)# $(SIM)/datasets-summary.RData
data : $(CODA) $(CODA_COUNT)

$(SIM)/datasets-summary.RData : $(SIM)/datasets-summary.R $(CODA_COUNT)
	Rscript $<

define DATA_RULE
$(SIM)/data/data_$(data)-dim_$(dim)-seed_$(seed).RData : sim-common/data_$(data).R
	Rscript -e 'SIM = "$(SIM)"; GEN = "dim_$(dim)-seed_$(seed)"; source("$$<")'
endef
$(foreach data,$(L_data),$(foreach seed,$(L_seed_s),$(foreach dim,$(L_dim),$(eval $(DATA_RULE)))))

define COUNT_RULE
$(SIM)/data/count_$(count)-size_$(size)-$(coda).RData : sim-common/count_$(count).R $(SIM)/data/$(coda).RData
	Rscript -e 'SIM = "$(SIM)"; GEN = "size_$(size)-$(coda)"; source("$$<")'
endef
$(foreach count,$(L_count),$(foreach size,$(L_count_size_s),$(foreach coda,$(L_CODA),$(eval $(COUNT_RULE)))))

define REPLACEMENT_RULE
$(SIM)/data/replacement_$(replacement)-$(coda_count).RData : sim-common/replacement_$(replacement).R $(SIM)/data/$(coda_count).RData
	Rscript -e 'SIM = "$(SIM)"; GEN = "$(coda_count)"; source("$$<")' > $(SIM)/data/replacement_$(replacement)-$(coda_count).log
endef
$(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(eval $(REPLACEMENT_RULE))))

define EVALUATE_RULE
$(SIM)/data/evaluate_$(evaluate)-$(replacement).rds : sim-common/evaluate_$(evaluate).R $(SIM)/data/$(replacement).RData
	Rscript -e 'SIM = "$(SIM)"; GEN = "$(replacement)"; source("$$<")'
endef
$(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(eval $(EVALUATE_RULE))))

overleaf/$(SIM)-fig%.pdf : sim-common/figure-%.R $(EVALUATE)
	Rscript -e 'SIM = "$(SIM)"; source("$<")'

clean :
	rm -f $(SIM)/data/*.RData
	rm -f $(SIM)/datasets-summary.RData

