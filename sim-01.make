L_data = lrnormal
L_seed = $(shell seq 1 5)
L_dim = 3 5 10 20
L_count =  uniform dim
L_uniform_size = $(shell seq 30 40 110)
L_dim_size = 50 100
L_replacement = lrnm-montecarlo lrnb-cond-1-hermite-new dm 
L_evaluate = stress
# frobenius stress1
# paired.distance paired.distance.in0 

####
####
L_seed_s = $(foreach seed,$(L_seed), $(shell printf '%05d' $(seed)))

L_CODA = $(foreach data,$(L_data),$(foreach dim,$(L_dim),$(foreach seed,$(L_seed_s),$(shell printf 'data_%s-dim_%d-seed_%s' $(data) $(dim) $(seed)))))
CODA = $(foreach coda,$(L_CODA),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_uniform_size_s = $(foreach size,$(L_uniform_size), $(shell printf '%05d' $(size)))
L_dim_size_s = $(foreach size,$(L_dim_size), $(shell printf '%05d' $(size)))
L_CODA_COUNT = $(foreach size,$(L_uniform_size_s),$(foreach coda,$(L_CODA),$(shell printf 'count_uniform-size_%s-%s' $(size) $(coda)))) \
               $(foreach size,$(L_dim_size_s),$(foreach coda,$(L_CODA),$(shell printf 'count_dim-size_%s-%s' $(size) $(coda))))
CODA_COUNT = $(foreach coda,$(L_CODA_COUNT),$(shell printf 'sim-01/data/%s.RData' $(coda)))

L_REPLACEMENT = $(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(shell printf 'replacement_%s-%s' $(replacement) $(coda_count))))
REPLACEMENT = $(foreach replacement,$(L_REPLACEMENT),$(shell printf 'sim-01/data/%s.RData' $(replacement)))

L_EVALUATE = $(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(shell printf 'evaluate_%s-%s' $(evaluate) $(replacement))))
EVALUATE = $(foreach evaluate,$(L_EVALUATE),$(shell printf 'sim-01/data/%s.rds' $(evaluate)))


all : $(CODA) $(CODA_COUNT) $(REPLACEMENT) $(EVALUATE) # sim-01/datasets-summary.RData
data : $(CODA) $(CODA_COUNT)

sim-01/datasets-summary.RData : sim-01/datasets-summary.R $(CODA_COUNT)
	Rscript $<

define DATA_RULE
sim-01/data/data_$(data)-dim_$(dim)-seed_$(seed).RData : sim-01/data_$(data).R
	Rscript -e 'GEN = "dim_$(dim)-seed_$(seed)"; source("$$<")'
endef
$(foreach data,$(L_data),$(foreach seed,$(L_seed_s),$(foreach dim,$(L_dim),$(eval $(DATA_RULE)))))

define UNIFORM_COUNT_RULE
sim-01/data/count_uniform-size_$(size)-$(coda).RData : sim-01/count_uniform.R sim-01/data/$(coda).RData
	Rscript -e 'GEN = "size_$(size)-$(coda)"; source("$$<")'
endef
$(foreach size,$(L_uniform_size_s),$(foreach coda,$(L_CODA),$(eval $(UNIFORM_COUNT_RULE))))

define DIM_COUNT_RULE
sim-01/data/count_dim-size_$(size)-$(coda).RData : sim-01/count_dim.R sim-01/data/$(coda).RData
	Rscript -e 'GEN = "size_$(size)-$(coda)"; source("$$<")'
endef
$(foreach size,$(L_dim_size_s),$(foreach coda,$(L_CODA),$(eval $(DIM_COUNT_RULE))))

define REPLACEMENT_RULE
sim-01/data/replacement_$(replacement)-$(coda_count).RData : sim-01/replacement_$(replacement).R sim-01/data/$(coda_count).RData
	Rscript -e 'GEN = "$(coda_count)"; source("$$<")' > sim-01/data/replacement_$(replacement)-$(coda_count).log
endef
$(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(eval $(REPLACEMENT_RULE))))

define EVALUATE_RULE
sim-01/data/evaluate_$(evaluate)-$(replacement).rds : sim-01/evaluate_$(evaluate).R sim-01/data/$(replacement).RData
	Rscript -e 'GEN = "$(replacement)"; source("$$<")'
endef
$(foreach evaluate,$(L_evaluate),$(foreach replacement,$(L_REPLACEMENT),$(eval $(EVALUATE_RULE))))

clean :
	rm -f sim-01/data/*.RData
	rm -f sim-01/datasets-summary.RData

