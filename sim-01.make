OUTPUT = sim-01

L_data = mvtnorm
L_seed = 1
L_count =  uniform
L_size = 100 10
L_replacement = gbm

L_CODA = $(foreach data,$(L_data),$(foreach seed,$(L_seed),$(shell printf 'data_%s-seed_%05d' $(data) $(seed))))
CODA = $(foreach coda,$(L_CODA),$(shell printf '$(OUTPUT)/%s.RData' $(coda)))

L_CODA_COUNT = $(foreach count,$(L_count),$(foreach size,$(L_size),$(foreach coda,$(L_CODA),$(shell printf 'count_%s-size_%05d-%s' $(count) $(size) $(coda)))))
CODA_COUNT = $(foreach coda,$(L_CODA_COUNT),$(shell printf '$(OUTPUT)/%s.RData' $(coda)))

L_REPLACEMENT = $(foreach replacement,$(L_replacement),$(foreach coda_count,$(L_CODA_COUNT),$(shell printf 'replacement_%s-%s' $(replacement) $(coda_count))))
REPLACEMENT = $(foreach replacement,$(L_REPLACEMENT),$(shell printf '$(OUTPUT)/%s.RData' $(replacement)))


all : $(CODA) $(CODA_COUNT) $(REPLACEMENT)

$(OUTPUT)/data_mvtnorm-%.RData : $(OUTPUT)/data_mvtnorm.R 
	Rscript -e 'GEN = "$*"; source("$<")'

$(OUTPUT)/count_uniform-%.RData : $(OUTPUT)/count_uniform.R $(CODA)
	Rscript -e 'GEN = "$*"; source("$<")'

$(OUTPUT)/replacement_gbm-%.RData : $(OUTPUT)/replacement_gbm.R $(CODA_COUNT)
	Rscript -e 'GEN = "$*"; source("$<")'
	
clean :
	rm $(OUTPUT)/*.RData
	