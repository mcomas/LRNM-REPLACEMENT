SIM = sim-01a sim-01b sim-02a
SIM_BUILDS = $(foreach sim,$(SIM), $(shell printf '.%s_build' $(sim)))

BUILDS = $(SIM_BUILDS)

all : $(BUILDS)

.%_build : FORCE 
	$(MAKE) -f $*.make
	date > $@

FORCE: