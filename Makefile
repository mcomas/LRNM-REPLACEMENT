SIM = sim-01a sim-01b #sim-02a sim-03a
SIM_BUILDS = $(foreach sim,$(SIM), $(shell printf '.%s_build' $(sim)))

BUILDS = $(SIM_BUILDS)

all : $(BUILDS)

.sim-01a_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-01a DATA=lrnormal-prop80 COUNT=uniform
	date > $@

.sim-01b_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-01b DATA=lrnormal-prop80 COUNT=dim
	date > $@

.sim-02a_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-02a DATA=lrskew COUNT=dim
	date > $@

.sim-03a_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-03a DATA=mixt2 COUNT=dim
	date > $@

.%_build : FORCE 
	$(MAKE) -f $*.make
	date > $@

FORCE: