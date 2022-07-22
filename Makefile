SIM = sim-01a sim-01b sim-02a sim-03a
SIM_BUILDS = $(foreach sim,$(SIM), $(shell printf '.%s_build' $(sim)))

BUILDS = $(SIM_BUILDS)

all : $(BUILDS)

.sim-01a_build : FORCE
	$(MAKE) -f sim-common.make -e COUNT=uniform
	$(MAKE) overleaf/sim-01a-fig01b.pdf -f sim-common.make -e COUNT=uniform

.sim-01b_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-01b
	$(MAKE) overleaf/sim-01b-fig01b.pdf -f sim-common.make -e SIM=sim-01b

.sim-02a_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-02a DATA=lrskew

.sim-03a_build : FORCE
	$(MAKE) -f sim-common.make -e SIM=sim-03a DATA=mixt2

.%_build : FORCE 
	$(MAKE) -f $*.make
	date > $@

FORCE: