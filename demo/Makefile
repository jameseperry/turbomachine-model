JULIA ?= julia --project=..
SCRIPTS := ../scripts

INTERPOLATION := bicubic

SWEEP_N_OMEGA ?= 100
SWEEP_N_PT_IN ?= 100
DIAG_PT_IN ?= 101325
DIAG_TT_IN ?= 288.15
COMPRESSOR_DIAG_OMEGA ?= 800
COMPRESSOR_DIAG_MDOT ?= 16
TURBINE_DIAG_OMEGA ?= 600
TURBINE_DIAG_MDOT ?= 8

USER_AXIAL_TOMLS  := $(wildcard axial/*/*.toml)
USER_MAP_TOMLS := $(wildcard map/*/*.toml)

GENERATED_AXIAL_TOMLS := axial/compressor/demo.toml axial/turbine/demo.toml
GENERATED_MAP_TOMLS := $(patsubst axial/%.toml,map/%.toml,$(GENERATED_AXIAL_TOMLS))

ALL_AXIAL_TOMLS := $(sort $(USER_AXIAL_TOMLS) $(GENERATED_AXIAL_TOMLS))
ALL_MAP_TOMLS := $(sort $(USER_MAP_TOMLS) $(GENERATED_MAP_TOMLS))

PERFORMANCE_MAP_PLOTS := $(ALL_MAP_TOMLS:.toml=.performance_map.png)
OPERATING_SWEEP_PLOTS := $(ALL_MAP_TOMLS:.toml=.operating_sweep.png)
OPERATING_SWEEP_CSV := $(ALL_MAP_TOMLS:.toml=.operating_sweep.csv)
AXIAL_DIAGNOSTIC_TXT := $(ALL_AXIAL_TOMLS:.toml=.operating_point.txt)

.PHONY: all plots sweeps diagnostics clean
.SECONDARY: $(GENERATED_MAP_TOMLS)

all: $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS) $(AXIAL_DIAGNOSTIC_TXT)

performance-plots: $(PERFORMANCE_MAP_PLOTS)

sweep-plots: $(OPERATING_SWEEP_PLOTS)

diagnostics: $(AXIAL_DIAGNOSTIC_TXT)

axial/compressor/demo.toml: $(SCRIPTS)/create_demo_object.jl
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/create_demo_object.jl axial_compressor_model $@

axial/turbine/demo.toml: $(SCRIPTS)/create_demo_object.jl
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/create_demo_object.jl axial_turbine_model $@

axial/compressor/%.operating_point.txt: axial/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/diagnose_axial_operating_point.jl $< \
		--pt-in $(DIAG_PT_IN) \
		--tt-in $(DIAG_TT_IN) \
		--omega $(COMPRESSOR_DIAG_OMEGA) \
		--mdot $(COMPRESSOR_DIAG_MDOT) \
		--output $@

axial/turbine/%.operating_point.txt: axial/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/diagnose_axial_operating_point.jl $< \
		--pt-in $(DIAG_PT_IN) \
		--tt-in $(DIAG_TT_IN) \
		--omega $(TURBINE_DIAG_OMEGA) \
		--mdot $(TURBINE_DIAG_MDOT) \
		--output $@

map/compressor/%.toml: axial/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_performance_map_from_axial_machine_model.jl $< $@ \
		--interpolation $(INTERPOLATION)

map/turbine/%.toml: axial/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_performance_map_from_axial_machine_model.jl $< $@ \
		--interpolation $(INTERPOLATION) \
		--n-speed 100 \
		--n-flow 100

map/compressor/%.performance_map.png: map/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --output $@

map/turbine/%.performance_map.png: map/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --output $@

map/compressor/%.operating_sweep.png: map/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_operating_sweep.jl $< \
		--machine-kind compressor \
		--track-all-branches \
		--n-omega $(SWEEP_N_OMEGA) \
		--n-pr $(SWEEP_N_PT_IN) \
		--pr-min 0.5 \
		--pr-max 10.0 \
		--csv $(@:.operating_sweep.png=.operating_sweep.csv) \
		--output $@

map/turbine/%.operating_sweep.png: map/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_operating_sweep.jl $< \
		--machine-kind turbine \
		--track-all-branches \
		--n-omega $(SWEEP_N_OMEGA) \
		--n-pr $(SWEEP_N_PT_IN) \
		--pr-min 0.3 \
		--pr-max 1.0 \
		--csv $(@:.operating_sweep.png=.operating_sweep.csv) \
		--output $@

clean:
	rm -f $(GENERATED_AXIAL_TOMLS) $(GENERATED_MAP_TOMLS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS) $(OPERATING_SWEEP_CSV) $(AXIAL_DIAGNOSTIC_TXT)
