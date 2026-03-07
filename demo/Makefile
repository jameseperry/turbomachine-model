JULIA ?= julia --project=..
SCRIPTS := ../scripts

INTERPOLATION := bicubic

TURBINE_SWEEP_N_OMEGA ?= 80
TURBINE_SWEEP_N_PT_IN ?= 60
TURBINE_SWEEP_PT_OUT ?= 101325.0

USER_AXIAL_TOMLS  := $(wildcard axial/*/*.toml)
USER_MAP_TOMLS := $(wildcard map/*/*.toml)

GENERATED_AXIAL_TOMLS := axial/compressor/demo.toml axial/turbine/demo.toml
GENERATED_MAP_TOMLS := $(patsubst axial/%.toml,map/%.toml,$(GENERATED_AXIAL_TOMLS))

ALL_MAP_TOMLS := $(sort $(USER_MAP_TOMLS) $(GENERATED_MAP_TOMLS))

PERFORMANCE_MAP_PLOTS := $(ALL_MAP_TOMLS:.toml=.performance_map.png)
OPERATING_SWEEP_PLOTS := $(ALL_MAP_TOMLS:.toml=.operating_sweep.png)
OPERATING_SWEEP_CSV := $(ALL_MAP_TOMLS:.toml=.operating_sweep.csv)

.PHONY: all plots sweeps clean
.SECONDARY: $(GENERATED_MAP_TOMLS)

all: $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)

performance-plots: $(PERFORMANCE_MAP_PLOTS)

sweep-plots: $(OPERATING_SWEEP_PLOTS)

axial/compressor/demo.toml: $(SCRIPTS)/create_demo_object.jl
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_meanline_model $@

axial/turbine/demo.toml: $(SCRIPTS)/create_demo_object.jl
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_meanline_turbine_model $@

map/compressor/%.toml: axial/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ \
		--output-format dimensional \
		--interpolation $(INTERPOLATION)

map/turbine/%.toml: axial/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_turbine_map_from_meanline_model.jl $< $@ \
		--interpolation $(INTERPOLATION) \
		--sampling-strategy scan

map/compressor/%.performance_map.png: map/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --output $@

map/turbine/%.performance_map.png: map/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --output $@

map/compressor/%.operating_sweep.png: map/compressor/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_compressor_operating_sweep.jl $< \
		--target-pr 10.0 \
		--allow-windmill \
		--branch all \
		--csv $(@:.operating_sweep.png=.operating_sweep.csv) \
		--output $@

map/turbine/%.operating_sweep.png: map/turbine/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_turbine_operating_sweep.jl $< \
		--mode fixed_pt_out_grid \
		--n-points $(TURBINE_SWEEP_N_OMEGA) \
		--n-pt-in $(TURBINE_SWEEP_N_PT_IN) \
		--pt-out $(TURBINE_SWEEP_PT_OUT) \
		--csv $(@:.operating_sweep.png=.operating_sweep.csv) \
		--output $@

clean:
	rm -f $(GENERATED_AXIAL_TOMLS) $(GENERATED_MAP_TOMLS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS) $(OPERATING_SWEEP_CSV)
