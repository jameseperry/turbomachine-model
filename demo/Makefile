JULIA ?= julia --project=..
SCRIPTS := ../scripts
INTERPOLATIONS := bilinear bicubic

MEANLINE_TOMLS := $(wildcard compressor/meanline/*.toml)
MEANLINE_TOMLS += compressor/meanline/demo_compressor_meanline_model.toml

DIM_MAP_TOMLS := $(wildcard compressor/map/dim/*/*.toml)
ND_MAP_TOMLS := $(wildcard compressor/map/nd/*/*.toml)

MEANLINE_ND_MAPS := $(foreach interp,$(INTERPOLATIONS),$(patsubst compressor/meanline/%.toml,compressor/map/nd/$(interp)/%.toml,$(MEANLINE_TOMLS)))
MEANLINE_DIM_MAPS := $(foreach interp,$(INTERPOLATIONS),$(patsubst compressor/meanline/%.toml,compressor/map/dim/$(interp)/%.toml,$(MEANLINE_TOMLS)))

ALL_DIM_MAP_TOMLS := $(sort $(DIM_MAP_TOMLS) $(MEANLINE_DIM_MAPS))
ALL_ND_MAP_TOMLS := $(sort $(ND_MAP_TOMLS) $(MEANLINE_ND_MAPS))

ALL_MAP_TOMLS := $(sort $(ALL_DIM_MAP_TOMLS) $(ALL_ND_MAP_TOMLS))

PERFORMANCE_MAP_PLOTS := $(ALL_MAP_TOMLS:.toml=_performance_map.png)
OPERATING_SWEEP_PLOTS := $(ALL_MAP_TOMLS:.toml=_operating_sweep.png)

.PHONY: all plots sweeps clean

all: $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)

plots: $(PERFORMANCE_MAP_PLOTS)

sweeps: $(OPERATING_SWEEP_PLOTS)

compressor/meanline/demo_compressor_meanline_model.toml:
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_meanline_model $@

compressor/map/nd/bilinear/%.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format nondimensional --interpolation bilinear

compressor/map/nd/bicubic/%.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format nondimensional --interpolation bicubic

compressor/map/dim/bilinear/%.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format dimensional --interpolation bilinear

compressor/map/dim/bicubic/%.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format dimensional --interpolation bicubic

%_performance_map.png: %.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --kind compressor --output $@

%_operating_sweep.png: %.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_compressor_operating_sweep.jl $< --target-pr 10.0 --branch all --csv $*_operating_sweep.csv --output $@

clean:
	rm -f $(MEANLINE_ND_MAPS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)
