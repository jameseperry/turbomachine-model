JULIA ?= julia --project=..
SCRIPTS := ../scripts
INTERPOLATIONS := bilinear bicubic

MEANLINE_TOMLS := $(wildcard compressor/meanline/*.toml)
DIM_MAP_TOMLS := $(wildcard compressor/map/dim/*/*.toml)
ND_MAP_TOMLS := $(wildcard compressor/map/nd/*/*.toml)

MEANLINE_ND_MAPS := $(foreach interp,$(INTERPOLATIONS),$(patsubst compressor/meanline/%.toml,compressor/map/nd/$(interp)/%_map_from_meanline.toml,$(MEANLINE_TOMLS)))
ALL_MAP_TOMLS := $(sort $(DIM_MAP_TOMLS) $(ND_MAP_TOMLS) $(MEANLINE_ND_MAPS))

PERFORMANCE_MAP_PLOTS := $(ALL_MAP_TOMLS:.toml=_performance_map.png)
OPERATING_SWEEP_PLOTS := $(ALL_MAP_TOMLS:.toml=_operating_sweep.png)

.PHONY: all sources plots sweeps clean

all: $(MEANLINE_ND_MAPS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)

sources:
	mkdir -p demo/compressor/map/dim/bilinear
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_tabulated demo/compressor/map/dim/bilinear/demo_tabulated_compressor_map.toml --interpolation bilinear --group compressor_map

	mkdir -p demo/compressor/map/dim/bicubic
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_tabulated demo/compressor/map/dim/bicubic/demo_tabulated_compressor_map.toml --interpolation bicubic --group compressor_map

	mkdir -p demo/compressor/map/nd/bilinear
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_nondimensional_tabulated demo/compressor/map/nd/bilinear/demo_tabulated_nondimensional_compressor_map.toml --interpolation bilinear --group compressor_map

	mkdir -p demo/compressor/map/nd/bicubic
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_nondimensional_tabulated demo/compressor/map/nd/bicubic/demo_tabulated_nondimensional_compressor_map.toml --interpolation bicubic --group compressor_map

	mkdir -p demo/compressor/meanline
	$(JULIA) $(SCRIPTS)/create_demo_object.jl compressor_meanline_model demo/compressor/meanline/demo_compressor_meanline_model.toml --group compressor_meanline_model


plots: $(PERFORMANCE_MAP_PLOTS)

sweeps: $(OPERATING_SWEEP_PLOTS)

compressor/map/nd/bilinear/%_map_from_meanline.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format nondimensional --interpolation bilinear

compressor/map/nd/bicubic/%_map_from_meanline.toml: compressor/meanline/%.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/generate_compressor_map_from_meanline_model.jl $< $@ --input-group compressor_meanline_model --output-group compressor_map --output-format nondimensional --interpolation bicubic

%_performance_map.png: %.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_performance_map.jl $< --kind compressor --output $@

%_operating_sweep.png: %.toml
	@mkdir -p $(@D)
	$(JULIA) $(SCRIPTS)/plot_compressor_operating_sweep.jl $< --target-pr 10.0 --branch all --csv $*_operating_sweep.csv --output $@

clean:
	rm -f $(MEANLINE_ND_MAPS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)
