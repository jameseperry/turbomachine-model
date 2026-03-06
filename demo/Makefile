JULIA ?= julia --project=..
SCRIPTS := ../scripts
MAP_VARIANTS := dim/bilinear dim/bicubic nd/bilinear nd/bicubic

USER_MEANLINE_TOMLS := $(wildcard compressor/meanline/*.toml)
USER_MAP_TOMLS := $(wildcard compressor/map/*/*/*.toml)

GENERATED_MEANLINE_TOMLS := compressor/meanline/demo_compressor_meanline_model.toml
GENERATED_MAP_TOMLS := $(foreach variant,$(MAP_VARIANTS),$(patsubst compressor/meanline/%.toml,compressor/map/$(variant)/%.toml,$(GENERATED_MEANLINE_TOMLS)))

ALL_MAP_TOMLS := $(sort $(USER_MAP_TOMLS) $(GENERATED_MAP_TOMLS))

PERFORMANCE_MAP_PLOTS := $(ALL_MAP_TOMLS:.toml=_performance_map.png)
OPERATING_SWEEP_PLOTS := $(ALL_MAP_TOMLS:.toml=_operating_sweep.png)
OPERATING_SWEEP_CSV := $(ALL_MAP_TOMLS:.toml=_operating_sweep.csv)

.PHONY: all plots sweeps clean

all: $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS)

plots: $(PERFORMANCE_MAP_PLOTS)

sweeps: $(OPERATING_SWEEP_PLOTS)

compressor/meanline/demo_compressor_meanline_model.toml: $(SCRIPTS)/create_demo_object.jl
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
	rm -f $(GENERATED_MEANLINE_TOMLS) $(GENERATED_MAP_TOMLS) $(PERFORMANCE_MAP_PLOTS) $(OPERATING_SWEEP_PLOTS) $(OPERATING_SWEEP_CSV)
