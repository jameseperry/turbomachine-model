module Compressor

using ...Fluids

include("api.jl")
include("map_tabulated.jl")
include("map_tabulated_nd.jl")
include("map_analytic.jl")
include("meanline_model.jl")
include("map_conversion.jl")
include("map_io.jl")
include("spec.jl")
include("design.jl")
include("residuals.jl")
include("operating_point.jl")

export AbstractCompressorPerformanceMap
export TabulatedCompressorPerformanceMap
export NonDimensionalTabulatedCompressorPerformanceMap
export AnalyticCompressorPerformanceMap
export CompressorMeanlineModel
export AbstractRowAeroModel, RotorAeroModel, StatorAeroModel
export RowAeroInput, RowAeroOutput
export CompressorRow
export tabulate_compressor_meanline_model
export compressor_performance_map_from_stagnation
export performance_map_domain
export read_performance_map_toml
export to_nondimensional_tabulated_compressor_map, to_tabulated_compressor_map
export demo_tabulated_compressor_performance_map
export demo_nondimensional_tabulated_compressor_performance_map
export demo_analytic_compressor_performance_map
export demo_compressor_meanline_model
export demo_compressor_spec, demo_compressor_design
export write_toml, read_toml
export CompressorSpec, CompressorDesign
export compile_compressor_spec, compile_compressor_map
export compressor_residuals
export compressor_residuals_scaled
export compressor_residual_scales
export solve_compressor_operating_point
export compressor_pr_roots
export solve_compressor_operating_sweep

end # module Compressor
