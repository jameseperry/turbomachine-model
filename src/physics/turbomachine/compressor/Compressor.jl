module Compressor

using ...Fluids

include("map_tabulated.jl")
include("map_analytic.jl")
include("spec.jl")
include("design.jl")
include("residuals.jl")

export AbstractCompressorPerformanceMap
export TabulatedCompressorPerformanceMap
export AnalyticCompressorPerformanceMap
export corrected_speed, corrected_flow
export compressor_performance_map, compressor_performance_map_from_stagnation
export demo_compressor_performance_map
export write_hdf5, read_hdf5
export write_toml, read_toml
export CompressorSpec, CompressorDesign
export compile_compressor_spec, compile_compressor_map
export mdot_surge, mdot_choke
export compressor_residuals
export compressor_residuals_scaled
export compressor_residual_scales
export solve_compressor_operating_point

end # module Compressor
