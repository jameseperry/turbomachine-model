module Turbine

using ...Fluids

include("map_tabulated.jl")
include("residuals.jl")

export AbstractTurbinePerformanceMap
export TabulatedTurbinePerformanceMap
export corrected_speed, corrected_flow
export turbine_performance_map, turbine_performance_map_from_stagnation
export demo_tabulated_turbine_performance_map
export write_toml, read_toml
export turbine_residuals
export turbine_residuals_scaled
export turbine_residual_scales
export solve_turbine_operating_point

end # module Turbine
