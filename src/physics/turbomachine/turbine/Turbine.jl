module Turbine

using ...Fluids

include("api.jl")
include("map_tabulated.jl")
include("meanline_model.jl")
include("residuals.jl")

export AbstractTurbinePerformanceMap
export TabulatedTurbinePerformanceMap
export TurbineMeanlineModel
export corrected_speed, corrected_flow
export turbine_performance_map, turbine_performance_map_from_stagnation
export performance_map_domain
export demo_tabulated_turbine_performance_map
export tabulate_turbine_meanline_model
export demo_turbine_meanline_model
export write_toml, read_toml
export turbine_residuals
export turbine_residuals_scaled
export turbine_residual_scales
export solve_turbine_operating_point

end # module Turbine
