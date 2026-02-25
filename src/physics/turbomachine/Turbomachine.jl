module Turbomachine

using ..Fluids

include("tabulated_performance_map.jl")
include("analytic_performance_map.jl")
include("residuals.jl")

export AbstractPerformanceMap, TabulatedPerformanceMap, AnalyticPerformanceMap
export corrected_speed, corrected_flow
export performance_map, performance_map_from_stagnation
export demo_compressor_performance_map, demo_turbine_performance_map
export turbomachine_residuals
export turbomachine_residuals_scaled
export turbomachine_residual_scales
export solve_turbomachine_operating_point

end # module Turbomachine
