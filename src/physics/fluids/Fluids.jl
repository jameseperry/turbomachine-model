module Fluids

include("stagnation_relations.jl")
include("flow_energy_relations.jl")
include("turbomachine_performance_map.jl")

export static_temperature_from_total, total_temperature_from_static
export static_pressure_from_total, total_pressure_from_static
export static_enthalpy_from_total, total_enthalpy_from_static
export PerformanceMap, corrected_speed, corrected_flow
export map_pr_eta, map_pr_eta_from_stagnation
export demo_compressor_map, demo_turbine_map

end # module Fluids
