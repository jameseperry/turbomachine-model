module Fluids

include("properties/composition.jl")
include("properties/ideal_gas.jl")
include("properties/air.jl")
include("properties/steam.jl")
include("properties/fluid_laws.jl")
include("stagnation_relations.jl")
include("flow_energy_relations.jl")
include("turbomachine_performance_map.jl")

export AbstractComposition, AirComposition, SteamComposition, IdealGasComposition
export FluidLaws, FluidLawsRegistry
export real_fluid_laws, ideal_fluid_laws
export enthalpy_from_temperature, temperature_from_enthalpy
export density_from_pressure_temperature, density_from_pressure_enthalpy
export speed_of_sound_from_temperature, speed_of_sound_from_enthalpy
export static_temperature_from_total, total_temperature_from_static
export static_pressure_from_total, total_pressure_from_static
export static_enthalpy_from_total, total_enthalpy_from_static
export velocity_from_massflow, velocity_from_ph_mdot
export PerformanceMap, corrected_speed, corrected_flow
export map_pr_eta, map_pr_eta_from_stagnation
export demo_compressor_map, demo_turbine_map

end # module Fluids
