module Fluids

include("EOS/equations_of_state.jl")
include("EOS/ideal_gas.jl")
include("EOS/air.jl")
include("EOS/steam.jl")
include("stagnation_relations.jl")
include("flow_energy_relations.jl")

export AbstractEOS, AirEOS, SteamEOS, IdealGasEOS
export EquationsOfStateRegistry
export real_EOS, ideal_EOS
export enthalpy_from_temperature, temperature_from_enthalpy
export density_from_pressure_temperature, density_from_pressure_enthalpy
export speed_of_sound_from_temperature, speed_of_sound_from_enthalpy
export temperature, entropy, density, speed_of_sound
export heat_capacity_cp, dynamic_viscosity, thermal_conductivity
export phase, enthalpy_from_pressure_entropy, isentropic_enthalpy
export static_temperature_from_total, total_temperature_from_static
export static_pressure_from_total, total_pressure_from_static
export static_enthalpy_from_total, total_enthalpy_from_static
export velocity_from_massflow, velocity_from_ph_mdot

end # module Fluids
