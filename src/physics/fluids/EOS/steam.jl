"""
Steam composition marker.

Placeholder type for future steam-property implementation.
"""
struct SteamComposition <: AbstractComposition
end

function temperature_from_enthalpy(::SteamComposition, enthalpy::Real)
    error("temperature_from_enthalpy for SteamComposition is not implemented")
end

function enthalpy_from_temperature(::SteamComposition, temperature::Real)
    error("enthalpy_from_temperature for SteamComposition is not implemented")
end

function density_from_pressure_temperature(
    ::SteamComposition,
    pressure::Real,
    temperature::Real,
)
    error("density_from_pressure_temperature for SteamComposition is not implemented")
end

function density_from_pressure_enthalpy(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("density_from_pressure_enthalpy for SteamComposition is not implemented")
end

function speed_of_sound_from_temperature(
    ::SteamComposition,
    temperature::Real,
)
    error("speed_of_sound_from_temperature for SteamComposition is not implemented")
end

function speed_of_sound_from_enthalpy(
    ::SteamComposition,
    enthalpy::Real,
)
    error("speed_of_sound_from_enthalpy for SteamComposition is not implemented")
end

function temperature(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("temperature for SteamComposition is not implemented")
end

function entropy(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("entropy for SteamComposition is not implemented")
end

function speed_of_sound(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("speed_of_sound for SteamComposition is not implemented")
end

function heat_capacity_cp(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("heat_capacity_cp for SteamComposition is not implemented")
end

function dynamic_viscosity(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("dynamic_viscosity for SteamComposition is not implemented")
end

function thermal_conductivity(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("thermal_conductivity for SteamComposition is not implemented")
end

function phase(
    ::SteamComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("phase for SteamComposition is not implemented")
end

function enthalpy_from_pressure_entropy(
    ::SteamComposition,
    pressure::Real,
    entropy::Real,
)
    error("enthalpy_from_pressure_entropy for SteamComposition is not implemented")
end
