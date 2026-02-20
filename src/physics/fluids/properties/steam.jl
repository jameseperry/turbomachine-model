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
