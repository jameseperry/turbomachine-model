"""
Steam EOS marker.

Placeholder type for future steam-property implementation.
"""
struct SteamEOS <: AbstractEOS
end

function temperature_from_enthalpy(::SteamEOS, enthalpy::Real)
    error("temperature_from_enthalpy for SteamEOS is not implemented")
end

function enthalpy_from_temperature(::SteamEOS, temperature::Real)
    error("enthalpy_from_temperature for SteamEOS is not implemented")
end

function density_from_pressure_temperature(
    ::SteamEOS,
    pressure::Real,
    temperature::Real,
)
    error("density_from_pressure_temperature for SteamEOS is not implemented")
end

function density_from_pressure_enthalpy(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("density_from_pressure_enthalpy for SteamEOS is not implemented")
end

function speed_of_sound_from_temperature(
    ::SteamEOS,
    temperature::Real,
)
    error("speed_of_sound_from_temperature for SteamEOS is not implemented")
end

function speed_of_sound_from_enthalpy(
    ::SteamEOS,
    enthalpy::Real,
)
    error("speed_of_sound_from_enthalpy for SteamEOS is not implemented")
end

function temperature(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("temperature for SteamEOS is not implemented")
end

function entropy(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("entropy for SteamEOS is not implemented")
end

function speed_of_sound(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("speed_of_sound for SteamEOS is not implemented")
end

function heat_capacity_cp(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("heat_capacity_cp for SteamEOS is not implemented")
end

function dynamic_viscosity(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("dynamic_viscosity for SteamEOS is not implemented")
end

function thermal_conductivity(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("thermal_conductivity for SteamEOS is not implemented")
end

function phase(
    ::SteamEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("phase for SteamEOS is not implemented")
end

function enthalpy_from_pressure_entropy(
    ::SteamEOS,
    pressure::Real,
    entropy::Real,
)
    error("enthalpy_from_pressure_entropy for SteamEOS is not implemented")
end
