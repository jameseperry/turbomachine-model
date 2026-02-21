"""
Air EOS marker.

Placeholder type for future air-property implementation.
"""
struct AirEOS <: AbstractEOS
end

function temperature_from_enthalpy(::AirEOS, enthalpy::Real)
    error("temperature_from_enthalpy for AirEOS is not implemented")
end

function enthalpy_from_temperature(::AirEOS, temperature::Real)
    error("enthalpy_from_temperature for AirEOS is not implemented")
end

function density_from_pressure_temperature(
    ::AirEOS,
    pressure::Real,
    temperature::Real,
)
    error("density_from_pressure_temperature for AirEOS is not implemented")
end

function density_from_pressure_enthalpy(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("density_from_pressure_enthalpy for AirEOS is not implemented")
end

function speed_of_sound_from_temperature(
    ::AirEOS,
    temperature::Real,
)
    error("speed_of_sound_from_temperature for AirEOS is not implemented")
end

function speed_of_sound_from_enthalpy(
    ::AirEOS,
    enthalpy::Real,
)
    error("speed_of_sound_from_enthalpy for AirEOS is not implemented")
end

function temperature(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("temperature for AirEOS is not implemented")
end

function entropy(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("entropy for AirEOS is not implemented")
end

function speed_of_sound(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("speed_of_sound for AirEOS is not implemented")
end

function heat_capacity_cp(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("heat_capacity_cp for AirEOS is not implemented")
end

function dynamic_viscosity(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("dynamic_viscosity for AirEOS is not implemented")
end

function thermal_conductivity(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("thermal_conductivity for AirEOS is not implemented")
end

function phase(
    ::AirEOS,
    pressure::Real,
    enthalpy::Real,
)
    error("phase for AirEOS is not implemented")
end

function enthalpy_from_pressure_entropy(
    ::AirEOS,
    pressure::Real,
    entropy::Real,
)
    error("enthalpy_from_pressure_entropy for AirEOS is not implemented")
end
