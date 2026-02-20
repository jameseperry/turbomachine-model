"""
Air composition marker.

Placeholder type for future air-property implementation.
"""
struct AirComposition <: AbstractComposition
end

function temperature_from_enthalpy(::AirComposition, enthalpy::Real)
    error("temperature_from_enthalpy for AirComposition is not implemented")
end

function enthalpy_from_temperature(::AirComposition, temperature::Real)
    error("enthalpy_from_temperature for AirComposition is not implemented")
end

function density_from_pressure_temperature(
    ::AirComposition,
    pressure::Real,
    temperature::Real,
)
    error("density_from_pressure_temperature for AirComposition is not implemented")
end

function density_from_pressure_enthalpy(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("density_from_pressure_enthalpy for AirComposition is not implemented")
end

function speed_of_sound_from_temperature(
    ::AirComposition,
    temperature::Real,
)
    error("speed_of_sound_from_temperature for AirComposition is not implemented")
end

function speed_of_sound_from_enthalpy(
    ::AirComposition,
    enthalpy::Real,
)
    error("speed_of_sound_from_enthalpy for AirComposition is not implemented")
end

function temperature(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("temperature for AirComposition is not implemented")
end

function entropy(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("entropy for AirComposition is not implemented")
end

function speed_of_sound(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("speed_of_sound for AirComposition is not implemented")
end

function heat_capacity_cp(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("heat_capacity_cp for AirComposition is not implemented")
end

function dynamic_viscosity(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("dynamic_viscosity for AirComposition is not implemented")
end

function thermal_conductivity(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("thermal_conductivity for AirComposition is not implemented")
end

function phase(
    ::AirComposition,
    pressure::Real,
    enthalpy::Real,
)
    error("phase for AirComposition is not implemented")
end

function enthalpy_from_pressure_entropy(
    ::AirComposition,
    pressure::Real,
    entropy::Real,
)
    error("enthalpy_from_pressure_entropy for AirComposition is not implemented")
end
