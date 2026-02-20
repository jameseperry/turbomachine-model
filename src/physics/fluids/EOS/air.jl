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
