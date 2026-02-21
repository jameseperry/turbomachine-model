"""
Core EOS interface and registry helpers.
"""

"""
Base type for all equation-of-state models.
"""
abstract type AbstractEOS end

_not_implemented(eos::AbstractEOS, fn::Symbol) =
    error("$fn not implemented for EOS type $(nameof(typeof(eos)))")

function enthalpy_from_temperature(
    eos::AbstractEOS,
    temperature::Real,
)
    _not_implemented(eos, :enthalpy_from_temperature)
end

function temperature_from_enthalpy(
    eos::AbstractEOS,
    enthalpy::Real,
)
    _not_implemented(eos, :temperature_from_enthalpy)
end

function density_from_pressure_temperature(
    eos::AbstractEOS,
    pressure::Real,
    temperature::Real,
)
    _not_implemented(eos, :density_from_pressure_temperature)
end

function density_from_pressure_enthalpy(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :density_from_pressure_enthalpy)
end

function speed_of_sound_from_temperature(
    eos::AbstractEOS,
    temperature::Real,
)
    _not_implemented(eos, :speed_of_sound_from_temperature)
end

function speed_of_sound_from_enthalpy(
    eos::AbstractEOS,
    enthalpy::Real,
)
    _not_implemented(eos, :speed_of_sound_from_enthalpy)
end

function temperature(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :temperature)
end

function entropy(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :entropy)
end

function density(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :density)
end

function speed_of_sound(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :speed_of_sound)
end

function heat_capacity_cp(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :heat_capacity_cp)
end

function dynamic_viscosity(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :dynamic_viscosity)
end

function thermal_conductivity(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :thermal_conductivity)
end

function phase(
    eos::AbstractEOS,
    pressure::Real,
    enthalpy::Real,
)
    _not_implemented(eos, :phase)
end

function enthalpy_from_pressure_entropy(
    eos::AbstractEOS,
    pressure::Real,
    entropy::Real,
)
    _not_implemented(eos, :enthalpy_from_pressure_entropy)
end

"""
Default isentropic enthalpy helper for any EOS model.

Computes:
1. `entropy_1 = entropy(eos, pressure_1, enthalpy_1)`
2. `h_2s = enthalpy_from_pressure_entropy(eos, pressure_2, entropy_1)`
"""
function isentropic_enthalpy(
    eos::AbstractEOS,
    pressure_1::Real,
    enthalpy_1::Real,
    pressure_2::Real,
)
    entropy_1 = entropy(eos, pressure_1, enthalpy_1)
    return enthalpy_from_pressure_entropy(eos, pressure_2, entropy_1)
end

"""
EOS registry mapping composition ids to concrete EOS models.
"""
const EquationsOfStateRegistry = Dict{Symbol,AbstractEOS}

"""
Build a registry using placeholder real-fluid EOS models.

`AirEOS` and `SteamEOS` currently raise `not implemented` errors until their
property functions are filled in.
"""
function real_EOS()
    return EquationsOfStateRegistry(
        :air => AirEOS(),
        :steam => SteamEOS(),
    )
end

"""
Build a registry using ideal-gas EOS models for air and steam.
"""
function ideal_EOS()
    air = IdealGasEOS(:air; gas_constant=287.05, gamma=1.4)
    steam = IdealGasEOS(:steam; gas_constant=461.5, gamma=1.33)
    return EquationsOfStateRegistry(
        :air => air,
        :steam => steam,
    )
end
