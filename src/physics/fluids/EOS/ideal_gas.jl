"""
Ideal-gas composition parameters.

Fields:
- `id`: composition identifier (e.g. `:air`, `:steam`).
- `gas_constant`: specific gas constant `R` in J/(kg*K).
- `specific_heat_cp`: constant-pressure specific heat `cp` in J/(kg*K).
- `specific_heat_cv`: constant-volume specific heat `cv` in J/(kg*K).
- `gamma`: heat-capacity ratio `cp/cv`.
"""
struct IdealGasComposition <: AbstractComposition
    id::Symbol
    gas_constant::Float64
    specific_heat_cp::Float64
    specific_heat_cv::Float64
    gamma::Float64
end

function IdealGasComposition(
    id::Symbol;
    gas_constant::Real,
    specific_heat_cp::Union{Nothing,Real}=nothing,
    specific_heat_cv::Union{Nothing,Real}=nothing,
    gamma::Union{Nothing,Real}=nothing,
)
    gas_constant > 0 || error("gas_constant must be > 0")

    if gamma !== nothing
        g = gamma
        g > 1 || error("gamma must be > 1")
        cv = gas_constant / (g - 1)
        cp = g * cv
        return IdealGasComposition(
            id,
            Float64(gas_constant),
            Float64(cp),
            Float64(cv),
            Float64(g),
        )
    end

    specific_heat_cp === nothing && error("specific_heat_cp is required when gamma is not provided")
    specific_heat_cv === nothing && error("specific_heat_cv is required when gamma is not provided")
    cp = specific_heat_cp
    cv = specific_heat_cv
    cp > 0 || error("specific_heat_cp must be > 0")
    cv > 0 || error("specific_heat_cv must be > 0")

    # Ideal-gas identity: cp - cv = R.
    isapprox(cp - cv, gas_constant; rtol=1e-6) ||
        error("ideal-gas identity violated: cp - cv must equal gas_constant")

    g = cp / cv
    g > 1 || error("gamma must be > 1")

    return IdealGasComposition(
        id,
        Float64(gas_constant),
        Float64(cp),
        Float64(cv),
        Float64(g),
    )
end

"""
Specific enthalpy from temperature: `h = cp * T`.
"""
function enthalpy_from_temperature(comp::IdealGasComposition, temperature::Real)
    temperature > 0 || error("temperature must be > 0")
    return comp.specific_heat_cp * temperature
end

"""
Temperature from specific enthalpy: `T = h / cp`.
"""
function temperature_from_enthalpy(comp::IdealGasComposition, enthalpy::Real)
    enthalpy > 0 || error("enthalpy must be > 0")
    return enthalpy / comp.specific_heat_cp
end

"""
Density from pressure and temperature: `rho = p / (R * T)`.
"""
function density_from_pressure_temperature(
    comp::IdealGasComposition,
    pressure::Real,
    temperature::Real,
)
    pressure > 0 || error("pressure must be > 0")
    temperature > 0 || error("temperature must be > 0")
    return pressure / (comp.gas_constant * temperature)
end

"""
Density from pressure and specific enthalpy.
"""
function density_from_pressure_enthalpy(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    temperature = temperature_from_enthalpy(comp, enthalpy)
    return density_from_pressure_temperature(comp, pressure, temperature)
end

"""
Speed of sound from temperature: `a = sqrt(gamma * R * T)`.
"""
function speed_of_sound_from_temperature(
    comp::IdealGasComposition,
    temperature::Real,
)
    temperature > 0 || error("temperature must be > 0")
    return sqrt(comp.gamma * comp.gas_constant * temperature)
end

"""
Speed of sound from specific enthalpy.
"""
function speed_of_sound_from_enthalpy(
    comp::IdealGasComposition,
    enthalpy::Real,
)
    temperature = temperature_from_enthalpy(comp, enthalpy)
    return speed_of_sound_from_temperature(comp, temperature)
end
