"""
Ideal-gas composition parameters.

Fields:
- `id`: composition identifier (e.g. `:air`, `:steam`).
- `gas_constant`: specific gas constant `R` in J/(kg*K).
- `specific_heat_cp`: constant-pressure specific heat `cp` in J/(kg*K).
- `specific_heat_cv`: constant-volume specific heat `cv` in J/(kg*K).
- `gamma`: heat-capacity ratio `cp/cv`.
- `dynamic_viscosity_ref`: reference dynamic viscosity in Pa*s.
- `dynamic_viscosity_ref_temperature`: reference temperature for viscosity in K.
- `sutherland_constant`: Sutherland constant for viscosity model in K.
- `prandtl_number`: Prandtl number for conductivity estimate.
- `pressure_reference`: reference pressure for entropy in Pa.
- `temperature_reference`: reference temperature for entropy in K.
- `entropy_reference`: reference specific entropy in J/(kg*K).
"""
struct IdealGasComposition <: AbstractComposition
    id::Symbol
    gas_constant::Float64
    specific_heat_cp::Float64
    specific_heat_cv::Float64
    gamma::Float64
    dynamic_viscosity_ref::Float64
    dynamic_viscosity_ref_temperature::Float64
    sutherland_constant::Float64
    prandtl_number::Float64
    pressure_reference::Float64
    temperature_reference::Float64
    entropy_reference::Float64
end

function IdealGasComposition(
    id::Symbol;
    gas_constant::Real,
    specific_heat_cp::Union{Nothing,Real}=nothing,
    specific_heat_cv::Union{Nothing,Real}=nothing,
    gamma::Union{Nothing,Real}=nothing,
    dynamic_viscosity_ref::Real=1.716e-5,
    dynamic_viscosity_ref_temperature::Real=273.15,
    sutherland_constant::Real=110.4,
    prandtl_number::Real=0.71,
    pressure_reference::Real=101_325.0,
    temperature_reference::Real=300.0,
    entropy_reference::Real=0.0,
)
    gas_constant > 0 || error("gas_constant must be > 0")
    dynamic_viscosity_ref > 0 || error("dynamic_viscosity_ref must be > 0")
    dynamic_viscosity_ref_temperature > 0 || error("dynamic_viscosity_ref_temperature must be > 0")
    sutherland_constant > 0 || error("sutherland_constant must be > 0")
    prandtl_number > 0 || error("prandtl_number must be > 0")
    pressure_reference > 0 || error("pressure_reference must be > 0")
    temperature_reference > 0 || error("temperature_reference must be > 0")

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
            Float64(dynamic_viscosity_ref),
            Float64(dynamic_viscosity_ref_temperature),
            Float64(sutherland_constant),
            Float64(prandtl_number),
            Float64(pressure_reference),
            Float64(temperature_reference),
            Float64(entropy_reference),
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
        Float64(dynamic_viscosity_ref),
        Float64(dynamic_viscosity_ref_temperature),
        Float64(sutherland_constant),
        Float64(prandtl_number),
        Float64(pressure_reference),
        Float64(temperature_reference),
        Float64(entropy_reference),
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

function density(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    return density_from_pressure_enthalpy(comp, pressure, enthalpy)
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

function temperature(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    pressure > 0 || error("pressure must be > 0")
    return temperature_from_enthalpy(comp, enthalpy)
end

function entropy(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    pressure > 0 || error("pressure must be > 0")
    temperature = temperature_from_enthalpy(comp, enthalpy)
    return comp.entropy_reference +
           comp.specific_heat_cp * log(temperature / comp.temperature_reference) -
           comp.gas_constant * log(pressure / comp.pressure_reference)
end

function speed_of_sound(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    temp = temperature(comp, pressure, enthalpy)
    return speed_of_sound_from_temperature(comp, temp)
end

function heat_capacity_cp(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    pressure > 0 || error("pressure must be > 0")
    enthalpy > 0 || error("enthalpy must be > 0")
    return comp.specific_heat_cp
end

function dynamic_viscosity(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    temp = temperature(comp, pressure, enthalpy)
    t_ref = comp.dynamic_viscosity_ref_temperature
    s = comp.sutherland_constant
    return comp.dynamic_viscosity_ref *
           (temp / t_ref)^(3 / 2) *
           (t_ref + s) / (temp + s)
end

function thermal_conductivity(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    mu = dynamic_viscosity(comp, pressure, enthalpy)
    cp = heat_capacity_cp(comp, pressure, enthalpy)
    return mu * cp / comp.prandtl_number
end

function phase(
    comp::IdealGasComposition,
    pressure::Real,
    enthalpy::Real,
)
    pressure > 0 || error("pressure must be > 0")
    enthalpy > 0 || error("enthalpy must be > 0")
    return :gas
end

function enthalpy_from_pressure_entropy(
    comp::IdealGasComposition,
    pressure::Real,
    entropy::Real,
)
    pressure > 0 || error("pressure must be > 0")
    temperature = comp.temperature_reference * exp(
        (entropy - comp.entropy_reference +
         comp.gas_constant * log(pressure / comp.pressure_reference)) /
        comp.specific_heat_cp,
    )
    return enthalpy_from_temperature(comp, temperature)
end
