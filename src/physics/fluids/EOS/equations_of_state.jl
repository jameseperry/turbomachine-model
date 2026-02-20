"""
Bundle of equation-of-state functions for a specific composition.

Each member is a callable function:
- `temperature(p, h)`
- `entropy(p, h)`
- `density(p, h)`
- `speed_of_sound(p, h)`
- `heat_capacity_cp(p, h)`
- `dynamic_viscosity(p, h)`
- `thermal_conductivity(p, h)`
- `phase(p, h)`
- `enthalpy_from_pressure_entropy(p, s)`
"""
struct EquationsOfState{
    FTemperature,
    FEntropy,
    FDensity,
    FSpeedOfSound,
    FHeatCapacityCp,
    FDynamicViscosity,
    FThermalConductivity,
    FPhase,
    FEnthalpyFromPressureEntropy,
}
    temperature::FTemperature
    entropy::FEntropy
    density::FDensity
    speed_of_sound::FSpeedOfSound
    heat_capacity_cp::FHeatCapacityCp
    dynamic_viscosity::FDynamicViscosity
    thermal_conductivity::FThermalConductivity
    phase::FPhase
    enthalpy_from_pressure_entropy::FEnthalpyFromPressureEntropy
end

"""
Registry mapping composition ids to equation-of-state bundles.
"""
const EquationsOfStateRegistry = Dict{Symbol,EquationsOfState}

"""
Create an `EquationsOfState` bundle from a composition object.
"""
function EquationsOfState(composition::AbstractComposition)
    return EquationsOfState(
        (pressure, enthalpy) -> temperature(composition, pressure, enthalpy),
        (pressure, enthalpy) -> entropy(composition, pressure, enthalpy),
        (pressure, enthalpy) -> density(composition, pressure, enthalpy),
        (pressure, enthalpy) -> speed_of_sound(composition, pressure, enthalpy),
        (pressure, enthalpy) -> heat_capacity_cp(composition, pressure, enthalpy),
        (pressure, enthalpy) -> dynamic_viscosity(composition, pressure, enthalpy),
        (pressure, enthalpy) -> thermal_conductivity(composition, pressure, enthalpy),
        (pressure, enthalpy) -> phase(composition, pressure, enthalpy),
        (pressure, entropy) -> enthalpy_from_pressure_entropy(composition, pressure, entropy),
    )
end

"""
Build a registry using the placeholder real-fluid composition types.

`AirComposition` and `SteamComposition` methods currently raise
`not implemented` errors until their property functions are filled in.
"""
function real_EOS()
    return EquationsOfStateRegistry(
        :air => EquationsOfState(AirComposition()),
        :steam => EquationsOfState(SteamComposition()),
    )
end

"""
Build a registry using ideal-gas compositions for air and steam.
"""
function ideal_EOS()
    air = IdealGasComposition(:air; gas_constant=287.05, gamma=1.4)
    steam = IdealGasComposition(:steam; gas_constant=461.5, gamma=1.33)
    return EquationsOfStateRegistry(
        :air => EquationsOfState(air),
        :steam => EquationsOfState(steam),
    )
end

function _get_equation_of_state(registry::EquationsOfStateRegistry, composition::Symbol)
    haskey(registry, composition) ||
        error("unknown fluid composition symbol: $composition")
    return registry[composition]
end

function temperature(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.temperature(pressure, enthalpy)
end

function entropy(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.entropy(pressure, enthalpy)
end

function density(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.density(pressure, enthalpy)
end

function speed_of_sound(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.speed_of_sound(pressure, enthalpy)
end

function heat_capacity_cp(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.heat_capacity_cp(pressure, enthalpy)
end

function dynamic_viscosity(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.dynamic_viscosity(pressure, enthalpy)
end

function thermal_conductivity(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.thermal_conductivity(pressure, enthalpy)
end

function phase(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.phase(pressure, enthalpy)
end

function enthalpy_from_pressure_entropy(
    registry::EquationsOfStateRegistry,
    composition::Symbol,
    pressure::Real,
    entropy::Real,
)
    equation_of_state = _get_equation_of_state(registry, composition)
    return equation_of_state.enthalpy_from_pressure_entropy(pressure, entropy)
end
