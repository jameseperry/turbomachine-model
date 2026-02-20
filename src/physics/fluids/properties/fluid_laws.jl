"""
Bundle of fluid-property law functions for a specific composition.

Each member is a callable function:
- `enthalpy_from_temperature(T)`
- `temperature_from_enthalpy(h)`
- `density_from_pressure_temperature(p, T)`
- `density_from_pressure_enthalpy(p, h)`
- `speed_of_sound_from_temperature(T)`
- `speed_of_sound_from_enthalpy(h)`
"""
struct FluidLaws{
    FEnthalpyFromTemperature,
    FTemperatureFromEnthalpy,
    FDensityFromPressureTemperature,
    FDensityFromPressureEnthalpy,
    FSpeedOfSoundFromTemperature,
    FSpeedOfSoundFromEnthalpy,
}
    enthalpy_from_temperature::FEnthalpyFromTemperature
    temperature_from_enthalpy::FTemperatureFromEnthalpy
    density_from_pressure_temperature::FDensityFromPressureTemperature
    density_from_pressure_enthalpy::FDensityFromPressureEnthalpy
    speed_of_sound_from_temperature::FSpeedOfSoundFromTemperature
    speed_of_sound_from_enthalpy::FSpeedOfSoundFromEnthalpy
end

"""
Registry mapping composition ids to fluid-law bundles.
"""
const FluidLawsRegistry = Dict{Symbol,FluidLaws}

"""
Create a `FluidLaws` bundle from a composition object.
"""
function FluidLaws(composition::AbstractComposition)
    return FluidLaws(
        temperature -> enthalpy_from_temperature(composition, temperature),
        enthalpy -> temperature_from_enthalpy(composition, enthalpy),
        (pressure, temperature) -> density_from_pressure_temperature(composition, pressure, temperature),
        (pressure, enthalpy) -> density_from_pressure_enthalpy(composition, pressure, enthalpy),
        temperature -> speed_of_sound_from_temperature(composition, temperature),
        enthalpy -> speed_of_sound_from_enthalpy(composition, enthalpy),
    )
end

"""
Build a registry using the placeholder real-fluid composition types.

`AirComposition` and `SteamComposition` methods currently raise
`not implemented` errors until their property functions are filled in.
"""
function real_fluid_laws()
    return FluidLawsRegistry(
        :air => FluidLaws(AirComposition()),
        :steam => FluidLaws(SteamComposition()),
    )
end

"""
Build a registry using ideal-gas compositions for air and steam.
"""
function ideal_fluid_laws()
    air = IdealGasComposition(:air; gas_constant=287.05, gamma=1.4)
    steam = IdealGasComposition(:steam; gas_constant=461.5, gamma=1.33)
    return FluidLawsRegistry(
        :air => FluidLaws(air),
        :steam => FluidLaws(steam),
    )
end

function _get_laws(registry::FluidLawsRegistry, composition::Symbol)
    haskey(registry, composition) ||
        error("unknown fluid composition symbol: $composition")
    return registry[composition]
end

function enthalpy_from_temperature(
    registry::FluidLawsRegistry,
    composition::Symbol,
    temperature::Real,
)
    laws = _get_laws(registry, composition)
    return laws.enthalpy_from_temperature(temperature)
end

function temperature_from_enthalpy(
    registry::FluidLawsRegistry,
    composition::Symbol,
    enthalpy::Real,
)
    laws = _get_laws(registry, composition)
    return laws.temperature_from_enthalpy(enthalpy)
end

function density_from_pressure_temperature(
    registry::FluidLawsRegistry,
    composition::Symbol,
    pressure::Real,
    temperature::Real,
)
    laws = _get_laws(registry, composition)
    return laws.density_from_pressure_temperature(pressure, temperature)
end

function density_from_pressure_enthalpy(
    registry::FluidLawsRegistry,
    composition::Symbol,
    pressure::Real,
    enthalpy::Real,
)
    laws = _get_laws(registry, composition)
    return laws.density_from_pressure_enthalpy(pressure, enthalpy)
end

function speed_of_sound_from_temperature(
    registry::FluidLawsRegistry,
    composition::Symbol,
    temperature::Real,
)
    laws = _get_laws(registry, composition)
    return laws.speed_of_sound_from_temperature(temperature)
end

function speed_of_sound_from_enthalpy(
    registry::FluidLawsRegistry,
    composition::Symbol,
    enthalpy::Real,
)
    laws = _get_laws(registry, composition)
    return laws.speed_of_sound_from_enthalpy(enthalpy)
end
