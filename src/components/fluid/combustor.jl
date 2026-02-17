"""
Lumped combustor component.

Fields:
- `dp_frac`: fractional total-pressure drop across the combustor.
- `fuel_lhv`: fuel lower heating value (J/kg).
- `init`: named-tuple of component-specific initial condition values.
"""
struct Combustor <: AbstractComponent
    dp_frac::Float64
    fuel_lhv::Float64
    init::NamedTuple
end

function port_specs(::Combustor)
    Dict(
        :inlet => FLUID_PORT,
        :outlet => FLUID_PORT,
    )
end

required_ports(::Combustor) = [:inlet, :outlet]
