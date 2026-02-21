"""
Lumped combustor component.

Fields:
- `dp_frac`: fractional total-pressure drop across the combustor.
- `fuel_lhv`: fuel lower heating value (J/kg).
- `init`: named-tuple of component-specific initial condition values.
- `port_list`: component ports.
- `variable_list`: component variables.
"""
struct Combustor <: AbstractComponent
    dp_frac::Float64
    fuel_lhv::Float64
    init::NamedTuple
    port_list::Vector{ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function Combustor(
    dp_frac::Real,
    fuel_lhv::Real,
    init::NamedTuple=NamedTuple(),
)
    port_list = ComponentPort[
        ComponentPort(shape_id=:FluidPort, name=:inlet),
        ComponentPort(shape_id=:FluidPort, name=:outlet),
    ]

    variable_list = ComponentVariable[
        ComponentVariable(
            id=:inlet_pt,
            unit=:Pa,
            mappings=[(port=:inlet, label=:pt)],
        ),
        ComponentVariable(
            id=:inlet_ht,
            unit=:J_per_kg,
            mappings=[(port=:inlet, label=:ht)],
        ),
        ComponentVariable(
            id=:composition,
            unit=:species,
            mappings=[(port=:inlet, label=:composition), (port=:outlet, label=:composition)],
        ),
        ComponentVariable(
            id=:mdot,
            unit=:kg_per_s,
            mappings=[(port=:inlet, label=:mdot), (port=:outlet, label=:mdot)],
        ),
        ComponentVariable(
            id=:outlet_pt,
            unit=:Pa,
            mappings=[(port=:outlet, label=:pt)],
        ),
        ComponentVariable(
            id=:outlet_ht,
            unit=:J_per_kg,
            mappings=[(port=:outlet, label=:ht)],
        ),
    ]

    return Combustor(
        Float64(dp_frac),
        Float64(fuel_lhv),
        init,
        port_list,
        variable_list,
    )
end

ports(c::Combustor) = c.port_list

function variables(c::Combustor, sim_type::Symbol)
    if sim_type === :steady || sim_type === :transient
        return c.variable_list
    end
    error("unsupported sim_type: $sim_type (expected :steady or :transient)")
end
