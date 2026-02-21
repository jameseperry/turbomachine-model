"""
Lumped combustor component.

Fields:
- `dp_frac`: fractional total-pressure drop across the combustor.
- `fuel_lhv`: fuel lower heating value (J/kg).
- `init`: named-tuple of component-specific initial condition values.
- `port_map`: component ports by name.
- `variable_list`: component variables.
"""
struct Combustor <: AbstractComponent
    dp_frac::Float64
    fuel_lhv::Float64
    init::NamedTuple
    port_map::Dict{Symbol,ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function Combustor(
    dp_frac::Real,
    fuel_lhv::Real,
    init::NamedTuple=NamedTuple(),
)
    port_map = Dict(
        :inlet => ComponentPort(
            :FluidPort;
            variable_ids=Dict(
                :pt => :inlet_pt,
                :ht => :inlet_ht,
                :composition => :inlet_composition,
                :mdot => :inlet_mdot,
            ),
        ),
        :outlet => ComponentPort(
            :FluidPort;
            variable_ids=Dict(
                :pt => :outlet_pt,
                :ht => :outlet_ht,
                :composition => :outlet_composition,
                :mdot => :outlet_mdot,
            ),
        ),
    )

    variable_list = ComponentVariable[
        ComponentVariable(:inlet_pt, :Pa),
        ComponentVariable(:inlet_ht, :J_per_kg),
        ComponentVariable(:inlet_composition, :species),
        ComponentVariable(:inlet_mdot, :kg_per_s),
        ComponentVariable(:outlet_pt, :Pa),
        ComponentVariable(:outlet_ht, :J_per_kg),
        ComponentVariable(:outlet_composition, :species),
        ComponentVariable(:outlet_mdot, :kg_per_s),
    ]

    return Combustor(Float64(dp_frac), Float64(fuel_lhv), init, port_map, variable_list)
end

ports(c::Combustor) = c.port_map

variables(c::Combustor) = c.variable_list
