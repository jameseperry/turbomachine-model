"""
Lumped plenum volume component with one inlet and one outlet fluid port.

Fields:
- `volume`: plenum control-volume size in m^3.
- `init`: named-tuple of component-specific initial condition values.
- `port_map`: component ports by name.
- `variable_list`: component variables.
"""
struct Plenum <: AbstractComponent
    volume::Float64
    init::NamedTuple
    port_map::Dict{Symbol,ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function Plenum(;
    volume::Real,
    init::NamedTuple=NamedTuple(),
)
    volume_f = Float64(volume)
    volume_f > 0.0 || error("Plenum.volume must be > 0")
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

    return Plenum(volume_f, init, port_map, variable_list)
end

ports(c::Plenum) = c.port_map

variables(c::Plenum) = c.variable_list
