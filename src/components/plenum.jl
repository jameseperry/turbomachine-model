"""
Lumped plenum volume component with one inlet and one outlet fluid port.

Fields:
- `volume`: plenum control-volume size in m^3.
- `init`: named-tuple of component-specific initial condition values.
- `port_list`: component ports.
- `steady_variable_list`: component variables for steady simulations.
- `transient_variable_list`: component variables for transient simulations.
"""
struct Plenum <: AbstractComponent
    volume::Float64
    init::NamedTuple
    port_list::Vector{ComponentPort}
    steady_variable_list::Vector{ComponentVariable}
    transient_variable_list::Vector{ComponentVariable}
end

function Plenum(;
    volume::Real,
    init::NamedTuple=NamedTuple(),
)
    volume_f = Float64(volume)
    volume_f > 0.0 || error("Plenum.volume must be > 0")
    port_list = ComponentPort[
        ComponentPort(shape_id=:FluidPort, name=:inlet),
        ComponentPort(shape_id=:FluidPort, name=:outlet),
    ]

    steady_variable_list = ComponentVariable[
        ComponentVariable(
            id=:pt,
            unit=:Pa,
            mappings=[(port=:inlet, label=:pt), (port=:outlet, label=:pt)],
        ),
        ComponentVariable(
            id=:ht,
            unit=:J_per_kg,
            mappings=[(port=:inlet, label=:ht), (port=:outlet, label=:ht)],
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
    ]

    transient_variable_list = ComponentVariable[
        ComponentVariable(
            id=:pt,
            unit=:Pa,
            mappings=[(port=:inlet, label=:pt), (port=:outlet, label=:pt)],
        ),
        ComponentVariable(
            id=:ht,
            unit=:J_per_kg,
            mappings=[(port=:inlet, label=:ht), (port=:outlet, label=:ht)],
        ),
        ComponentVariable(
            id=:composition,
            unit=:species,
            mappings=[(port=:inlet, label=:composition), (port=:outlet, label=:composition)],
        ),
        ComponentVariable(
            id=:inlet_mdot,
            unit=:kg_per_s,
            mappings=[(port=:inlet, label=:mdot)],
        ),
        ComponentVariable(
            id=:outlet_mdot,
            unit=:kg_per_s,
            mappings=[(port=:outlet, label=:mdot)],
        ),
        ComponentVariable(
            id=:retained_mass,
            unit=:kg,
        ),
    ]

    return Plenum(volume_f, init, port_list, steady_variable_list, transient_variable_list)
end

ports(c::Plenum) = c.port_list

function variables(c::Plenum, sim_type::Symbol)
    if sim_type === :steady
        return c.steady_variable_list
    elseif sim_type === :transient
        return c.transient_variable_list
    end
    error("unsupported sim_type: $sim_type (expected :steady or :transient)")
end
