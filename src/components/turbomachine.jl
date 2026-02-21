"""
Unified compressor/turbine component.

Fields:
- `mode`: operation mode, `:compressor` or `:turbine`.
- `performance_map`: referenced compressor/turbine performance map data.
- `eta_guess`: nominal isentropic efficiency guess used for initialization.
- `init`: named-tuple of component-specific initial condition values.
- `port_list`: component ports.
- `steady_variable_list`: component variables for steady simulations.
- `transient_variable_list`: component variables for transient simulations.
"""
struct Turbomachine <: AbstractComponent
    mode::Symbol
    performance_map::PerformanceMap
    eta_guess::Float64
    init::NamedTuple
    port_list::Vector{ComponentPort}
    steady_variable_list::Vector{ComponentVariable}
    transient_variable_list::Vector{ComponentVariable}
end

function Turbomachine(;
    mode::Symbol,
    performance_map::PerformanceMap,
    eta_guess::Real,
    init::NamedTuple=NamedTuple(),
)
    mode in (:compressor, :turbine) ||
        error("Turbomachine.mode must be :compressor or :turbine")
    eta_f = Float64(eta_guess)
    0.0 < eta_f <= 1.0 ||
        error("Turbomachine.eta_guess must be in (0, 1]")

    port_list = ComponentPort[
        ComponentPort(shape_id=:FluidPort, name=:inlet),
        ComponentPort(shape_id=:FluidPort, name=:outlet),
        ComponentPort(shape_id=:ShaftPort, name=:shaft),
    ]

    steady_variable_list = ComponentVariable[
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
        ComponentVariable(
            id=:shaft_omega,
            unit=:rad_per_s,
            mappings=[(port=:shaft, label=:omega)],
        ),
        ComponentVariable(
            id=:shaft_tau,
            unit=:N_m,
            mappings=[(port=:shaft, label=:tau)],
        ),
    ]

    transient_variable_list = ComponentVariable[
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
            id=:inlet_mdot,
            unit=:kg_per_s,
            mappings=[(port=:inlet, label=:mdot)],
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
        ComponentVariable(
            id=:outlet_mdot,
            unit=:kg_per_s,
            mappings=[(port=:outlet, label=:mdot)],
        ),
        ComponentVariable(
            id=:shaft_omega,
            unit=:rad_per_s,
            mappings=[(port=:shaft, label=:omega)],
        ),
        ComponentVariable(
            id=:shaft_tau,
            unit=:N_m,
            mappings=[(port=:shaft, label=:tau)],
        ),
    ]

    return Turbomachine(
        mode,
        performance_map,
        eta_f,
        init,
        port_list,
        steady_variable_list,
        transient_variable_list,
    )
end

ports(c::Turbomachine) = c.port_list

function variables(c::Turbomachine, sim_type::Symbol)
    if sim_type === :steady
        return c.steady_variable_list
    elseif sim_type === :transient
        return c.transient_variable_list
    end
    error("unsupported sim_type: $sim_type (expected :steady or :transient)")
end
