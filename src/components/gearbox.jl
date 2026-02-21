"""
Idealized two-port gearbox component.

Fields:
- `ratio`: speed ratio `omega_out / omega_in`.
- `efficiency`: mechanical efficiency in `(0, 1]`.
- `init`: named-tuple of component-specific initial condition values.
- `port_list`: component ports.
- `variable_list`: component variables.
"""
struct Gearbox <: AbstractComponent
    ratio::Float64
    efficiency::Float64
    init::NamedTuple
    port_list::Vector{ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function Gearbox(;
    ratio::Real,
    efficiency::Real=1.0,
    init::NamedTuple=NamedTuple(),
)
    ratio_f = Float64(ratio)
    efficiency_f = Float64(efficiency)
    ratio_f > 0.0 || error("Gearbox.ratio must be > 0")
    0.0 < efficiency_f <= 1.0 || error("Gearbox.efficiency must be in (0, 1]")
    port_list = ComponentPort[
        ComponentPort(shape_id=:ShaftPort, name=:input),
        ComponentPort(shape_id=:ShaftPort, name=:output),
    ]

    variable_list = ComponentVariable[
        ComponentVariable(
            id=:input_omega,
            unit=:rad_per_s,
            mappings=[(port=:input, label=:omega)],
        ),
        ComponentVariable(
            id=:input_tau,
            unit=:N_m,
            mappings=[(port=:input, label=:tau)],
        ),
        ComponentVariable(
            id=:output_omega,
            unit=:rad_per_s,
            mappings=[(port=:output, label=:omega)],
        ),
        ComponentVariable(
            id=:output_tau,
            unit=:N_m,
            mappings=[(port=:output, label=:tau)],
        ),
    ]

    return Gearbox(ratio_f, efficiency_f, init, port_list, variable_list)
end

ports(c::Gearbox) = c.port_list

function variables(c::Gearbox, sim_type::Symbol)
    if sim_type === :steady || sim_type === :transient
        return c.variable_list
    end
    error("unsupported sim_type: $sim_type (expected :steady or :transient)")
end
