"""
Rigid inertial shaft/spool component with configurable attachment count.

Fields:
- `J`: lumped rotational inertia (kg*m^2).
- `damping`: linear viscous damping coefficient about the shaft axis.
- `n_ports`: number of shaft attachment ports to expose.
- `init`: named-tuple of component-specific initial condition values (for example `omega`).
- `port_list`: component ports.
- `variable_list`: component variables.
"""
struct InertialShaft <: AbstractComponent
    J::Float64
    damping::Float64
    n_ports::Int
    init::NamedTuple
    port_list::Vector{ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function InertialShaft(;
    J::Real,
    damping::Real=0.0,
    n_ports::Integer=2,
    init::NamedTuple=NamedTuple(),
)
    J_f = Float64(J)
    damping_f = Float64(damping)
    n_ports_i = Int(n_ports)
    J_f >= 0.0 || error("InertialShaft.J must be >= 0")
    damping_f >= 0.0 || error("InertialShaft.damping must be >= 0")
    n_ports_i >= 1 || error("InertialShaft.n_ports must be >= 1")
    port_list = ComponentPort[]
    omega_mappings = NamedTuple{(:port,:label),Tuple{Symbol,Symbol}}[]
    variable_list = ComponentVariable[]

    for i in 1:n_ports_i
        name = Symbol("shaft", i)
        tau_id = Symbol(string(name), "_tau")
        push!(port_list, ComponentPort(shape_id=:ShaftPort, name=name))
        push!(omega_mappings, (port=name, label=:omega))
        push!(
            variable_list,
            ComponentVariable(
                id=tau_id,
                unit=:N_m,
                mappings=[(port=name, label=:tau)],
            ),
        )
    end

    push!(
        variable_list,
        ComponentVariable(
            id=:omega,
            unit=:rad_per_s,
            mappings=omega_mappings,
        ),
    )

    return InertialShaft(J_f, damping_f, n_ports_i, init, port_list, variable_list)
end

ports(c::InertialShaft) = c.port_list

function variables(c::InertialShaft, sim_type::Symbol)
    if sim_type === :steady || sim_type === :transient
        return c.variable_list
    end
    error("unsupported sim_type: $sim_type (expected :steady or :transient)")
end
