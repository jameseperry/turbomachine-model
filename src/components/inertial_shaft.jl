"""
Rigid inertial shaft/spool component with configurable attachment count.

Fields:
- `J`: lumped rotational inertia (kg*m^2).
- `damping`: linear viscous damping coefficient about the shaft axis.
- `n_ports`: number of shaft attachment ports to expose.
- `init`: named-tuple of component-specific initial condition values (for example `omega`).
- `port_map`: component ports by name.
- `variable_list`: component variables.
"""
struct InertialShaft <: AbstractComponent
    J::Float64
    damping::Float64
    n_ports::Int
    init::NamedTuple
    port_map::Dict{Symbol,ComponentPort}
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
    port_map = Dict{Symbol,ComponentPort}()
    variable_list = ComponentVariable[]

    for i in 1:n_ports_i
        name = Symbol("shaft", i)
        omega_id = Symbol(string(name), "_omega")
        tau_id = Symbol(string(name), "_tau")
        port_map[name] = ComponentPort(
            :ShaftPort;
            variable_ids=Dict(
                :omega => omega_id,
                :tau => tau_id,
            ),
        )
        push!(variable_list, ComponentVariable(omega_id, :rad_per_s))
        push!(variable_list, ComponentVariable(tau_id, :N_m))
    end

    return InertialShaft(J_f, damping_f, n_ports_i, init, port_map, variable_list)
end

ports(c::InertialShaft) = c.port_map

variables(c::InertialShaft) = c.variable_list
