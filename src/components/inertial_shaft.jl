"""
Rigid inertial shaft/spool component with configurable attachment count.

Fields:
- `J`: lumped rotational inertia (kg*m^2).
- `damping`: linear viscous damping coefficient about the shaft axis.
- `n_ports`: number of shaft attachment ports to expose.
- `init`: named-tuple of component-specific initial condition values (for example `omega`).
"""
struct InertialShaft <: AbstractComponent
    J::Float64
    damping::Float64
    n_ports::Int
    init::NamedTuple
end

function InertialShaft(;
    J::Real,
    damping::Real=0.0,
    n_ports::Integer=2,
    init::NamedTuple=NamedTuple(),
)
    comp = InertialShaft(Float64(J), Float64(damping), Int(n_ports), init)
    validate(comp)
    return comp
end

_shaft_port_names(n_ports::Int) = [Symbol("shaft", i) for i in 1:n_ports]

function port_specs(c::InertialShaft)
    Dict(name => SHAFT_PORT for name in _shaft_port_names(c.n_ports))
end

required_ports(c::InertialShaft) = _shaft_port_names(c.n_ports)

function validate(c::InertialShaft)
    c.J >= 0.0 || error("InertialShaft.J must be >= 0")
    c.damping >= 0.0 || error("InertialShaft.damping must be >= 0")
    c.n_ports >= 1 || error("InertialShaft.n_ports must be >= 1")
    return c
end
