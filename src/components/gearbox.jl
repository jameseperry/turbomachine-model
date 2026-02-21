"""
Idealized two-port gearbox component.

Fields:
- `ratio`: speed ratio `omega_out / omega_in`.
- `efficiency`: mechanical efficiency in `(0, 1]`.
- `init`: named-tuple of component-specific initial condition values.
"""
struct Gearbox <: AbstractComponent
    ratio::Float64
    efficiency::Float64
    init::NamedTuple
end

function Gearbox(;
    ratio::Real,
    efficiency::Real=1.0,
    init::NamedTuple=NamedTuple(),
)
    comp = Gearbox(Float64(ratio), Float64(efficiency), init)
    validate(comp)
    return comp
end

function port_specs(::Gearbox)
    Dict(
        :input => SHAFT_PORT,
        :output => SHAFT_PORT,
    )
end

required_ports(::Gearbox) = [:input, :output]

function validate(c::Gearbox)
    c.ratio > 0.0 || error("Gearbox.ratio must be > 0")
    0.0 < c.efficiency <= 1.0 || error("Gearbox.efficiency must be in (0, 1]")
    return c
end
