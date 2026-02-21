"""
Lumped plenum volume component with one inlet and one outlet fluid port.

Fields:
- `volume`: plenum control-volume size in m^3.
- `init`: named-tuple of component-specific initial condition values.
"""
struct Plenum <: AbstractComponent
    volume::Float64
    init::NamedTuple
end

function Plenum(;
    volume::Real,
    init::NamedTuple=NamedTuple(),
)
    comp = Plenum(Float64(volume), init)
    validate(comp)
    return comp
end

function port_specs(::Plenum)
    Dict(
        :inlet => FLUID_PORT,
        :outlet => FLUID_PORT,
    )
end

required_ports(::Plenum) = [:inlet, :outlet]

function validate(c::Plenum)
    c.volume > 0.0 || error("Plenum.volume must be > 0")
    return c
end
