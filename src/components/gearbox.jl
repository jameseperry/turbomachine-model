"""
Idealized two-port gearbox component.

Fields:
- `ratio`: speed ratio `omega_out / omega_in`.
- `efficiency`: mechanical efficiency in `(0, 1]`.
- `init`: named-tuple of component-specific initial condition values.
- `port_map`: component ports by name.
- `variable_list`: component variables.
"""
struct Gearbox <: AbstractComponent
    ratio::Float64
    efficiency::Float64
    init::NamedTuple
    port_map::Dict{Symbol,ComponentPort}
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
    port_map = Dict(
        :input => ComponentPort(
            :ShaftPort;
            variable_ids=Dict(
                :omega => :input_omega,
                :tau => :input_tau,
            ),
        ),
        :output => ComponentPort(
            :ShaftPort;
            variable_ids=Dict(
                :omega => :output_omega,
                :tau => :output_tau,
            ),
        ),
    )

    variable_list = ComponentVariable[
        ComponentVariable(:input_omega, :rad_per_s),
        ComponentVariable(:input_tau, :N_m),
        ComponentVariable(:output_omega, :rad_per_s),
        ComponentVariable(:output_tau, :N_m),
    ]

    return Gearbox(ratio_f, efficiency_f, init, port_map, variable_list)
end

ports(c::Gearbox) = c.port_map

variables(c::Gearbox) = c.variable_list
