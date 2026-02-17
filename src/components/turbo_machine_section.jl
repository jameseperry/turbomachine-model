"""
Unified compressor/turbine component.

Fields:
- `mode`: operation mode, `:compressor` or `:turbine`.
- `performance_map`: referenced compressor/turbine performance map data.
- `eta_guess`: nominal isentropic efficiency guess used for initialization.
- `init`: named-tuple of component-specific initial condition values.
"""
struct TurboMachineSection <: AbstractComponent
    mode::Symbol
    performance_map::PerformanceMap
    eta_guess::Float64
    init::NamedTuple
end

function TurboMachineSection(;
    mode::Symbol,
    performance_map::PerformanceMap,
    eta_guess::Real,
    init::NamedTuple=NamedTuple(),
)
    comp = TurboMachineSection(mode, performance_map, Float64(eta_guess), init)
    validate(comp)
    return comp
end

function port_specs(::TurboMachineSection)
    Dict(
        :inlet => FLUID_PORT,
        :outlet => FLUID_PORT,
        :shaft => SHAFT_PORT,
    )
end

required_ports(::TurboMachineSection) = [:inlet, :outlet, :shaft]

function validate(c::TurboMachineSection)
    c.mode in (:compressor, :turbine) ||
        error("TurboMachineSection.mode must be :compressor or :turbine")
    0.0 < c.eta_guess <= 1.0 ||
        error("TurboMachineSection.eta_guess must be in (0, 1]")
    return c
end
