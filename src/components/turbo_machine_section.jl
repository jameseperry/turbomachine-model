"""
Unified compressor/turbine component.

Fields:
- `mode`: operation mode, `:compressor` or `:turbine`.
- `performance_map`: referenced compressor/turbine performance map data.
- `eta_guess`: nominal isentropic efficiency guess used for initialization.
- `init`: named-tuple of component-specific initial condition values.
- `port_map`: component ports by name.
- `variable_list`: component variables.
"""
struct TurboMachineSection <: AbstractComponent
    mode::Symbol
    performance_map::PerformanceMap
    eta_guess::Float64
    init::NamedTuple
    port_map::Dict{Symbol,ComponentPort}
    variable_list::Vector{ComponentVariable}
end

function TurboMachineSection(;
    mode::Symbol,
    performance_map::PerformanceMap,
    eta_guess::Real,
    init::NamedTuple=NamedTuple(),
)
    mode in (:compressor, :turbine) ||
        error("TurboMachineSection.mode must be :compressor or :turbine")
    eta_f = Float64(eta_guess)
    0.0 < eta_f <= 1.0 ||
        error("TurboMachineSection.eta_guess must be in (0, 1]")

    port_map = Dict(
        :inlet => ComponentPort(
            :FluidPort;
            variable_ids=Dict(
                :pt => :inlet_pt,
                :ht => :inlet_ht,
                :composition => :inlet_composition,
                :mdot => :inlet_mdot,
            ),
        ),
        :outlet => ComponentPort(
            :FluidPort;
            variable_ids=Dict(
                :pt => :outlet_pt,
                :ht => :outlet_ht,
                :composition => :outlet_composition,
                :mdot => :outlet_mdot,
            ),
        ),
        :shaft => ComponentPort(
            :ShaftPort;
            variable_ids=Dict(
                :omega => :shaft_omega,
                :tau => :shaft_tau,
            ),
        ),
    )
    variable_list = ComponentVariable[
        ComponentVariable(:inlet_pt, :Pa),
        ComponentVariable(:inlet_ht, :J_per_kg),
        ComponentVariable(:inlet_composition, :species),
        ComponentVariable(:inlet_mdot, :kg_per_s),
        ComponentVariable(:outlet_pt, :Pa),
        ComponentVariable(:outlet_ht, :J_per_kg),
        ComponentVariable(:outlet_composition, :species),
        ComponentVariable(:outlet_mdot, :kg_per_s),
        ComponentVariable(:shaft_omega, :rad_per_s),
        ComponentVariable(:shaft_tau, :N_m),
    ]

    return TurboMachineSection(mode, performance_map, eta_f, init, port_map, variable_list)
end

ports(c::TurboMachineSection) = c.port_map

variables(c::TurboMachineSection) = c.variable_list
