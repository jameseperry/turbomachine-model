"""
Turbine performance map API.
"""

abstract type AbstractTurbinePerformanceMap end

"""
Evaluate a turbine map in corrected coordinates.

Inputs:
- `omega_corr`: corrected shaft speed coordinate.
- `pr_turb`: turbine pressure ratio coordinate (`Pt_in / Pt_out`).

Returns:
- `mdot_corr`: corrected mass flow.
- `eta`: turbine adiabatic efficiency.
"""
function turbine_performance_map(
    map::AbstractTurbinePerformanceMap,
    omega_corr::Real,
    pr_turb::Real,
)
    error("turbine_performance_map not implemented for $(typeof(map))")
end

"""
Evaluate a turbine map from physical stagnation-state operating conditions.

Inputs:
- `omega`: shaft speed [rad/s].
- `Pt_in`: inlet total pressure [Pa].
- `Pt_out`: outlet total pressure [Pa].
- `Tt_in`: inlet total temperature [K].

Returns at minimum:
- `PR_turb = Pt_in / Pt_out`
- `mdot`: physical mass flow [kg/s]
- `eta`: turbine adiabatic efficiency
- optional map-coordinate diagnostics (implementation-defined), for example
  `omega_corr` and `mdot_corr`.
"""
function turbine_performance_map_from_stagnation(
    map::AbstractTurbinePerformanceMap,
    omega::Real,
    Pt_in::Real,
    Pt_out::Real,
    Tt_in::Real,
)
    error("turbine_performance_map_from_stagnation not implemented for $(typeof(map))")
end

"""
Recommended corrected-coordinate plotting/usage domain for a turbine map.

Returns:
- `omega_corr = (omega_min, omega_max)`
- `pr_turb = (pr_min, pr_max)`
"""
performance_map_domain(map::AbstractTurbinePerformanceMap) =
    error("performance_map_domain not implemented for $(typeof(map))")
