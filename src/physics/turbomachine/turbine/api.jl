"""
Turbine performance map API.
"""

abstract type AbstractTurbinePerformanceMap end

"""
Recommended corrected-coordinate plotting/usage domain.

Returns:
- `omega_corr = (omega_min, omega_max)`
- `pr_turb = (pr_min, pr_max)`
"""
performance_map_domain(map::AbstractTurbinePerformanceMap) =
    error("performance_map_domain not implemented for $(typeof(map))")
