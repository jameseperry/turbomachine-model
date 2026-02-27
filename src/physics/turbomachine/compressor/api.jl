"""
Compressor performance map API.
"""

abstract type AbstractCompressorPerformanceMap end

"""
Recommended corrected-coordinate plotting/usage domain.

Returns:
- `omega_corr = (omega_min, omega_max)`
- `mdot_corr = (mdot_min, mdot_max)`
- `mdot_corr_flow_range = (surge, choke)` where each entry is a function
  `omega_corr -> mdot_corr`
"""
performance_map_domain(map::AbstractCompressorPerformanceMap) =
    error("performance_map_domain not implemented for $(typeof(map))")
