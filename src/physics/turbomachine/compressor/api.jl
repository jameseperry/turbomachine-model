"""
Compressor performance map API.
"""

import ..AbstractCompressorPerformanceMap

"""
Evaluate a compressor map using physical operating conditions.

Returns:
- `PR`: total-pressure ratio (`Pt_out / Pt_in`)
- `eta`: adiabatic efficiency
- optional diagnostics (implementation-defined), for example map-coordinate values.
"""
performance_from_stagnation(
    map::AbstractCompressorPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
) = error("performance_from_stagnation not implemented for $(typeof(map))")

"""
Physical operating domain for a compressor map at a given inlet stagnation state.

Returns:
- `omega = (omega_min, omega_max)` in rad/s
- `mdot = (mdot_min, mdot_max)` in kg/s
- `mdot_flow_range = (surge, choke)` where each entry is `omega -> mdot`
"""
performance_map_domain(
    map::AbstractCompressorPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
) = error("performance_map_domain not implemented for $(typeof(map))")
