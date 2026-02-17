"""
Ideal-gas utility functions for converting between total (stagnation)
and static flow properties.

These helpers are intentionally simple and use scalar values.
"""

"""
Static temperature from total temperature, Mach number, and gamma.
"""
static_temperature_from_total(Tt::Real, mach::Real, gamma::Real) =
    Tt / (1 + 0.5 * (gamma - 1) * mach^2)

"""
Total temperature from static temperature, Mach number, and gamma.
"""
total_temperature_from_static(T::Real, mach::Real, gamma::Real) =
    T * (1 + 0.5 * (gamma - 1) * mach^2)

"""
Static pressure from total pressure, Mach number, and gamma.
"""
static_pressure_from_total(Pt::Real, mach::Real, gamma::Real) =
    Pt / (1 + 0.5 * (gamma - 1) * mach^2)^(gamma / (gamma - 1))

"""
Total pressure from static pressure, Mach number, and gamma.
"""
total_pressure_from_static(P::Real, mach::Real, gamma::Real) =
    P * (1 + 0.5 * (gamma - 1) * mach^2)^(gamma / (gamma - 1))

"""
Static enthalpy from total enthalpy and velocity.
"""
static_enthalpy_from_total(ht::Real, velocity::Real) =
    ht - 0.5 * velocity^2

"""
Total enthalpy from static enthalpy and velocity.
"""
total_enthalpy_from_static(h::Real, velocity::Real) =
    h + 0.5 * velocity^2
