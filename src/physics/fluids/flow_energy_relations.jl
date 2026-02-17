"""
Simple energy relations used in 1D flow models.
"""

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
