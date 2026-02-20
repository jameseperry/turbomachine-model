"""
Simple energy relations used in 1D flow models.
"""

"""
Flow velocity from mass flow rate, density, and cross-sectional area.

Uses `V = mdot / (density * A)`. Signed `mdot` returns signed velocity.
"""
function velocity_from_massflow(mdot::Real, density::Real, area::Real)
    density > 0 || error("density must be > 0")
    area > 0 || error("area must be > 0")
    return mdot / (density * area)
end

"""
Flow velocity from `(p, h, mdot, A)` using a density closure `rho_from_ph`.

`rho_from_ph` must accept `(p, h)` and return a positive density.
"""
function velocity_from_ph_mdot(
    p::Real,
    h::Real,
    mdot::Real,
    area::Real,
    rho_from_ph::Function,
)
    rho = rho_from_ph(p, h)
    return velocity_from_massflow(mdot, rho, area)
end

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
