"""
Simple energy relations used in 1D flow models.
"""

@inline _flow_primal_value(x::Real) = hasfield(typeof(x), :value) ? getfield(x, :value) : x
@inline _flow_positive_finite(x::Real) = isfinite(_flow_primal_value(x)) && _flow_primal_value(x) > 0

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
Solve for axial velocity from stagnation state, mass flow, area, and tangential
velocity.

This inverts the continuity relation

`mdot = rho(p, h) * Vx * area`

where the static state `(p, h)` is implied by the stagnation state
`(pt_t, ht_t)` and the full velocity vector `(Vx, Vtheta)`.
"""
function velocity_from_stagnation_massflow(
    eos::AbstractEOS,
    mdot::Real,
    pt_t::Real,
    ht_t::Real,
    area::Real,
    Vtheta::Real;
    prefer::Symbol=:low,
    prior_roots::AbstractVector{<:Real}=Float64[],
)
    mdot_f = Float64(mdot)
    pt_t_f = Float64(pt_t)
    ht_t_f = Float64(ht_t)
    area_f = Float64(area)
    Vtheta_f = Float64(Vtheta)

    (_flow_positive_finite(mdot_f) &&
     _flow_positive_finite(pt_t_f) &&
     _flow_positive_finite(ht_t_f) &&
     _flow_positive_finite(area_f)) || return (converged=false, Vx=NaN, roots=Float64[])

    s_t = try
        entropy(eos, pt_t_f, ht_t_f)
    catch
        NaN
    end
    isfinite(_flow_primal_value(s_t)) || return (converged=false, Vx=NaN, roots=Float64[])

    kinetic_margin = 2 * ht_t_f - Vtheta_f^2
    kinetic_margin > 0.0 || return (converged=false, Vx=NaN, roots=Float64[])
    x_hi = sqrt(kinetic_margin) * (1 - 1e-10)

    function residual(Vx)
        Vx > 0.0 || return NaN
        h = static_enthalpy_from_total(ht_t_f, hypot(Vx, Vtheta_f))
        p = try
            pressure_from_enthalpy_entropy(eos, h, s_t)
        catch
            NaN
        end
        (_flow_positive_finite(p) && _flow_positive_finite(h)) || return NaN
        rho = try
            density(eos, p, h)
        catch
            NaN
        end
        _flow_positive_finite(rho) || return NaN
        return mdot_f - rho * Vx * area_f
    end

    roots = bracket_bisect_roots(
        residual,
        (1e-10, x_hi);
        n_scan=401,
        root_tol=1e-9,
        max_bisect_iters=80,
        prior_roots=prior_roots,
        dedupe_atol=1e-8,
    )
    isempty(roots) && return (converged=false, Vx=NaN, roots=roots)
    return (converged=true, Vx=(prefer == :high ? last(roots) : first(roots)), roots=roots)
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
