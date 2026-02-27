"""
Analytic compressor performance map implementation.
"""

using TOML
import ....Utility: write_toml, read_toml

"""
Analytic compressor performance map on physical shaft speed and corrected flow.

This parametric model returns pressure ratio (`PR`) and adiabatic efficiency (`eta`)
as functions of physical shaft speed (`omega`, rad/s) and corrected mass flow.
"""
Base.@kwdef struct AnalyticCompressorPerformanceMap{T<:Real} <: AbstractCompressorPerformanceMap
    # Surge boundary: mdot_surge(N) = ms0 + ms1*(N-1) + ms2*(N-1)^2
    ms0::T = T(0.55)
    ms1::T = T(0.10)
    ms2::T = T(0.00)

    # Choke boundary: mdot_choke(N) = mc0 + mc1*(N-1) + mc2*(N-1)^2
    mc0::T = T(1.10)
    mc1::T = T(0.12)
    mc2::T = T(0.00)

    # Smooth saturation sharpness (larger => sharper)
    sat_k::T = T(40.0)

    # Pressure-ratio model: PR = 1 + Pi(N) * S_pi(u)
    Pi_max::T = T(1.6)
    pr_speed_exp::T = T(2.0)

    # Beta-bump shape parameters
    alpha::T = T(2.2)
    beta::T = T(2.6)

    # Boundary penalties for PR
    eps_s_pr::T = T(0.35)
    del_s_pr::T = T(0.10)
    eps_c_pr::T = T(0.18)
    del_c_pr::T = T(0.12)

    # Efficiency model
    eta_max::T = T(0.88)
    eta_speed_quad::T = T(0.08)
    u0::T = T(0.52)
    u1::T = T(0.05)
    A0::T = T(0.18)
    A_speed::T = T(0.40)
    Ds_eta::T = T(0.05)
    del_s_eta::T = T(0.12)
    Dc_eta::T = T(0.03)
    del_c_eta::T = T(0.10)

    # Output clamps
    eta_min::T = T(0.50)
    eta_max_clip::T = T(0.92)

    # Reference stagnation conditions for flow correction
    Tt_ref::T = T(288.15)
    Pt_ref::T = T(101_325.0)

    # Internal speed normalization and recommended physical speed domain.
    omega_ref::T = T(1_000.0)
    omega_norm_min::T = T(0.6)
    omega_norm_max::T = T(1.0)
end

const _ANALYTIC_MAP_FIELDS = fieldnames(AnalyticCompressorPerformanceMap{Float64})

"""Analytic maps use physical shaft speed directly (`omega` in rad/s)."""
corrected_speed(omega::Real, Tt_in::Real, map::AnalyticCompressorPerformanceMap) = omega

@inline function _omega_norm(map::AnalyticCompressorPerformanceMap{T}, omega::Real) where {T<:Real}
    U = promote_type(typeof(omega), T)
    omega_ref = max(U(map.omega_ref), eps(U))
    return U(omega) / omega_ref
end

@inline function _analytic_softplus(x::T, k::T) where {T<:Real}
    z = k * x
    if z > T(40)
        return x
    elseif z < T(-40)
        return zero(T)
    else
        return log1p(exp(z)) / k
    end
end

@inline function _analytic_softplus(x::Real, k::Real)
    T = promote_type(typeof(x), typeof(k))
    return _analytic_softplus(T(x), T(k))
end

@inline function _analytic_smooth_saturate(x::Real, k::Real)
    T = promote_type(typeof(x), typeof(k))
    xT = T(x)
    kT = T(k)
    return xT - _analytic_softplus(xT - one(T), kT) + _analytic_softplus(-xT, kT)
end

@inline function mdot_surge(map::AnalyticCompressorPerformanceMap, omega::Real)
    N = _omega_norm(map, omega)
    d = N - one(N)
    return map.ms0 + map.ms1 * d + map.ms2 * d * d
end

@inline function mdot_choke(map::AnalyticCompressorPerformanceMap, omega::Real)
    N = _omega_norm(map, omega)
    d = N - one(N)
    return map.mc0 + map.mc1 * d + map.mc2 * d * d
end

@inline function _normalized_flow_u(
    map::AnalyticCompressorPerformanceMap,
    omega::Real,
    mdot_corr::Real,
)
    ms = mdot_surge(map, omega)
    mc = mdot_choke(map, omega)
    delta_m = mc - ms
    denom = abs(delta_m) > 1e-12 ? delta_m : oftype(delta_m, 1e-12)
    x = (mdot_corr - ms) / denom
    return _analytic_smooth_saturate(x, map.sat_k)
end

@inline function _beta_bump_raw(u::Real, alpha::Real, beta::Real)
    T = promote_type(typeof(u), typeof(alpha), typeof(beta))
    uT = T(u)
    alphaT = T(alpha)
    betaT = T(beta)
    u_eps = clamp(uT, T(1e-12), T(1 - 1e-12))
    return u_eps^(alphaT - one(T)) * (one(T) - u_eps)^(betaT - one(T))
end

function _beta_bump_norm(map::AnalyticCompressorPerformanceMap{T}) where {T<:Real}
    alpha = map.alpha
    beta = map.beta

    if alpha > one(T) && beta > one(T)
        u_star = (alpha - one(T)) / (alpha + beta - T(2))
        return _beta_bump_raw(u_star, alpha, beta)
    end

    us = range(T(1e-6), T(1 - 1e-6), length=1001)
    return maximum(_beta_bump_raw(u, alpha, beta) for u in us)
end

@inline function _pr_shape(map::AnalyticCompressorPerformanceMap{T}, u::Real) where {T<:Real}
    U = promote_type(typeof(u), T)
    uU = U(u)
    bump = _beta_bump_raw(uU, U(map.alpha), U(map.beta)) /
           max(U(_beta_bump_norm(map)), eps(U))
    p_surge = one(U) - U(map.eps_s_pr) * exp(-uU / U(map.del_s_pr))
    p_choke = one(U) - U(map.eps_c_pr) * exp(-(one(U) - uU) / U(map.del_c_pr))
    return bump * p_surge * p_choke
end

@inline function _pr_speed_scale(map::AnalyticCompressorPerformanceMap{T}, omega::Real) where {T<:Real}
    U = promote_type(typeof(omega), T)
    N = _omega_norm(map, omega)
    N_pos = max(N, zero(U))
    return U(map.Pi_max) * (N_pos^U(map.pr_speed_exp))
end

"""
Evaluate a compressor map at physical shaft speed and corrected mass flow.

Returns named tuple `(PR, eta)` where:
- `PR` is total-pressure ratio (`Pt_out/Pt_in`)
- `eta` is adiabatic efficiency
"""
function compressor_performance_map(
    map::AnalyticCompressorPerformanceMap{T},
    omega::Real,
    mdot_corr::Real,
) where {T<:Real}
    U = promote_type(typeof(omega), typeof(mdot_corr), T)
    N = _omega_norm(map, omega)
    u = U(_normalized_flow_u(map, omega, mdot_corr))

    Pi = _pr_speed_scale(map, omega)
    PR = one(U) + Pi * _pr_shape(map, u)

    eta_peak = U(map.eta_max) - U(map.eta_speed_quad) * (N - one(U))^2
    u_star = U(map.u0) + U(map.u1) * (N - one(U))
    A = U(map.A0) * (one(U) + U(map.A_speed) * (N - one(U))^2)

    eta = eta_peak - A * (u - u_star)^2
    eta -= U(map.Ds_eta) * exp(-u / U(map.del_s_eta))
    eta -= U(map.Dc_eta) * exp(-(one(U) - u) / U(map.del_c_eta))
    eta = clamp(eta, U(map.eta_min), U(map.eta_max_clip))

    return (PR=PR, eta=eta)
end

function performance_map_domain(map::AnalyticCompressorPerformanceMap)
    omega_min = min(map.omega_norm_min, map.omega_norm_max) * map.omega_ref
    omega_max = max(map.omega_norm_min, map.omega_norm_max) * map.omega_ref
    smin = mdot_surge(map, omega_min)
    smax = mdot_surge(map, omega_max)
    cmin = mdot_choke(map, omega_min)
    cmax = mdot_choke(map, omega_max)
    return (
        omega_corr=(omega_min, omega_max),
        mdot_corr=(min(smin, smax), max(cmin, cmax)),
        mdot_corr_flow_range=(
            surge=(omega_corr -> mdot_surge(map, omega_corr)),
            choke=(omega_corr -> mdot_choke(map, omega_corr)),
        ),
    )
end

function write_toml(
    map::AnalyticCompressorPerformanceMap,
    path::AbstractString;
    group::AbstractString="compressor_analytic_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_analytic_performance_map"
    node["format_version"] = 1

    for field in _ANALYTIC_MAP_FIELDS
        node[String(field)] = Float64(getfield(map, field))
    end

    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{AnalyticCompressorPerformanceMap{T}},
    path::AbstractString;
    group::AbstractString="compressor_analytic_map",
) where {T<:Real}
    data = TOML.parsefile(path)
    node = _find_group(data, group)

    defaults = AnalyticCompressorPerformanceMap{T}()
    kwargs = (; (
        field => (haskey(node, String(field)) ? T(node[String(field)]) : getfield(defaults, field))
        for field in _ANALYTIC_MAP_FIELDS
    )...)
    return AnalyticCompressorPerformanceMap{T}(; kwargs...)
end

function read_toml(
    ::Type{AnalyticCompressorPerformanceMap},
    path::AbstractString;
    group::AbstractString="compressor_analytic_map",
)
    return read_toml(AnalyticCompressorPerformanceMap{Float64}, path; group=group)
end

"""Demo analytic compressor map for development/testing."""
function demo_analytic_compressor_performance_map()
    return AnalyticCompressorPerformanceMap{Float64}()
end
