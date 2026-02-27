"""
Analytic compressor performance map implementation.
"""

using TOML
import ....Utility: write_toml, read_toml

"""
Analytic compressor performance map on corrected coordinates.

This parametric model returns pressure ratio (`PR`) and adiabatic efficiency (`eta`)
as functions of corrected speed and corrected mass flow.
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

    # Reference stagnation conditions for corrected normalization
    Tt_ref::T = T(288.15)
    Pt_ref::T = T(101_325.0)
end

const _ANALYTIC_MAP_FIELDS = fieldnames(AnalyticCompressorPerformanceMap{Float64})

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

@inline _analytic_smooth_saturate(x::T, k::T) where {T<:Real} =
    x - _analytic_softplus(x - one(T), k) + _analytic_softplus(-x, k)

@inline function mdot_surge(map::AnalyticCompressorPerformanceMap, omega_corr::Real)
    N = omega_corr
    d = N - 1
    return map.ms0 + map.ms1 * d + map.ms2 * d * d
end

@inline function mdot_choke(map::AnalyticCompressorPerformanceMap, omega_corr::Real)
    N = omega_corr
    d = N - 1
    return map.mc0 + map.mc1 * d + map.mc2 * d * d
end

@inline function _normalized_flow_u(
    map::AnalyticCompressorPerformanceMap,
    omega_corr::Real,
    mdot_corr::Real,
)
    ms = mdot_surge(map, omega_corr)
    mc = mdot_choke(map, omega_corr)
    denom = max(mc - ms, eps(promote_type(typeof(ms), typeof(mc), typeof(mdot_corr))))
    x = (mdot_corr - ms) / denom
    return _analytic_smooth_saturate(x, map.sat_k)
end

@inline function _beta_bump_raw(u::T, alpha::T, beta::T) where {T<:Real}
    u_eps = clamp(u, T(1e-12), T(1 - 1e-12))
    return u_eps^(alpha - one(T)) * (one(T) - u_eps)^(beta - one(T))
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

@inline function _pr_shape(map::AnalyticCompressorPerformanceMap{T}, u::T) where {T<:Real}
    bump = _beta_bump_raw(u, map.alpha, map.beta) / max(_beta_bump_norm(map), eps(T))
    p_surge = one(T) - map.eps_s_pr * exp(-u / map.del_s_pr)
    p_choke = one(T) - map.eps_c_pr * exp(-(one(T) - u) / map.del_c_pr)
    return bump * p_surge * p_choke
end

@inline function _pr_speed_scale(map::AnalyticCompressorPerformanceMap{T}, omega_corr::Real) where {T<:Real}
    N = T(omega_corr)
    N_pos = max(N, zero(T))
    return map.Pi_max * (N_pos^map.pr_speed_exp)
end

"""
Evaluate a compressor map at corrected coordinates.

Returns named tuple `(PR, eta)` where:
- `PR` is total-pressure ratio (`Pt_out/Pt_in`)
- `eta` is adiabatic efficiency
"""
function compressor_performance_map(
    map::AnalyticCompressorPerformanceMap{T},
    omega_corr::Real,
    mdot_corr::Real,
) where {T<:Real}
    N = T(omega_corr)
    u = T(_normalized_flow_u(map, omega_corr, mdot_corr))

    Pi = _pr_speed_scale(map, omega_corr)
    PR = one(T) + Pi * _pr_shape(map, u)

    eta_peak = map.eta_max - map.eta_speed_quad * (N - one(T))^2
    u_star = map.u0 + map.u1 * (N - one(T))
    A = map.A0 * (one(T) + map.A_speed * (N - one(T))^2)

    eta = eta_peak - A * (u - u_star)^2
    eta -= map.Ds_eta * exp(-u / map.del_s_eta)
    eta -= map.Dc_eta * exp(-(one(T) - u) / map.del_c_eta)
    eta = clamp(eta, map.eta_min, map.eta_max_clip)

    return (PR=PR, eta=eta)
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

    for field in _ANALYTIC_MAP_FIELDS
        key = String(field)
        haskey(node, key) || error("missing TOML key $(key)")
    end

    kwargs = (; (field => T(node[String(field)]) for field in _ANALYTIC_MAP_FIELDS)...)
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
