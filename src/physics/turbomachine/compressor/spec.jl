"""
Behavioral compressor envelope for compiling to `AnalyticCompressorPerformanceMap`.
"""

using TOML
import ....Utility: write_toml, read_toml

Base.@kwdef struct CompressorSpec{T<:Real}
    pr_design::T = 3.0
    eta_design::T = 0.88
    flow_range::T = 0.55
    surge_margin::T = 0.60
    choke_sharpness::T = 0.40
    speed_sensitivity::T = 0.50
end

const _COMPRESSOR_SPEC_FIELDS = fieldnames(CompressorSpec{Float64})

@inline _clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))
@inline _lerp(a::T, b::T, t::T) where {T<:Real} = a + (b - a) * t

"""
    compile_compressor_map(spec::CompressorSpec; kind=:axial) -> AnalyticCompressorPerformanceMap

Compile behavioral envelope parameters into analytic map coefficients.
"""
function compile_compressor_map(spec::CompressorSpec{T}; kind::Symbol=:axial) where {T<:Real}
    is_axial = (kind == :axial)
    is_cent = (kind == :centrifugal)
    is_axial || is_cent || error("kind must be :axial or :centrifugal")

    pr_design = max(spec.pr_design, T(1.05))
    eta_design = clamp(spec.eta_design, T(0.65), T(0.92))
    flow_range = clamp(spec.flow_range, T(0.15), T(0.85))
    surge_margin = _clamp01(spec.surge_margin)
    choke_sharpness = _clamp01(spec.choke_sharpness)
    speed_sensitivity = _clamp01(spec.speed_sensitivity)

    ms0 = T(1.0) - flow_range / T(2)
    mc0 = T(1.0) + flow_range / T(2)

    ms1 = is_axial ? _lerp(T(0.06), T(0.14), T(1) - surge_margin) :
                     _lerp(T(0.05), T(0.12), T(1) - surge_margin)
    mc1 = is_axial ? _lerp(T(0.08), T(0.16), choke_sharpness) :
                     _lerp(T(0.10), T(0.20), choke_sharpness)
    ms2 = T(0.0)
    mc2 = T(0.0)

    Pi_max = pr_design - T(1.0)
    pmin, pmax = is_axial ? (T(1.8), T(2.2)) : (T(2.4), T(3.0))
    pr_speed_exp = _lerp(pmin, pmax, speed_sensitivity)

    peak_u = is_axial ? _lerp(T(0.55), T(0.48), choke_sharpness) :
                        _lerp(T(0.52), T(0.44), choke_sharpness)
    kappa_base = is_axial ? T(6.0) : T(8.0)
    kappa = kappa_base +
            _lerp(T(0.0), T(6.0), T(1) - surge_margin) +
            _lerp(T(0.0), T(4.0), choke_sharpness)

    alpha = peak_u * (kappa - T(2)) + T(1)
    beta = (T(1) - peak_u) * (kappa - T(2)) + T(1)

    eps_s_pr = is_axial ? _lerp(T(0.50), T(0.20), surge_margin) :
                          _lerp(T(0.60), T(0.25), surge_margin)
    del_s_pr = is_axial ? _lerp(T(0.06), T(0.14), surge_margin) :
                          _lerp(T(0.05), T(0.12), surge_margin)

    eps_c_pr = is_axial ? _lerp(T(0.12), T(0.35), choke_sharpness) :
                          _lerp(T(0.20), T(0.45), choke_sharpness)
    del_c_pr = is_axial ? _lerp(T(0.14), T(0.06), choke_sharpness) :
                          _lerp(T(0.10), T(0.05), choke_sharpness)

    eta_max = eta_design
    eta_speed_quad = is_axial ? T(0.06) + T(0.04) * speed_sensitivity :
                                T(0.08) + T(0.06) * speed_sensitivity

    u0 = is_axial ? T(0.52) : T(0.48)
    u1 = is_axial ? _lerp(T(0.02), T(0.08), speed_sensitivity) :
                    _lerp(T(0.05), T(0.12), speed_sensitivity)

    width_base = is_axial ? T(0.18) : T(0.26)
    A0 = clamp(width_base * (T(0.55) / flow_range), T(0.10), T(0.40))
    A_speed = is_axial ? T(0.35) : T(0.55)

    Ds_eta = is_axial ? _lerp(T(0.07), T(0.03), surge_margin) :
                        _lerp(T(0.09), T(0.04), surge_margin)
    del_s_eta = is_axial ? _lerp(T(0.08), T(0.14), surge_margin) :
                           _lerp(T(0.07), T(0.12), surge_margin)

    Dc_eta = is_axial ? _lerp(T(0.02), T(0.06), choke_sharpness) :
                        _lerp(T(0.03), T(0.08), choke_sharpness)
    del_c_eta = is_axial ? _lerp(T(0.12), T(0.07), choke_sharpness) :
                           _lerp(T(0.10), T(0.05), choke_sharpness)

    sat_k = is_axial ? T(50.0) : T(40.0)

    return AnalyticCompressorPerformanceMap{T}(
        ms0=ms0, ms1=ms1, ms2=ms2,
        mc0=mc0, mc1=mc1, mc2=mc2,
        sat_k=sat_k,
        Pi_max=Pi_max, pr_speed_exp=pr_speed_exp,
        alpha=alpha, beta=beta,
        eps_s_pr=eps_s_pr, del_s_pr=del_s_pr,
        eps_c_pr=eps_c_pr, del_c_pr=del_c_pr,
        eta_max=eta_max, eta_speed_quad=eta_speed_quad,
        u0=u0, u1=u1,
        A0=A0, A_speed=A_speed,
        Ds_eta=Ds_eta, del_s_eta=del_s_eta,
        Dc_eta=Dc_eta, del_c_eta=del_c_eta,
        eta_min=T(0.50), eta_max_clip=T(0.92),
        Tt_ref=T(288.15), Pt_ref=T(101_325.0),
    )
end

function write_toml(
    spec::CompressorSpec,
    path::AbstractString;
    group::AbstractString="compressor_spec",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_spec"
    node["format_version"] = 1

    for field in _COMPRESSOR_SPEC_FIELDS
        node[String(field)] = Float64(getfield(spec, field))
    end

    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{CompressorSpec{T}},
    path::AbstractString;
    group::AbstractString="compressor_spec",
) where {T<:Real}
    data = TOML.parsefile(path)
    node = _find_group(data, group)

    for field in _COMPRESSOR_SPEC_FIELDS
        key = String(field)
        haskey(node, key) || error("missing TOML key $(key)")
    end

    kwargs = (; (field => T(node[String(field)]) for field in _COMPRESSOR_SPEC_FIELDS)...)
    return CompressorSpec{T}(; kwargs...)
end

function read_toml(
    ::Type{CompressorSpec},
    path::AbstractString;
    group::AbstractString="compressor_spec",
)
    return read_toml(CompressorSpec{Float64}, path; group=group)
end

"""Demo compressor spec for development/testing."""
function demo_compressor_spec()
    return CompressorSpec{Float64}()
end
