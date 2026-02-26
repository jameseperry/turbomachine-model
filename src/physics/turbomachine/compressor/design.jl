"""
Physical-ish compressor design knobs for compiling to `CompressorSpec`.
"""

Base.@kwdef struct CompressorDesign{T<:Real}
    kind::Symbol = :axial
    stage_count::Int = 8
    stage_loading::T = 0.55
    tip_mach_design::T = 0.50
    diffusion_aggressiveness::T = 0.50
    clearance_fraction::T = 0.20
    diffuser_quality::T = 0.60
    variable_geometry::T = 0.30
    reynolds_quality::T = 0.80
end

"""
    compile_compressor_spec(design::CompressorDesign) -> CompressorSpec

Compile physical-ish design knobs to a behavioral `CompressorSpec`.
"""
function compile_compressor_spec(design::CompressorDesign{T}) where {T<:Real}
    kind = design.kind
    is_axial = (kind == :axial)
    is_cent = (kind == :centrifugal)
    is_axial || is_cent || error("design.kind must be :axial or :centrifugal")

    psi = _clamp01(design.stage_loading)
    M = _clamp01(design.tip_mach_design)
    D = _clamp01(design.diffusion_aggressiveness)
    C = _clamp01(design.clearance_fraction)
    Qd = _clamp01(design.diffuser_quality)
    VG = _clamp01(design.variable_geometry)
    Rq = _clamp01(design.reynolds_quality)

    nstg = is_axial ? max(design.stage_count, 1) : 1

    base_pr = is_axial ? T(2.0) : T(2.4)
    stg_gain = is_axial ? T(0.12) * log(T(nstg)) : T(0.0)
    load_gain = _lerp(T(0.3), T(2.0), psi)
    diff_gain = _lerp(T(0.0), T(0.6), D)
    clear_loss = _lerp(T(0.0), T(0.6), C)
    pr_design = max(base_pr + stg_gain + load_gain + diff_gain - clear_loss, T(1.1))

    base_eta = is_axial ? T(0.89) : T(0.85)
    qdiff = is_cent ? _lerp(T(-0.03), T(0.03), Qd) : _lerp(T(-0.01), T(0.01), Qd)
    qrey = _lerp(T(-0.04), T(0.02), Rq)
    qclr = -_lerp(T(0.00), T(0.06), C)
    qdif = -_lerp(T(0.00), T(0.03), D)
    qvg = is_axial ? _lerp(T(0.00), T(0.015), VG) : _lerp(T(0.00), T(0.005), VG)
    eta_design = clamp(base_eta + qdiff + qrey + qclr + qdif + qvg, T(0.70), T(0.92))

    base_range = is_axial ? T(0.60) : T(0.38)
    flow_range = base_range
    flow_range -= _lerp(T(0.00), T(0.25), psi)
    flow_range -= _lerp(T(0.00), T(0.15), D)
    flow_range -= is_cent ? _lerp(T(0.00), T(0.10), M) : _lerp(T(0.00), T(0.05), M)
    flow_range += is_axial ? _lerp(T(0.00), T(0.15), VG) : _lerp(T(0.00), T(0.05), VG)
    flow_range -= _lerp(T(0.00), T(0.05), C)
    flow_range = clamp(flow_range, T(0.15), T(0.85))

    surge_margin = T(0.55)
    surge_margin += is_axial ? _lerp(T(0.00), T(0.30), VG) : _lerp(T(0.00), T(0.15), VG)
    surge_margin -= _lerp(T(0.00), T(0.25), psi)
    surge_margin -= _lerp(T(0.00), T(0.20), D)
    surge_margin -= is_cent ? T(0.05) : T(0.00)
    surge_margin = _clamp01(surge_margin)

    choke_sharpness = T(0.45)
    choke_sharpness += is_cent ? _lerp(T(0.00), T(0.35), M) : _lerp(T(0.00), T(0.20), M)
    choke_sharpness += is_cent ? _lerp(T(0.00), T(0.25), one(T) - Qd) :
                                 _lerp(T(0.00), T(0.10), one(T) - Qd)
    choke_sharpness += _lerp(T(0.00), T(0.10), psi)
    choke_sharpness = _clamp01(choke_sharpness)

    speed_sensitivity = is_axial ? T(0.45) : T(0.65)
    speed_sensitivity += _lerp(T(-0.10), T(0.25), M)
    speed_sensitivity += _lerp(T(0.00), T(0.05), psi)
    speed_sensitivity = _clamp01(speed_sensitivity)

    return CompressorSpec{T}(
        pr_design=pr_design,
        eta_design=eta_design,
        flow_range=flow_range,
        surge_margin=surge_margin,
        choke_sharpness=choke_sharpness,
        speed_sensitivity=speed_sensitivity,
    )
end

"""
    compile_compressor_map(design::CompressorDesign) -> AnalyticCompressorPerformanceMap

Compile a design directly to an analytic compressor performance map.
"""
function compile_compressor_map(design::CompressorDesign{T}) where {T<:Real}
    spec = compile_compressor_spec(design)
    return compile_compressor_map(spec; kind=design.kind)
end
