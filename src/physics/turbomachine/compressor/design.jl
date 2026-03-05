"""
Physical-ish compressor design knobs for compiling to `CompressorSpec`.
"""

using TOML
import ....Utility: write_toml, read_toml

"""
Physically grounded compressor design inputs.

This is a sketch/schema for a physics-first design interface and is not yet wired
into `compile_compressor_spec` / `compile_compressor_map`.
"""
Base.@kwdef struct PhysicalCompressorDesign{T<:Real}
    # Working-fluid identifier used to resolve EOS properties (cp, gamma, rho, etc).
    composition::Symbol = :air

    # Design inlet total state (station 0->1), used with EOS to compute static state.
    pt_in_design::T = 101_325.0          # inlet total pressure [Pa]
    tt_in_design::T = 288.15             # inlet total temperature [K]

    # Design shaft speed (physical, unnormalized).
    omega_design::T = 1_000.0            # shaft speed [rad/s]

    # Inlet annulus geometry at rotor leading edge.
    #  Flow area is A1 = pi * (r_tip_inlet^2 - r_hub_inlet^2).
    r_tip_inlet::T = 0.25                # inlet tip radius [m]
    r_hub_inlet::T = 0.15                # inlet hub radius [m]

    # Number of aerodynamic stages (rotor+stator pairs).
    stage_count::Int = 8

    # Meanline design coefficients.
    #  phi = Vx / U, where Vx is axial absolute inlet velocity and
    #  U = omega*r is blade peripheral speed at the reference radius.
    #  psi = Delta h0 / U^2, where Delta h0 is the stage stagnation-enthalpy rise
    #  (h0_out - h0_in) and U^2 is the blade-speed specific-energy scale.
    phi_design::T = 0.55                 # flow coefficient at design point [-]
    psi_design::T = 0.45                 # stage loading coefficient at design point [-]
    reaction_design::T = 0.50            # stage reaction at design point [-]

    # Aerodynamic incidence/deviation targets near design.
    alpha1_design::T = 0.0               # inlet absolute flow angle [rad]
    beta1_metal_design::T = 0.0          # rotor inlet metal angle [rad]
    deviation_model_gain::T = 1.0        # rotor/stator deviation sensitivity gain [-]

    # Tip-speed and compressibility targets.
    tip_mach_design::T = 0.75            # rotor tip Mach at design point [-]
    max_tip_mach::T = 0.95               # allowable maximum rotor tip Mach [-]

    # Operability margins used to shape surge/choke boundaries.
    surge_margin_design::T = 0.20        # fractional standoff from surge at design [-]
    choke_margin_design::T = 0.10        # fractional standoff from choke at design [-]

    # Loss-model coefficients (lumped, meanline-level).
    profile_loss_coeff::T = 0.030        # blade profile-loss coefficient [-]
    incidence_loss_coeff::T = 0.020      # off-design incidence-loss coefficient [-]
    secondary_loss_coeff::T = 0.015      # secondary-flow loss coefficient [-]
    clearance_loss_coeff::T = 0.020      # tip-clearance loss coefficient [-]

    # Geometry/clearance quality.
    tip_clearance_fraction::T = 0.015    # tip clearance / blade height [-]
    surface_roughness_factor::T = 1.0    # roughness multiplier on profile losses [-]

    # Non-aerodynamic efficiencies used for shaft-power translation.
    mechanical_efficiency::T = 0.99      # bearing/seal mechanical efficiency [-]
    diffuser_pressure_recovery::T = 0.85 # diffuser pressure-recovery effectiveness [-]

    # Reynolds-number anchor for loss corrections (if enabled).
    reynolds_reference::T = 5.0e5        # reference Reynolds number [-]
end

Base.@kwdef struct CompressorDesign{T<:Real}
    # Compressor family/type assumption used by the empirical compilation logic.
    kind::Symbol = :axial

    # Number of aerodynamic stages (mainly relevant for axial machines).
    stage_count::Int = 8

    # Stage loading coefficient proxy; larger implies more work per stage (-).
    stage_loading::T = 0.55

    # Rotor tip Mach number at the design point; influences speed/choke tendencies (-).
    tip_mach_design::T = 0.50

    # Aggregate measure of diffusion aggressiveness in blade rows and passages (-).
    diffusion_aggressiveness::T = 0.50

    # Effective tip-clearance severity fraction; larger implies more leakage/losses (-).
    clearance_fraction::T = 0.20

    # Diffuser quality factor capturing pressure-recovery quality (-).
    diffuser_quality::T = 0.60

    # Variable-geometry authority/usage level (e.g., IGV/VSV effectiveness) (-).
    variable_geometry::T = 0.30

    # Reynolds-number quality factor capturing viscous-performance quality (-).
    reynolds_quality::T = 0.80
end

const _COMPRESSOR_DESIGN_REAL_FIELDS = (
    :stage_loading,
    :tip_mach_design,
    :diffusion_aggressiveness,
    :clearance_fraction,
    :diffuser_quality,
    :variable_geometry,
    :reynolds_quality,
)

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

function write_toml(
    design::CompressorDesign,
    path::AbstractString;
    group::AbstractString="compressor_design",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_design"
    node["format_version"] = 1
    node["kind"] = String(design.kind)
    node["stage_count"] = design.stage_count
    for field in _COMPRESSOR_DESIGN_REAL_FIELDS
        node[String(field)] = Float64(getfield(design, field))
    end

    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{CompressorDesign{T}},
    path::AbstractString;
    group::AbstractString="compressor_design",
) where {T<:Real}
    data = TOML.parsefile(path)
    node = _find_group(data, group)

    haskey(node, "kind") || error("missing TOML key kind")
    haskey(node, "stage_count") || error("missing TOML key stage_count")
    for field in _COMPRESSOR_DESIGN_REAL_FIELDS
        key = String(field)
        haskey(node, key) || error("missing TOML key $(key)")
    end

    kind = Symbol(String(node["kind"]))
    stage_count = Int(node["stage_count"])
    kwargs = (; (field => T(node[String(field)]) for field in _COMPRESSOR_DESIGN_REAL_FIELDS)...)
    return CompressorDesign{T}(;
        kind=kind,
        stage_count=stage_count,
        kwargs...,
    )
end

function read_toml(
    ::Type{CompressorDesign},
    path::AbstractString;
    group::AbstractString="compressor_design",
)
    return read_toml(CompressorDesign{Float64}, path; group=group)
end

"""Demo compressor design for development/testing."""
function demo_compressor_design()
    return CompressorDesign{Float64}()
end
