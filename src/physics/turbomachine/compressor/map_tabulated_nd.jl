"""
Tabulated compressor performance map implementation (non-dimensional form).
"""

import ..NondimensionalPerformanceMap
const NonDimensionalTabulatedCompressorPerformanceMap = NondimensionalPerformanceMap

_m_tip_grid(map::NondimensionalPerformanceMap) = table_xgrid(map.pr_map)
_phi_in_grid(map::NondimensionalPerformanceMap) = table_ygrid(map.pr_map)
_pr_table(map::NondimensionalPerformanceMap) = table_values(map.pr_map)
_eta_table(map::NondimensionalPerformanceMap) = table_values(map.eta_map)
_interpolation_kind(map::NondimensionalPerformanceMap) = table_interpolation(map.pr_map)

function _interp_linear_1d_clamped(xgrid::AbstractVector{<:Real}, ygrid::AbstractVector{<:Real}, x::Real)
    x1 = first(xgrid)
    x2 = last(xgrid)
    xc = clamp(x, x1, x2)
    i_hi = searchsortedfirst(xgrid, xc)
    i_hi <= 1 && return ygrid[1]
    i_hi > length(xgrid) && return ygrid[end]
    i_lo = i_hi - 1
    xl = xgrid[i_lo]
    xr = xgrid[i_hi]
    t = (xc - xl) / (xr - xl)
    return (1 - t) * ygrid[i_lo] + t * ygrid[i_hi]
end

"""Map-coordinate speed for non-dimensional maps (`M_tip`)."""
function _map_speed_coordinate_from_stagnation(
    map::NondimensionalPerformanceMap,
    omega::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    Tt_in > 0 || error("Tt_in must be > 0")
    a0_in = sqrt(map.gamma * map.gas_constant * Tt_in)
    return omega * map.tip_radius_inlet / a0_in
end

"""Map-coordinate flow for non-dimensional maps (`phi_in`)."""
function _map_flow_coordinate_from_stagnation(
    map::NondimensionalPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    Tt_in > 0 || error("Tt_in must be > 0")
    Pt_in > 0 || error("Pt_in must be > 0")
    omega > 0 || error("omega must be > 0")
    rho0_in = Pt_in / (map.gas_constant * Tt_in)
    U_m1 = omega * map.mean_radius_inlet
    return mdot / (rho0_in * map.inlet_area * U_m1)
end

"""Physical mdot from map flow-coordinate (`phi_in`)."""
function _physical_mdot_from_map_flow_coordinate(
    map::NondimensionalPerformanceMap,
    omega::Real,
    map_flow::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    phi_in = map_flow
    rho0_in = Pt_in / (map.gas_constant * Tt_in)
    U_m1 = omega * map.mean_radius_inlet
    return phi_in * rho0_in * map.inlet_area * U_m1
end

"""Low-level map evaluation in map coordinates."""
function _compressor_performance_map(
    map::NondimensionalPerformanceMap,
    speed_coord::Real,
    flow_coord::Real,
)
    m_tip = speed_coord
    phi_in = flow_coord
    PR = table_evaluate(map.pr_map, m_tip, phi_in)
    eta = table_evaluate(map.eta_map, m_tip, phi_in)
    phi_s = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_surge, m_tip)
    phi_c = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_choke, m_tip)
    return (
        PR=PR,
        eta=eta,
        stall=(phi_in < phi_s),
        choke=(phi_in > phi_c),
        valid=(phi_s <= phi_in <= phi_c),
    )
end

"""
Evaluate a compressor map from physical values and local stagnation state.

Returns `(PR, eta, speed_coord, flow_coord, stall, choke, valid)`.
"""
function performance_from_stagnation(
    map::NondimensionalPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    speed_coord = _map_speed_coordinate_from_stagnation(map, omega, Tt_in, Pt_in)
    flow_coord = _map_flow_coordinate_from_stagnation(map, omega, mdot, Tt_in, Pt_in)
    vals = _compressor_performance_map(map, speed_coord, flow_coord)
    return (
        PR=vals.PR,
        eta=vals.eta,
        speed_coord=speed_coord,
        flow_coord=flow_coord,
        stall=vals.stall,
        choke=vals.choke,
        valid=vals.valid,
    )
end

"""
Physical operating domain for a non-dimensional map at a given inlet state.
"""
function performance_map_domain(
    map::NondimensionalPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
)
    Tt_in > 0 || error("Tt_in must be > 0")
    Pt_in > 0 || error("Pt_in must be > 0")

    a0_in = sqrt(map.gamma * map.gas_constant * Tt_in)
    omega_lo = first(_m_tip_grid(map)) * a0_in / map.tip_radius_inlet
    omega_hi = last(_m_tip_grid(map)) * a0_in / map.tip_radius_inlet

    mdot_vals = Float64[]
    for m_tip in _m_tip_grid(map)
        omega = m_tip * a0_in / map.tip_radius_inlet
        phi_s = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_surge, m_tip)
        phi_c = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_choke, m_tip)
        push!(mdot_vals, _physical_mdot_from_map_flow_coordinate(map, omega, phi_s, Tt_in, Pt_in))
        push!(mdot_vals, _physical_mdot_from_map_flow_coordinate(map, omega, phi_c, Tt_in, Pt_in))
    end

    return (
        omega=(omega_lo, omega_hi),
        mdot=(minimum(mdot_vals), maximum(mdot_vals)),
        mdot_flow_range=(
            surge=(omega -> begin
                m_tip = _map_speed_coordinate_from_stagnation(map, omega, Tt_in, Pt_in)
                phi_s = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_surge, m_tip)
                _physical_mdot_from_map_flow_coordinate(map, omega, phi_s, Tt_in, Pt_in)
            end),
            choke=(omega -> begin
                m_tip = _map_speed_coordinate_from_stagnation(map, omega, Tt_in, Pt_in)
                phi_c = _interp_linear_1d_clamped(_m_tip_grid(map), map.phi_choke, m_tip)
                _physical_mdot_from_map_flow_coordinate(map, omega, phi_c, Tt_in, Pt_in)
            end),
        ),
    )
end

function write_toml(
    map::NondimensionalPerformanceMap,
    path::AbstractString;
    group::AbstractString="compressor_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_nondimensional_performance_map"
    node["format_version"] = 1
    node["gamma"] = map.gamma
    node["gas_constant"] = map.gas_constant
    node["tip_radius_inlet"] = map.tip_radius_inlet
    node["mean_radius_inlet"] = map.mean_radius_inlet
    node["inlet_area"] = map.inlet_area
    node["phi_surge"] = Float64.(map.phi_surge)
    node["phi_choke"] = Float64.(map.phi_choke)
    node["pr_map"] = _table_map_to_toml_dict(map.pr_map)
    node["eta_map"] = _table_map_to_toml_dict(map.eta_map)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{NondimensionalPerformanceMap},
    path::AbstractString;
    group::AbstractString="compressor_map",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    haskey(node, "gamma") || error("missing TOML key gamma")
    haskey(node, "gas_constant") || error("missing TOML key gas_constant")
    haskey(node, "tip_radius_inlet") || error("missing TOML key tip_radius_inlet")
    haskey(node, "mean_radius_inlet") || error("missing TOML key mean_radius_inlet")
    haskey(node, "inlet_area") || error("missing TOML key inlet_area")
    haskey(node, "phi_surge") || error("missing TOML key phi_surge")
    haskey(node, "phi_choke") || error("missing TOML key phi_choke")
    haskey(node, "pr_map") || error("missing TOML key pr_map")
    haskey(node, "eta_map") || error("missing TOML key eta_map")

    gamma = Float64(node["gamma"])
    gas_constant = Float64(node["gas_constant"])
    tip_radius_inlet = Float64(node["tip_radius_inlet"])
    mean_radius_inlet = Float64(node["mean_radius_inlet"])
    inlet_area = Float64(node["inlet_area"])
    phi_surge = Float64.(node["phi_surge"])
    phi_choke = Float64.(node["phi_choke"])
    pr_map = _table_map_from_toml_dict(node["pr_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])

    return NondimensionalPerformanceMap(
        gamma,
        gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        pr_map,
        eta_map,
        phi_surge,
        phi_choke,
    )
end

"""Demo non-dimensional tabulated compressor map for development/testing."""
function demo_nondimensional_tabulated_compressor_performance_map(; interpolation::Symbol=:bilinear)
    NondimensionalPerformanceMap(
        1.4,
        287.05,
        0.22,
        0.18,
        0.060,
        [0.50, 0.70, 0.90],
        [0.30, 0.45, 0.60, 0.75],
        [
            1.20 1.34 1.42 1.45;
            1.35 1.55 1.68 1.73;
            1.50 1.78 1.98 2.05;
        ],
        [
            0.67 0.74 0.72 0.64;
            0.71 0.81 0.80 0.71;
            0.72 0.85 0.84 0.75;
        ];
        interpolation=interpolation,
        phi_surge=[0.33, 0.35, 0.38],
        phi_choke=[0.72, 0.76, 0.80],
    )
end
