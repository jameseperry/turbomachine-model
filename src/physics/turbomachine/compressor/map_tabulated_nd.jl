"""
Tabulated compressor performance map implementation (non-dimensional form).
"""

"""
Tabulated compressor map over non-dimensional coordinates.

Domain inputs:
- `M_tip = U_tip / a0_in`
- `phi_in = Vx_1 / U_m_1`

Outputs:
- `PR = Pt_out / Pt_in`
- `eta` total-to-total efficiency

Geometry fields are used to convert physical `(omega, mdot, Tt_in, Pt_in)` into
`(M_tip, phi_in)`.
"""
struct NonDimensionalTabulatedCompressorPerformanceMap{M<:AbstractTableMap} <: AbstractCompressorPerformanceMap
    gamma::Float64
    gas_constant::Float64
    tip_radius_inlet::Float64
    mean_radius_inlet::Float64
    inlet_area::Float64
    pr_map::M
    eta_map::M
    phi_surge::Vector{Float64}
    phi_choke::Vector{Float64}
end

function NonDimensionalTabulatedCompressorPerformanceMap(
    gamma::Real,
    gas_constant::Real,
    tip_radius_inlet::Real,
    mean_radius_inlet::Real,
    inlet_area::Real,
    pr_map::M,
    eta_map::M,
    phi_surge::Vector{<:Real},
    phi_choke::Vector{<:Real},
) where {M<:AbstractTableMap}
    gamma > 1 || error("gamma must be > 1")
    gas_constant > 0 || error("gas_constant must be > 0")
    tip_radius_inlet > 0 || error("tip_radius_inlet must be > 0")
    mean_radius_inlet > 0 || error("mean_radius_inlet must be > 0")
    inlet_area > 0 || error("inlet_area must be > 0")
    table_xgrid(pr_map) == table_xgrid(eta_map) || error("pr_map/eta_map x grids must match")
    table_ygrid(pr_map) == table_ygrid(eta_map) || error("pr_map/eta_map y grids must match")

    m_grid = table_xgrid(pr_map)
    length(phi_surge) == length(m_grid) || error("phi_surge length must match M_tip grid length")
    length(phi_choke) == length(m_grid) || error("phi_choke length must match M_tip grid length")

    phi_surge_f = Float64.(phi_surge)
    phi_choke_f = Float64.(phi_choke)
    all(phi_surge_f .<= phi_choke_f) || error("phi_surge must be <= phi_choke at every M_tip grid point")

    return NonDimensionalTabulatedCompressorPerformanceMap(
        Float64(gamma),
        Float64(gas_constant),
        Float64(tip_radius_inlet),
        Float64(mean_radius_inlet),
        Float64(inlet_area),
        pr_map,
        eta_map,
        phi_surge_f,
        phi_choke_f,
    )
end

function NonDimensionalTabulatedCompressorPerformanceMap(
    gamma::Real,
    gas_constant::Real,
    tip_radius_inlet::Real,
    mean_radius_inlet::Real,
    inlet_area::Real,
    m_tip_grid::Vector{<:Real},
    phi_in_grid::Vector{<:Real},
    pr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real};
    interpolation::Symbol,
    phi_surge::Union{Nothing,Vector{<:Real}}=nothing,
    phi_choke::Union{Nothing,Vector{<:Real}}=nothing,
)
    length(m_tip_grid) >= 2 || error("m_tip_grid must have at least 2 points")
    length(phi_in_grid) >= 2 || error("phi_in_grid must have at least 2 points")
    issorted(m_tip_grid) || error("m_tip_grid must be sorted ascending")
    issorted(phi_in_grid) || error("phi_in_grid must be sorted ascending")
    size(pr_table) == (length(m_tip_grid), length(phi_in_grid)) ||
        error("pr_table size must match (length(m_tip_grid), length(phi_in_grid))")
    size(eta_table) == (length(m_tip_grid), length(phi_in_grid)) ||
        error("eta_table size must match (length(m_tip_grid), length(phi_in_grid))")

    m_tip_grid_f = Float64.(m_tip_grid)
    phi_in_grid_f = Float64.(phi_in_grid)
    pr_table_f = Float64.(pr_table)
    eta_table_f = Float64.(eta_table)
    pr_map = interpolation_map(interpolation, m_tip_grid_f, phi_in_grid_f, pr_table_f)
    eta_map = interpolation_map(interpolation, m_tip_grid_f, phi_in_grid_f, eta_table_f)

    phi_surge_f = isnothing(phi_surge) ? fill(first(phi_in_grid_f), length(m_tip_grid_f)) : Float64.(phi_surge)
    phi_choke_f = isnothing(phi_choke) ? fill(last(phi_in_grid_f), length(m_tip_grid_f)) : Float64.(phi_choke)

    return NonDimensionalTabulatedCompressorPerformanceMap(
        gamma,
        gas_constant,
        tip_radius_inlet,
        mean_radius_inlet,
        inlet_area,
        pr_map,
        eta_map,
        phi_surge_f,
        phi_choke_f,
    )
end

_m_tip_grid(map::NonDimensionalTabulatedCompressorPerformanceMap) = table_xgrid(map.pr_map)
_phi_in_grid(map::NonDimensionalTabulatedCompressorPerformanceMap) = table_ygrid(map.pr_map)
_pr_table(map::NonDimensionalTabulatedCompressorPerformanceMap) = table_values(map.pr_map)
_eta_table(map::NonDimensionalTabulatedCompressorPerformanceMap) = table_values(map.eta_map)
_interpolation_kind(map::NonDimensionalTabulatedCompressorPerformanceMap) = table_interpolation(map.pr_map)

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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
function compressor_performance_map_from_stagnation(
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    map::NonDimensionalTabulatedCompressorPerformanceMap,
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
    ::Type{NonDimensionalTabulatedCompressorPerformanceMap},
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

    return NonDimensionalTabulatedCompressorPerformanceMap(
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
    NonDimensionalTabulatedCompressorPerformanceMap(
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
