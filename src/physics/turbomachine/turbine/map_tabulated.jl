"""
Tabulated turbine performance map implementation (dimensional/corrected-flow form).
"""

using TOML
using ....Utility: AbstractTableMap, interpolation_map, table_evaluate
using ....Utility: linear_evaluate
using ....Utility: table_xgrid, table_ygrid, table_values, table_interpolation
import ....Utility: write_toml, read_toml

"""
Tabulated turbine performance map on corrected coordinates.

Fields:
- `Tt_ref`: reference total temperature for corrected normalization.
- `Pt_ref`: reference total pressure for corrected normalization.
- `mdot_corr_map`: interpolant for corrected mass flow.
- `eta_map`: interpolant for adiabatic efficiency.
- `pr_turb_min`: per-speed lower turbine-pressure-ratio bound.
- `pr_turb_max`: per-speed upper turbine-pressure-ratio bound.
"""
struct TabulatedTurbinePerformanceMap{M<:AbstractTableMap} <: AbstractTurbinePerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    mdot_corr_map::M
    eta_map::M
    pr_turb_min::Vector{Float64}
    pr_turb_max::Vector{Float64}
end

function TabulatedTurbinePerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    mdot_corr_map::M,
    eta_map::M,
    pr_turb_min::Vector{<:Real},
    pr_turb_max::Vector{<:Real},
) where {M<:AbstractTableMap}
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    table_xgrid(mdot_corr_map) == table_xgrid(eta_map) || error("mdot_corr_map/eta_map x grids must match")
    table_ygrid(mdot_corr_map) == table_ygrid(eta_map) || error("mdot_corr_map/eta_map y grids must match")
    xgrid = table_xgrid(mdot_corr_map)
    length(pr_turb_min) == length(xgrid) || error("pr_turb_min length must match omega grid length")
    length(pr_turb_max) == length(xgrid) || error("pr_turb_max length must match omega grid length")
    pr_min = Float64.(pr_turb_min)
    pr_max = Float64.(pr_turb_max)
    all(pr_min .<= pr_max) || error("pr_turb_min must be <= pr_turb_max at every omega grid point")
    return TabulatedTurbinePerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        mdot_corr_map,
        eta_map,
        pr_min,
        pr_max,
    )
end

function TabulatedTurbinePerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    mdot_corr_map::M,
    eta_map::M,
) where {M<:AbstractTableMap}
    ygrid = table_ygrid(mdot_corr_map)
    xgrid = table_xgrid(mdot_corr_map)
    pr_turb_min = fill(Float64(first(ygrid)), length(xgrid))
    pr_turb_max = fill(Float64(last(ygrid)), length(xgrid))
    return TabulatedTurbinePerformanceMap(
        Tt_ref,
        Pt_ref,
        mdot_corr_map,
        eta_map,
        pr_turb_min,
        pr_turb_max,
    )
end

function TabulatedTurbinePerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    omega_corr_grid::Vector{<:Real},
    pr_turb_grid::Vector{<:Real},
    mdot_corr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
    ;
    interpolation::Symbol,
    pr_turb_min::Union{Nothing,Vector{<:Real}}=nothing,
    pr_turb_max::Union{Nothing,Vector{<:Real}}=nothing,
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(omega_corr_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    length(pr_turb_grid) >= 2 || error("pr_turb_grid must have at least 2 points")
    issorted(omega_corr_grid) || error("omega_corr_grid must be sorted ascending")
    issorted(pr_turb_grid) || error("pr_turb_grid must be sorted ascending")
    size(mdot_corr_table) == (length(omega_corr_grid), length(pr_turb_grid)) ||
        error("mdot_corr_table size must match (length(omega_corr_grid), length(pr_turb_grid))")
    size(eta_table) == (length(omega_corr_grid), length(pr_turb_grid)) ||
        error("eta_table size must match (length(omega_corr_grid), length(pr_turb_grid))")

    omega_corr_grid_f = Float64.(omega_corr_grid)
    pr_turb_grid_f = Float64.(pr_turb_grid)
    mdot_corr_table_f = Float64.(mdot_corr_table)
    eta_table_f = Float64.(eta_table)
    mdot_corr_map = interpolation_map(interpolation, omega_corr_grid_f, pr_turb_grid_f, mdot_corr_table_f)
    eta_map = interpolation_map(interpolation, omega_corr_grid_f, pr_turb_grid_f, eta_table_f)

    pr_min = isnothing(pr_turb_min) ?
        fill(first(pr_turb_grid_f), length(omega_corr_grid_f)) :
        Float64.(pr_turb_min)
    pr_max = isnothing(pr_turb_max) ?
        fill(last(pr_turb_grid_f), length(omega_corr_grid_f)) :
        Float64.(pr_turb_max)

    return TabulatedTurbinePerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        mdot_corr_map,
        eta_map,
        pr_min,
        pr_max,
    )
end

_omega_corr_grid(map::TabulatedTurbinePerformanceMap) = table_xgrid(map.mdot_corr_map)
_pr_turb_grid(map::TabulatedTurbinePerformanceMap) = table_ygrid(map.mdot_corr_map)
_mdot_corr_table(map::TabulatedTurbinePerformanceMap) = table_values(map.mdot_corr_map)
_eta_table(map::TabulatedTurbinePerformanceMap) = table_values(map.eta_map)
_interpolation_kind(map::TabulatedTurbinePerformanceMap) = table_interpolation(map.mdot_corr_map)

omega_corr_grid(map::TabulatedTurbinePerformanceMap) = _omega_corr_grid(map)
pr_turb_grid(map::TabulatedTurbinePerformanceMap) = _pr_turb_grid(map)
mdot_corr_table(map::TabulatedTurbinePerformanceMap) = _mdot_corr_table(map)
eta_table(map::TabulatedTurbinePerformanceMap) = _eta_table(map)
interpolation_kind(map::TabulatedTurbinePerformanceMap) = _interpolation_kind(map)

"""
Corrected shaft speed from physical speed and local total temperature.
"""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""
Corrected mass flow from physical flow and local total conditions.
"""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

_map_speed_coordinate_from_stagnation(
    map::TabulatedTurbinePerformanceMap,
    omega::Real,
    Tt_in::Real,
    Pt_in::Real,
) = corrected_speed(omega, Tt_in, map)

_map_pr_coordinate_from_stagnation(
    map::TabulatedTurbinePerformanceMap,
    Pt_in::Real,
    Pt_out::Real,
) = Pt_in / Pt_out

_physical_flow_from_map_flow_coordinate(
    map::TabulatedTurbinePerformanceMap,
    flow_coord::Real,
    Tt_in::Real,
    Pt_in::Real,
) = flow_coord * (Pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)

"""Low-level map evaluation in map coordinates."""
function _turbine_performance_map(
    map::TabulatedTurbinePerformanceMap,
    speed_coord::Real,
    pr_coord::Real,
)
    mdot_corr = table_evaluate(map.mdot_corr_map, speed_coord, pr_coord)
    eta = table_evaluate(map.eta_map, speed_coord, pr_coord)
    pr_min = linear_evaluate(_omega_corr_grid(map), map.pr_turb_min, speed_coord)
    pr_max = linear_evaluate(_omega_corr_grid(map), map.pr_turb_max, speed_coord)
    return (
        mdot_corr=mdot_corr,
        eta=eta,
        low_pr=(pr_coord < pr_min),
        high_pr=(pr_coord > pr_max),
        valid=(pr_min <= pr_coord <= pr_max),
    )
end

"""
Evaluate a turbine map at corrected coordinates.

Returns `(mdot_corr, eta, low_pr, high_pr, valid)`.
"""
function turbine_performance_map(
    map::TabulatedTurbinePerformanceMap,
    omega_corr::Real,
    pr_turb::Real,
)
    return _turbine_performance_map(map, omega_corr, pr_turb)
end

"""
Evaluate a turbine map from physical values and local stagnation state.

Returns `(omega_corr, PR_turb, mdot_corr, mdot, eta, low_pr, high_pr, valid)`.
"""
function turbine_performance_map_from_stagnation(
    map::TabulatedTurbinePerformanceMap,
    omega::Real,
    Pt_in::Real,
    Pt_out::Real,
    Tt_in::Real,
)
    Pt_out > 0 || error("Pt_out must be > 0")
    omega_corr = _map_speed_coordinate_from_stagnation(map, omega, Tt_in, Pt_in)
    PR_turb = _map_pr_coordinate_from_stagnation(map, Pt_in, Pt_out)
    vals = _turbine_performance_map(map, omega_corr, PR_turb)
    mdot = _physical_flow_from_map_flow_coordinate(map, vals.mdot_corr, Tt_in, Pt_in)
    return (
        omega_corr=omega_corr,
        PR_turb=PR_turb,
        mdot_corr=vals.mdot_corr,
        mdot=mdot,
        eta=vals.eta,
        low_pr=vals.low_pr,
        high_pr=vals.high_pr,
        valid=vals.valid,
    )
end

"""
Recommended corrected-coordinate operating domain for a tabulated turbine map.
"""
function performance_map_domain(map::TabulatedTurbinePerformanceMap)
    return (
        omega_corr=(first(_omega_corr_grid(map)), last(_omega_corr_grid(map))),
        pr_turb=(minimum(map.pr_turb_min), maximum(map.pr_turb_max)),
        pr_turb_range=(
            min=(omega -> linear_evaluate(_omega_corr_grid(map), map.pr_turb_min, omega)),
            max=(omega -> linear_evaluate(_omega_corr_grid(map), map.pr_turb_max, omega)),
        ),
    )
end

function _table_to_rows(table::AbstractMatrix{<:Real})
    return [Float64.(collect(view(table, i, :))) for i in 1:size(table, 1)]
end

function _rows_to_table(rows::Vector)
    length(rows) >= 1 || error("table must have at least one row")
    ncols = length(rows[1])
    ncols >= 1 || error("table rows must have at least one column")
    all(length(row) == ncols for row in rows) || error("table rows must have consistent lengths")
    table = Matrix{Float64}(undef, length(rows), ncols)
    for i in eachindex(rows)
        table[i, :] = Float64.(rows[i])
    end
    return table
end

function _table_map_to_toml_dict(map::AbstractTableMap)
    return Dict{String,Any}(
        "interpolation" => String(table_interpolation(map)),
        "xgrid" => Float64.(table_xgrid(map)),
        "ygrid" => Float64.(table_ygrid(map)),
        "table" => _table_to_rows(table_values(map)),
    )
end

function _table_map_from_toml_dict(data::Dict{String,Any})
    haskey(data, "interpolation") || error("missing TOML key interpolation")
    haskey(data, "xgrid") || error("missing TOML key xgrid")
    haskey(data, "ygrid") || error("missing TOML key ygrid")
    haskey(data, "table") || error("missing TOML key table")
    interpolation = Symbol(String(data["interpolation"]))
    xgrid = Float64.(data["xgrid"])
    ygrid = Float64.(data["ygrid"])
    table = _rows_to_table(data["table"])
    return interpolation_map(interpolation, xgrid, ygrid, table)
end

function _find_or_create_group!(data::Dict{String,Any}, group::AbstractString)
    isempty(group) && return data
    node = data
    for key in split(group, '.')
        if !haskey(node, key)
            node[key] = Dict{String,Any}()
        end
        child = node[key]
        child isa Dict || error("group path conflicts with non-table key $(key)")
        node = child
    end
    return node
end

function _find_group(data::Dict{String,Any}, group::AbstractString)
    isempty(group) && return data
    node = data
    for key in split(group, '.')
        haskey(node, key) || error("missing TOML group $(group)")
        child = node[key]
        child isa Dict || error("TOML group $(group) is not a table")
        node = child
    end
    return node
end

function write_toml(
    map::TabulatedTurbinePerformanceMap,
    path::AbstractString;
    group::AbstractString="turbine_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "turbine_performance_map"
    node["format_version"] = 2
    node["Tt_ref"] = map.Tt_ref
    node["Pt_ref"] = map.Pt_ref
    node["pr_turb_min"] = Float64.(map.pr_turb_min)
    node["pr_turb_max"] = Float64.(map.pr_turb_max)
    node["mdot_corr_map"] = _table_map_to_toml_dict(map.mdot_corr_map)
    node["eta_map"] = _table_map_to_toml_dict(map.eta_map)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{TabulatedTurbinePerformanceMap},
    path::AbstractString;
    group::AbstractString="turbine_map",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    haskey(node, "Tt_ref") || error("missing TOML key Tt_ref")
    haskey(node, "Pt_ref") || error("missing TOML key Pt_ref")
    haskey(node, "mdot_corr_map") || error("missing TOML key mdot_corr_map")
    haskey(node, "eta_map") || error("missing TOML key eta_map")
    Tt_ref = Float64(node["Tt_ref"])
    Pt_ref = Float64(node["Pt_ref"])
    pr_turb_min = haskey(node, "pr_turb_min") ? Float64.(node["pr_turb_min"]) : nothing
    pr_turb_max = haskey(node, "pr_turb_max") ? Float64.(node["pr_turb_max"]) : nothing
    mdot_corr_map = _table_map_from_toml_dict(node["mdot_corr_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])
    if isnothing(pr_turb_min) || isnothing(pr_turb_max)
        return TabulatedTurbinePerformanceMap(Tt_ref, Pt_ref, mdot_corr_map, eta_map)
    end
    return TabulatedTurbinePerformanceMap(
        Tt_ref,
        Pt_ref,
        mdot_corr_map,
        eta_map,
        pr_turb_min,
        pr_turb_max,
    )
end

"""Demo tabulated turbine map for development/testing."""
function demo_tabulated_turbine_performance_map(; interpolation::Symbol=:bilinear)
    TabulatedTurbinePerformanceMap(
        288.15,
        101_325.0,
        [0.6, 0.8, 1.0],
        [1.4, 1.8, 2.2],
        [
            10.0 12.0 13.5;
            12.5 14.5 16.0;
            14.0 16.5 18.5;
        ],
        [
            0.80 0.84 0.82;
            0.83 0.87 0.85;
            0.81 0.86 0.84;
        ];
        interpolation=interpolation,
        pr_turb_min=[1.4, 1.4, 1.4],
        pr_turb_max=[2.2, 2.2, 2.2],
    )
end
