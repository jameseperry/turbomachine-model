"""
Tabulated compressor performance map implementation (dimensional/corrected-flow form).
"""

using TOML
using ....Utility: AbstractTableMap, interpolation_map, table_evaluate
using ....Utility: linear_evaluate
using ....Utility: table_xgrid, table_ygrid, table_values, table_interpolation
import ....Utility: write_toml, read_toml

"""
Tabulated compressor performance map on physical speed and corrected flow.

Fields:
- `Tt_ref`: reference total temperature for flow correction.
- `Pt_ref`: reference total pressure for corrected normalization.
- `pr_map`: interpolant for total-pressure-ratio (`Pt_out/Pt_in`).
- `eta_map`: interpolant for adiabatic efficiency.
"""
struct TabulatedCompressorPerformanceMap{M<:AbstractTableMap} <: AbstractCompressorPerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    pr_map::M
    eta_map::M
    mdot_corr_surge::Vector{Float64}
    mdot_corr_choke::Vector{Float64}
end

function TabulatedCompressorPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    pr_map::M,
    eta_map::M,
    mdot_corr_surge::Vector{<:Real},
    mdot_corr_choke::Vector{<:Real},
) where {M<:AbstractTableMap}
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    table_xgrid(pr_map) == table_xgrid(eta_map) || error("pr_map/eta_map x grids must match")
    table_ygrid(pr_map) == table_ygrid(eta_map) || error("pr_map/eta_map y grids must match")
    xgrid = table_xgrid(pr_map)
    length(mdot_corr_surge) == length(xgrid) || error("mdot_corr_surge length must match omega grid length")
    length(mdot_corr_choke) == length(xgrid) || error("mdot_corr_choke length must match omega grid length")
    surge = Float64.(mdot_corr_surge)
    choke = Float64.(mdot_corr_choke)
    all(surge .<= choke) || error("mdot_corr_surge must be <= mdot_corr_choke at every omega grid point")
    return TabulatedCompressorPerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        pr_map,
        eta_map,
        surge,
        choke,
    )
end

function TabulatedCompressorPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    pr_map::M,
    eta_map::M,
) where {M<:AbstractTableMap}
    ygrid = table_ygrid(pr_map)
    xgrid = table_xgrid(pr_map)
    mdot_corr_surge = fill(Float64(first(ygrid)), length(xgrid))
    mdot_corr_choke = fill(Float64(last(ygrid)), length(xgrid))
    return TabulatedCompressorPerformanceMap(
        Tt_ref,
        Pt_ref,
        pr_map,
        eta_map,
        mdot_corr_surge,
        mdot_corr_choke,
    )
end

function TabulatedCompressorPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    omega_corr_grid::Vector{<:Real},
    mdot_corr_grid::Vector{<:Real},
    pr_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
    ;
    interpolation::Symbol,
    mdot_corr_surge::Union{Nothing,Vector{<:Real}}=nothing,
    mdot_corr_choke::Union{Nothing,Vector{<:Real}}=nothing,
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(omega_corr_grid) >= 2 || error("omega_corr_grid must have at least 2 points")
    length(mdot_corr_grid) >= 2 || error("mdot_corr_grid must have at least 2 points")
    issorted(omega_corr_grid) || error("omega_corr_grid must be sorted ascending")
    issorted(mdot_corr_grid) || error("mdot_corr_grid must be sorted ascending")
    size(pr_table) == (length(omega_corr_grid), length(mdot_corr_grid)) ||
        error("pr_table size must match (length(omega_corr_grid), length(mdot_corr_grid))")
    size(eta_table) == (length(omega_corr_grid), length(mdot_corr_grid)) ||
        error("eta_table size must match (length(omega_corr_grid), length(mdot_corr_grid))")

    omega_corr_grid_f = Float64.(omega_corr_grid)
    mdot_corr_grid_f = Float64.(mdot_corr_grid)
    pr_table_f = Float64.(pr_table)
    eta_table_f = Float64.(eta_table)
    pr_map = interpolation_map(interpolation, omega_corr_grid_f, mdot_corr_grid_f, pr_table_f)
    eta_map = interpolation_map(interpolation, omega_corr_grid_f, mdot_corr_grid_f, eta_table_f)

    surge = isnothing(mdot_corr_surge) ?
        fill(first(mdot_corr_grid_f), length(omega_corr_grid_f)) :
        Float64.(mdot_corr_surge)
    choke = isnothing(mdot_corr_choke) ?
        fill(last(mdot_corr_grid_f), length(omega_corr_grid_f)) :
        Float64.(mdot_corr_choke)

    return TabulatedCompressorPerformanceMap(Float64(Tt_ref), Float64(Pt_ref), pr_map, eta_map, surge, choke)
end

_omega_corr_grid(map::TabulatedCompressorPerformanceMap) = table_xgrid(map.pr_map)
_mdot_corr_grid(map::TabulatedCompressorPerformanceMap) = table_ygrid(map.pr_map)
_pr_table(map::TabulatedCompressorPerformanceMap) = table_values(map.pr_map)
_eta_table(map::TabulatedCompressorPerformanceMap) = table_values(map.eta_map)
_interpolation_kind(map::TabulatedCompressorPerformanceMap) = table_interpolation(map.pr_map)

"""Map-coordinate speed for tabulated corrected-flow maps (physical omega)."""
_map_speed_coordinate_from_stagnation(
    map::TabulatedCompressorPerformanceMap,
    omega::Real,
    Tt_in::Real,
    Pt_in::Real,
) = omega

"""Map-coordinate flow for tabulated corrected-flow maps (mdot_corr)."""
_map_flow_coordinate_from_stagnation(
    map::TabulatedCompressorPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
) = mdot * sqrt(Tt_in / map.Tt_ref) / (Pt_in / map.Pt_ref)

"""Physical mdot from map flow-coordinate for tabulated corrected-flow maps."""
function _physical_mdot_from_map_flow_coordinate(
    map::TabulatedCompressorPerformanceMap,
    omega::Real,
    map_flow::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    mdot_corr = map_flow
    return mdot_corr * (Pt_in / map.Pt_ref) / sqrt(Tt_in / map.Tt_ref)
end

"""Low-level map evaluation in map coordinates."""
function _compressor_performance_map(
    map::TabulatedCompressorPerformanceMap,
    speed_coord::Real,
    flow_coord::Real,
)
    PR = table_evaluate(map.pr_map, speed_coord, flow_coord)
    eta = table_evaluate(map.eta_map, speed_coord, flow_coord)
    mdot_s = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_surge, speed_coord)
    mdot_c = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_choke, speed_coord)
    return (
        PR=PR,
        eta=eta,
        stall=(flow_coord < mdot_s),
        choke=(flow_coord > mdot_c),
        valid=(mdot_s <= flow_coord <= mdot_c),
    )
end

"""
Evaluate a compressor map from physical values and local stagnation state.

Returns `(PR, eta, speed_coord, flow_coord, stall, choke, valid)`.
"""
function performance_from_stagnation(
    map::TabulatedCompressorPerformanceMap,
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
Physical operating domain for a tabulated corrected-flow map at a given inlet state.
"""
function performance_map_domain(
    map::TabulatedCompressorPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
)
    speed_lo = first(_omega_corr_grid(map))
    speed_hi = last(_omega_corr_grid(map))
    mdot_vals = Float64[]
    for omega in _omega_corr_grid(map)
        mdot_s = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_surge, omega)
        mdot_c = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_choke, omega)
        push!(mdot_vals, _physical_mdot_from_map_flow_coordinate(map, omega, mdot_s, Tt_in, Pt_in))
        push!(mdot_vals, _physical_mdot_from_map_flow_coordinate(map, omega, mdot_c, Tt_in, Pt_in))
    end
    return (
        omega=(speed_lo, speed_hi),
        mdot=(minimum(mdot_vals), maximum(mdot_vals)),
        mdot_flow_range=(
            surge=(omega -> begin
                mdot_s = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_surge, omega)
                _physical_mdot_from_map_flow_coordinate(map, omega, mdot_s, Tt_in, Pt_in)
            end),
            choke=(omega -> begin
                mdot_c = linear_evaluate(_omega_corr_grid(map), map.mdot_corr_choke, omega)
                _physical_mdot_from_map_flow_coordinate(map, omega, mdot_c, Tt_in, Pt_in)
            end),
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
    map::TabulatedCompressorPerformanceMap,
    path::AbstractString;
    group::AbstractString="compressor_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "compressor_performance_map"
    node["format_version"] = 2
    node["Tt_ref"] = map.Tt_ref
    node["Pt_ref"] = map.Pt_ref
    node["mdot_corr_surge"] = Float64.(map.mdot_corr_surge)
    node["mdot_corr_choke"] = Float64.(map.mdot_corr_choke)
    node["pr_map"] = _table_map_to_toml_dict(map.pr_map)
    node["eta_map"] = _table_map_to_toml_dict(map.eta_map)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{TabulatedCompressorPerformanceMap},
    path::AbstractString;
    group::AbstractString="compressor_map",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    haskey(node, "Tt_ref") || error("missing TOML key Tt_ref")
    haskey(node, "Pt_ref") || error("missing TOML key Pt_ref")
    haskey(node, "pr_map") || error("missing TOML key pr_map")
    haskey(node, "eta_map") || error("missing TOML key eta_map")
    Tt_ref = Float64(node["Tt_ref"])
    Pt_ref = Float64(node["Pt_ref"])
    mdot_corr_surge = haskey(node, "mdot_corr_surge") ? Float64.(node["mdot_corr_surge"]) : nothing
    mdot_corr_choke = haskey(node, "mdot_corr_choke") ? Float64.(node["mdot_corr_choke"]) : nothing
    pr_map = _table_map_from_toml_dict(node["pr_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])
    if isnothing(mdot_corr_surge) || isnothing(mdot_corr_choke)
        return TabulatedCompressorPerformanceMap(Tt_ref, Pt_ref, pr_map, eta_map)
    end
    return TabulatedCompressorPerformanceMap(
        Tt_ref,
        Pt_ref,
        pr_map,
        eta_map,
        mdot_corr_surge,
        mdot_corr_choke,
    )
end

"""Demo tabulated compressor map for development/testing."""
function demo_tabulated_compressor_performance_map(; interpolation::Symbol=:bilinear)
    TabulatedCompressorPerformanceMap(
        288.15,
        101_325.0,
        [600.0, 800.0, 1000.0],
        [12.0, 16.0, 20.0],
        [
            1.35 1.55 1.70;
            1.55 1.80 2.00;
            1.70 2.00 2.25;
        ],
        [
            0.74 0.78 0.75;
            0.78 0.83 0.80;
            0.76 0.82 0.79;
        ];
        interpolation=interpolation,
        mdot_corr_surge=[12.0, 12.0, 12.0],
        mdot_corr_choke=[20.0, 20.0, 20.0],
    )
end
