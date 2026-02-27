"""
Tabulated compressor performance map implementation.
"""

using TOML
using ....Utility: AbstractTableMap, interpolation_map, table_evaluate
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
end

function TabulatedCompressorPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    pr_map::M,
    eta_map::M,
) where {M<:AbstractTableMap}
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    table_xgrid(pr_map) == table_xgrid(eta_map) || error("pr_map/eta_map x grids must match")
    table_ygrid(pr_map) == table_ygrid(eta_map) || error("pr_map/eta_map y grids must match")
    return TabulatedCompressorPerformanceMap(Float64(Tt_ref), Float64(Pt_ref), pr_map, eta_map)
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

    return TabulatedCompressorPerformanceMap(Float64(Tt_ref), Float64(Pt_ref), pr_map, eta_map)
end

omega_corr_grid(map::TabulatedCompressorPerformanceMap) = table_xgrid(map.pr_map)
mdot_corr_grid(map::TabulatedCompressorPerformanceMap) = table_ygrid(map.pr_map)
pr_table(map::TabulatedCompressorPerformanceMap) = table_values(map.pr_map)
eta_table(map::TabulatedCompressorPerformanceMap) = table_values(map.eta_map)
interpolation_kind(map::TabulatedCompressorPerformanceMap) = table_interpolation(map.pr_map)
function performance_map_domain(map::TabulatedCompressorPerformanceMap)
    mdot_lo = first(mdot_corr_grid(map))
    mdot_hi = last(mdot_corr_grid(map))
    return (
        omega_corr=(first(omega_corr_grid(map)), last(omega_corr_grid(map))),
        mdot_corr=(mdot_lo, mdot_hi),
        mdot_corr_flow_range=(
            surge=(omega_corr -> mdot_lo),
            choke=(omega_corr -> mdot_hi),
        ),
    )
end

"""Tabulated maps currently use physical shaft speed directly (omega in rad/s)."""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) = omega

"""Corrected mass flow from physical flow and local total conditions."""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::TabulatedCompressorPerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractCompressorPerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

"""
Evaluate a compressor map at corrected coordinates.

Returns named tuple `(PR, eta)` where:
- `PR` is total-pressure ratio (`Pt_out/Pt_in`)
- `eta` is adiabatic efficiency
"""
function compressor_performance_map(
    map::TabulatedCompressorPerformanceMap,
    omega_corr::Real,
    mdot_corr::Real,
)
    PR = table_evaluate(map.pr_map, omega_corr, mdot_corr)
    eta = table_evaluate(map.eta_map, omega_corr, mdot_corr)
    return (PR=PR, eta=eta)
end

"""
Evaluate a compressor map from physical values and local stagnation state.

Returns `(omega_corr, mdot_corr, PR, eta)`.
"""
function compressor_performance_map_from_stagnation(
    map::AbstractCompressorPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    omega_corr = corrected_speed(omega, Tt_in, map)
    mdot_corr = corrected_flow(mdot, Tt_in, Pt_in, map)
    vals = compressor_performance_map(map, omega_corr, mdot_corr)
    return (omega_corr=omega_corr, mdot_corr=mdot_corr, PR=vals.PR, eta=vals.eta)
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
    node["format_version"] = 1
    node["Tt_ref"] = map.Tt_ref
    node["Pt_ref"] = map.Pt_ref
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
    pr_map = _table_map_from_toml_dict(node["pr_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])
    return TabulatedCompressorPerformanceMap(Tt_ref, Pt_ref, pr_map, eta_map)
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
    )
end
