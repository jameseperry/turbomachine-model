"""
Tabulated turbine performance map implementation.
"""

using HDF5
using TOML
using ....Utility: AbstractTableMap, interpolation_map, table_evaluate
using ....Utility: table_xgrid, table_ygrid, table_values, table_interpolation
import ....Utility: write_hdf5, read_hdf5, write_toml, read_toml

abstract type AbstractTurbinePerformanceMap end

"""
Tabulated turbine performance map on corrected coordinates.

Fields:
- `Tt_ref`: reference total temperature for corrected normalization.
- `Pt_ref`: reference total pressure for corrected normalization.
- `mdot_corr_map`: interpolant for corrected mass flow.
- `eta_map`: interpolant for adiabatic efficiency.
"""
struct TabulatedTurbinePerformanceMap{M<:AbstractTableMap} <: AbstractTurbinePerformanceMap
    Tt_ref::Float64
    Pt_ref::Float64
    mdot_corr_map::M
    eta_map::M
end

function TabulatedTurbinePerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    mdot_corr_map::M,
    eta_map::M,
) where {M<:AbstractTableMap}
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    table_xgrid(mdot_corr_map) == table_xgrid(eta_map) || error("mdot_corr_map/eta_map x grids must match")
    table_ygrid(mdot_corr_map) == table_ygrid(eta_map) || error("mdot_corr_map/eta_map y grids must match")
    return TabulatedTurbinePerformanceMap(Float64(Tt_ref), Float64(Pt_ref), mdot_corr_map, eta_map)
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

    return TabulatedTurbinePerformanceMap(Float64(Tt_ref), Float64(Pt_ref), mdot_corr_map, eta_map)
end

omega_corr_grid(map::TabulatedTurbinePerformanceMap) = table_xgrid(map.mdot_corr_map)
pr_turb_grid(map::TabulatedTurbinePerformanceMap) = table_ygrid(map.mdot_corr_map)
mdot_corr_table(map::TabulatedTurbinePerformanceMap) = table_values(map.mdot_corr_map)
eta_table(map::TabulatedTurbinePerformanceMap) = table_values(map.eta_map)
interpolation_kind(map::TabulatedTurbinePerformanceMap) = table_interpolation(map.mdot_corr_map)

"""Corrected shaft speed from physical speed and local total temperature."""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""Corrected mass flow from physical flow and local total conditions."""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

corrected_speed(omega::Real, Tt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_speed(omega, Tt_in, map.Tt_ref)

corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, map::AbstractTurbinePerformanceMap) =
    corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)

_physical_flow_from_corrected(
    mdot_corr::Real,
    Tt_in::Real,
    Pt_in::Real,
    Tt_ref::Real,
    Pt_ref::Real,
) = mdot_corr * (Pt_in / Pt_ref) / sqrt(Tt_in / Tt_ref)

"""
Evaluate a turbine map at corrected coordinates.

Returns named tuple `(mdot_corr, eta)` where:
- `mdot_corr` is corrected mass flow
- `eta` is adiabatic efficiency
"""
function turbine_performance_map(
    map::TabulatedTurbinePerformanceMap,
    omega_corr::Real,
    pr_turb::Real,
)
    mdot_corr = table_evaluate(map.mdot_corr_map, omega_corr, pr_turb)
    eta = table_evaluate(map.eta_map, omega_corr, pr_turb)
    return (mdot_corr=mdot_corr, eta=eta)
end

"""
Evaluate a turbine map from physical values and local stagnation state.

Returns `(omega_corr, PR_turb, mdot_corr, mdot, eta)` where
`PR_turb = Pt_in / Pt_out`.
"""
function turbine_performance_map_from_stagnation(
    map::AbstractTurbinePerformanceMap,
    omega::Real,
    Pt_in::Real,
    Pt_out::Real,
    Tt_in::Real,
)
    Pt_out > 0 || error("Pt_out must be > 0")
    PR_turb = Pt_in / Pt_out
    omega_corr = corrected_speed(omega, Tt_in, map)
    vals = turbine_performance_map(map, omega_corr, PR_turb)
    mdot = _physical_flow_from_corrected(vals.mdot_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
    return (
        omega_corr=omega_corr,
        PR_turb=PR_turb,
        mdot_corr=vals.mdot_corr,
        mdot=mdot,
        eta=vals.eta,
    )
end

function write_hdf5(
    parent::Union{HDF5.File, HDF5.Group},
    name::AbstractString,
    map::TabulatedTurbinePerformanceMap,
)
    haskey(parent, name) && error("HDF5 object $(name) already exists")
    group = create_group(parent, name)
    attrs(group)["Tt_ref"] = map.Tt_ref
    attrs(group)["Pt_ref"] = map.Pt_ref
    write_hdf5(group, "mdot_corr_map", map.mdot_corr_map)
    write_hdf5(group, "eta_map", map.eta_map)
    return nothing
end

function read_hdf5(
    ::Type{TabulatedTurbinePerformanceMap},
    parent::Union{HDF5.File, HDF5.Group},
    name::AbstractString,
)
    haskey(parent, name) || error("missing HDF5 object $(name)")
    group = parent[name]
    haskey(attrs(group), "Tt_ref") || error("missing attribute Tt_ref")
    haskey(attrs(group), "Pt_ref") || error("missing attribute Pt_ref")
    Tt_ref = Float64(attrs(group)["Tt_ref"])
    Pt_ref = Float64(attrs(group)["Pt_ref"])
    mdot_corr_map = read_hdf5(AbstractTableMap, group, "mdot_corr_map")
    eta_map = read_hdf5(AbstractTableMap, group, "eta_map")
    return TabulatedTurbinePerformanceMap(Tt_ref, Pt_ref, mdot_corr_map, eta_map)
end

function write_hdf5(
    map::TabulatedTurbinePerformanceMap,
    path::AbstractString;
    group::AbstractString="turbine_map",
)
    h5open(path, "w") do file
        write_hdf5(file, group, map)
    end
    return path
end

function read_hdf5(
    ::Type{TabulatedTurbinePerformanceMap},
    path::AbstractString;
    group::AbstractString="turbine_map",
)
    return h5open(path, "r") do file
        read_hdf5(TabulatedTurbinePerformanceMap, file, group)
    end
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
    node["format_version"] = 1
    node["Tt_ref"] = map.Tt_ref
    node["Pt_ref"] = map.Pt_ref
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
    mdot_corr_map = _table_map_from_toml_dict(node["mdot_corr_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])
    return TabulatedTurbinePerformanceMap(Tt_ref, Pt_ref, mdot_corr_map, eta_map)
end

"""Demo tabulated turbine map for development/testing."""
function demo_turbine_performance_map(; interpolation::Symbol=:bilinear)
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
    )
end
