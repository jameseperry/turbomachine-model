using TOML
using ...Utility: AbstractTableMap, interpolation_map, table_evaluate
using ...Utility: linear_evaluate
using ...Utility: table_xgrid, table_ygrid, table_values, table_interpolation
import ...Utility: write_toml, read_toml

"""
Shared tabulated performance map core.

The table storage uses corrected coordinates internally:
- corrected speed on the x-axis
- corrected flow on the y-axis

That corrected-coordinate representation is an implementation detail of
`TabulatedPerformanceMap`. Callers should evaluate the map using physical
`omega` and physical `mdot` via `performance_from_stagnation(...)`.
"""
struct TabulatedPerformanceMap{M<:AbstractTableMap}
    Tt_ref::Float64
    Pt_ref::Float64
    pressure_ratio_map::M
    eta_map::M
    flow_min::Vector{Float64}
    flow_max::Vector{Float64}
end

"""
Corrected shaft speed from physical speed and local total temperature.
"""
corrected_speed(omega::Real, Tt_in::Real, Tt_ref::Real) =
    omega / sqrt(Tt_in / Tt_ref)

"""
Physical shaft speed from corrected speed and local total temperature.
"""
physical_speed(omega_corr::Real, Tt_in::Real, Tt_ref::Real) =
    omega_corr * sqrt(Tt_in / Tt_ref)

"""
Corrected mass flow from physical flow and local total conditions.
"""
corrected_flow(mdot::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot * sqrt(Tt_in / Tt_ref) / (Pt_in / Pt_ref)

"""
Physical mass flow from corrected flow and local total conditions.
"""
physical_flow(mdot_corr::Real, Tt_in::Real, Pt_in::Real, Tt_ref::Real, Pt_ref::Real) =
    mdot_corr * (Pt_in / Pt_ref) / sqrt(Tt_in / Tt_ref)

function TabulatedPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    pressure_ratio_map::M,
    eta_map::M,
    flow_min::Vector{<:Real},
    flow_max::Vector{<:Real},
) where {M<:AbstractTableMap}
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    table_xgrid(pressure_ratio_map) == table_xgrid(eta_map) || error("pressure_ratio_map/eta_map x grids must match")
    table_ygrid(pressure_ratio_map) == table_ygrid(eta_map) || error("pressure_ratio_map/eta_map y grids must match")
    xgrid = table_xgrid(pressure_ratio_map)
    length(flow_min) == length(xgrid) || error("flow_min length must match speed grid length")
    length(flow_max) == length(xgrid) || error("flow_max length must match speed grid length")
    flow_min_f = Float64.(flow_min)
    flow_max_f = Float64.(flow_max)
    all(flow_min_f .<= flow_max_f) || error("flow_min must be <= flow_max at every speed grid point")
    return TabulatedPerformanceMap(
        Float64(Tt_ref),
        Float64(Pt_ref),
        pressure_ratio_map,
        eta_map,
        flow_min_f,
        flow_max_f,
    )
end

function TabulatedPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    pressure_ratio_map::M,
    eta_map::M,
) where {M<:AbstractTableMap}
    ygrid = table_ygrid(pressure_ratio_map)
    xgrid = table_xgrid(pressure_ratio_map)
    return TabulatedPerformanceMap(
        Tt_ref,
        Pt_ref,
        pressure_ratio_map,
        eta_map,
        fill(Float64(first(ygrid)), length(xgrid)),
        fill(Float64(last(ygrid)), length(xgrid)),
    )
end

function TabulatedPerformanceMap(
    Tt_ref::Real,
    Pt_ref::Real,
    speed_grid::Vector{<:Real},
    flow_grid::Vector{<:Real},
    pressure_ratio_table::Matrix{<:Real},
    eta_table::Matrix{<:Real},
    ;
    interpolation::Symbol,
    flow_min::Union{Nothing,Vector{<:Real}}=nothing,
    flow_max::Union{Nothing,Vector{<:Real}}=nothing,
)
    Tt_ref > 0 || error("Tt_ref must be > 0")
    Pt_ref > 0 || error("Pt_ref must be > 0")
    length(speed_grid) >= 2 || error("speed_grid must have at least 2 points")
    length(flow_grid) >= 2 || error("flow_grid must have at least 2 points")
    issorted(speed_grid) || error("speed_grid must be sorted ascending")
    issorted(flow_grid) || error("flow_grid must be sorted ascending")
    size(pressure_ratio_table) == (length(speed_grid), length(flow_grid)) ||
        error("pressure_ratio_table size must match (length(speed_grid), length(flow_grid))")
    size(eta_table) == (length(speed_grid), length(flow_grid)) ||
        error("eta_table size must match (length(speed_grid), length(flow_grid))")

    speed_grid_f = Float64.(speed_grid)
    flow_grid_f = Float64.(flow_grid)
    pressure_ratio_table_f = Float64.(pressure_ratio_table)
    eta_table_f = Float64.(eta_table)
    pressure_ratio_map = interpolation_map(interpolation, speed_grid_f, flow_grid_f, pressure_ratio_table_f)
    eta_map = interpolation_map(interpolation, speed_grid_f, flow_grid_f, eta_table_f)

    flow_min_f = isnothing(flow_min) ? fill(first(flow_grid_f), length(speed_grid_f)) : Float64.(flow_min)
    flow_max_f = isnothing(flow_max) ? fill(last(flow_grid_f), length(speed_grid_f)) : Float64.(flow_max)

    return TabulatedPerformanceMap(Tt_ref, Pt_ref, pressure_ratio_map, eta_map, flow_min_f, flow_max_f)
end

tabulated_speed_grid(map::TabulatedPerformanceMap) = table_xgrid(map.pressure_ratio_map)
tabulated_flow_grid(map::TabulatedPerformanceMap) = table_ygrid(map.pressure_ratio_map)
tabulated_pressure_ratio_table(map::TabulatedPerformanceMap) = table_values(map.pressure_ratio_map)
tabulated_eta_table(map::TabulatedPerformanceMap) = table_values(map.eta_map)
tabulated_interpolation_kind(map::TabulatedPerformanceMap) = table_interpolation(map.pressure_ratio_map)

function tabulated_evaluate(
    map::TabulatedPerformanceMap,
    speed_coord::Real,
    flow_coord::Real,
)
    pressure_ratio = table_evaluate(map.pressure_ratio_map, speed_coord, flow_coord)
    eta = table_evaluate(map.eta_map, speed_coord, flow_coord)
    flow_min = linear_evaluate(tabulated_speed_grid(map), map.flow_min, speed_coord)
    flow_max = linear_evaluate(tabulated_speed_grid(map), map.flow_max, speed_coord)
    return (
        pressure_ratio=pressure_ratio,
        eta=eta,
        low_flow=(flow_coord < flow_min),
        high_flow=(flow_coord > flow_max),
        valid=(flow_min <= flow_coord <= flow_max),
    )
end

function tabulated_domain(map::TabulatedPerformanceMap)
    return (
        speed=(first(tabulated_speed_grid(map)), last(tabulated_speed_grid(map))),
        flow=(minimum(map.flow_min), maximum(map.flow_max)),
        flow_range=(
            min=(speed -> linear_evaluate(tabulated_speed_grid(map), map.flow_min, speed)),
            max=(speed -> linear_evaluate(tabulated_speed_grid(map), map.flow_max, speed)),
        ),
    )
end

"""
Evaluate a shared tabulated performance map from physical operating conditions.

Internally the map converts `(omega, mdot)` to corrected coordinates before
performing the table lookup.
"""
function performance_from_stagnation(
    map::TabulatedPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
)
    speed_coord = corrected_speed(omega, Tt_in, map.Tt_ref)
    flow_coord = corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
    vals = tabulated_evaluate(map, speed_coord, flow_coord)
    return (
        pressure_ratio=vals.pressure_ratio,
        eta=vals.eta,
        speed_coord=speed_coord,
        flow_coord=flow_coord,
        low_flow=vals.low_flow,
        high_flow=vals.high_flow,
        valid=vals.valid,
    )
end

"""
Physical operating domain for a shared tabulated performance map at a given
inlet stagnation state.

The returned bounds are in physical `omega` and physical `mdot`, even though
the stored tables are on corrected coordinates internally.
"""
function performance_map_domain(
    map::TabulatedPerformanceMap,
    Tt_in::Real,
    Pt_in::Real,
)
    domain = tabulated_domain(map)
    omega_scale = sqrt(Tt_in / map.Tt_ref)
    omega_vals = Float64[]
    mdot_vals = Float64[]
    for omega_corr in tabulated_speed_grid(map)
        omega = physical_speed(omega_corr, Tt_in, map.Tt_ref)
        mdot_min_corr = domain.flow_range.min(omega_corr)
        mdot_max_corr = domain.flow_range.max(omega_corr)
        push!(omega_vals, omega)
        push!(mdot_vals, physical_flow(mdot_min_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref))
        push!(mdot_vals, physical_flow(mdot_max_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref))
    end
    return (
        omega=(minimum(omega_vals), maximum(omega_vals)),
        mdot=(minimum(mdot_vals), maximum(mdot_vals)),
        mdot_flow_range=(
            min=(omega -> begin
                omega_corr = corrected_speed(omega, Tt_in, map.Tt_ref)
                mdot_corr = domain.flow_range.min(omega_corr)
                physical_flow(mdot_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
            end),
            max=(omega -> begin
                omega_corr = corrected_speed(omega, Tt_in, map.Tt_ref)
                mdot_corr = domain.flow_range.max(omega_corr)
                physical_flow(mdot_corr, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
            end),
        ),
    )
end

function nearest_grid_indices(
    map::TabulatedPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real;
    n::Int=1,
)
    n >= 1 || error("n must be >= 1")
    speed_coord = corrected_speed(omega, Tt_in, map.Tt_ref)
    flow_coord = corrected_flow(mdot, Tt_in, Pt_in, map.Tt_ref, map.Pt_ref)
    xgrid = tabulated_speed_grid(map)
    ygrid = tabulated_flow_grid(map)
    xspan = max(last(xgrid) - first(xgrid), eps(Float64))
    yspan = max(last(ygrid) - first(ygrid), eps(Float64))
    pairs = Tuple{Float64,CartesianIndex{2}}[]
    for i in eachindex(xgrid), j in eachindex(ygrid)
        dx = (xgrid[i] - speed_coord) / xspan
        dy = (ygrid[j] - flow_coord) / yspan
        push!(pairs, (dx^2 + dy^2, CartesianIndex(i, j)))
    end
    sort!(pairs; by=first)
    return [idx for (_, idx) in Iterators.take(pairs, min(n, length(pairs)))]
end

nearest_grid_index(
    map::TabulatedPerformanceMap,
    omega::Real,
    mdot::Real,
    Tt_in::Real,
    Pt_in::Real,
) = first(nearest_grid_indices(map, omega, mdot, Tt_in, Pt_in; n=1))

"""Demo common tabulated map using the legacy compressor sample data."""
function demo_tabulated_performance_map_compressor(; interpolation::Symbol=:bilinear)
    return TabulatedPerformanceMap(
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
        flow_min=[12.0, 12.0, 12.0],
        flow_max=[20.0, 20.0, 20.0],
    )
end

"""Demo common tabulated map using the legacy turbine sample data."""
function demo_tabulated_performance_map_turbine(; interpolation::Symbol=:bilinear)
    return TabulatedPerformanceMap(
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
        flow_min=[1.4, 1.4, 1.4],
        flow_max=[2.2, 2.2, 2.2],
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
    map::TabulatedPerformanceMap,
    path::AbstractString;
    group::AbstractString="performance_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    node["format"] = "tabulated_performance_map"
    node["format_version"] = 1
    node["Tt_ref"] = map.Tt_ref
    node["Pt_ref"] = map.Pt_ref
    node["flow_min"] = Float64.(map.flow_min)
    node["flow_max"] = Float64.(map.flow_max)
    node["pressure_ratio_map"] = _table_map_to_toml_dict(map.pressure_ratio_map)
    node["eta_map"] = _table_map_to_toml_dict(map.eta_map)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{TabulatedPerformanceMap},
    path::AbstractString;
    group::AbstractString="performance_map",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    haskey(node, "Tt_ref") || error("missing TOML key Tt_ref")
    haskey(node, "Pt_ref") || error("missing TOML key Pt_ref")
    haskey(node, "pressure_ratio_map") || error("missing TOML key pressure_ratio_map")
    haskey(node, "eta_map") || error("missing TOML key eta_map")
    haskey(node, "flow_min") || error("missing TOML key flow_min")
    haskey(node, "flow_max") || error("missing TOML key flow_max")
    Tt_ref = Float64(node["Tt_ref"])
    Pt_ref = Float64(node["Pt_ref"])
    flow_min = Float64.(node["flow_min"])
    flow_max = Float64.(node["flow_max"])
    pressure_ratio_map = _table_map_from_toml_dict(node["pressure_ratio_map"])
    eta_map = _table_map_from_toml_dict(node["eta_map"])
    return TabulatedPerformanceMap(Tt_ref, Pt_ref, pressure_ratio_map, eta_map, flow_min, flow_max)
end
