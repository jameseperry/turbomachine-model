using TOML

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
        "format" => "table_map",
        "format_version" => 1,
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
    map::AbstractTableMap,
    path::AbstractString;
    group::AbstractString="table_map",
)
    data = Dict{String,Any}()
    node = _find_or_create_group!(data, group)
    merge!(node, _table_map_to_toml_dict(map))
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function read_toml(
    ::Type{AbstractTableMap},
    path::AbstractString;
    group::AbstractString="table_map",
)
    data = TOML.parsefile(path)
    node = _find_group(data, group)
    return _table_map_from_toml_dict(node)
end
