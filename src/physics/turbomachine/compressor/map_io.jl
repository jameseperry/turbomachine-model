"""
Compressor performance-map TOML loader helpers.
"""

using TOML

const _DEFAULT_COMPRESSOR_MAP_TOML_GROUPS = (
    "compressor_map",
    "compressor_analytic_map",
)

function _toml_group(data::Dict{String,Any}, group::AbstractString)
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

function _compressor_map_type_from_format(format::AbstractString)
    if format == "compressor_performance_map"
        return TabulatedCompressorPerformanceMap
    elseif format == "compressor_nondimensional_performance_map"
        return NonDimensionalTabulatedCompressorPerformanceMap
    elseif format == "compressor_analytic_performance_map"
        return AnalyticCompressorPerformanceMap
    end
    error("unsupported compressor map format $(format)")
end

"""
    read_performance_map_toml(path; group=nothing) -> AbstractCompressorPerformanceMap

Read a compressor performance-map TOML file and dispatch to the concrete map type
based on the TOML `format` field.

If `group` is omitted (`nothing`), this function probes default groups:
- `compressor_map`
- `compressor_analytic_map`
"""
function read_performance_map_toml(
    path::AbstractString;
    group::Union{Nothing,AbstractString}=nothing,
)
    data = TOML.parsefile(path)
    groups = isnothing(group) ? _DEFAULT_COMPRESSOR_MAP_TOML_GROUPS : (String(group),)

    errors = String[]
    for g in groups
        node = try
            _toml_group(data, g)
        catch e
            push!(errors, "group=$(g): $(sprint(showerror, e))")
            continue
        end

        if !haskey(node, "format")
            push!(errors, "group=$(g): missing TOML key format")
            continue
        end

        format = String(node["format"])
        map_type = try
            _compressor_map_type_from_format(format)
        catch e
            push!(errors, "group=$(g): $(sprint(showerror, e))")
            continue
        end

        return read_toml(map_type, path; group=g)
    end

    detail = isempty(errors) ? "" : "\n  " * join(errors, "\n  ")
    error("failed to load compressor performance map from $(path).$(detail)")
end

function read_toml(
    ::Type{AbstractCompressorPerformanceMap},
    path::AbstractString;
    group::Union{Nothing,AbstractString}=nothing,
)
    return read_performance_map_toml(path; group=group)
end
