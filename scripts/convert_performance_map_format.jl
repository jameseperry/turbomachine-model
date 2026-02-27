#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor
const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

function _parse_kind(raw::AbstractString)
    kind = Symbol(lowercase(raw))
    kind in (:compressor, :turbine) || error("unsupported kind=$(raw), expected compressor|turbine")
    return kind
end

function _infer_format(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".toml"
        return :toml
    elseif ext == ".h5" || ext == ".hdf5"
        return :hdf5
    end
    error("unsupported extension $(ext) for path $(path); expected .toml/.h5/.hdf5")
end

_default_group(::Val{:compressor}) = "compressor_map"
_default_group(::Val{:turbine}) = "turbine_map"

function _read_map(kind::Symbol, format::Symbol, path::AbstractString, group::AbstractString)
    if kind == :compressor
        if format == :toml
            return Compressor.read_toml(Compressor.TabulatedCompressorPerformanceMap, path; group=group)
        else
            return Compressor.read_hdf5(Compressor.TabulatedCompressorPerformanceMap, path; group=group)
        end
    end

    if format == :toml
        return Turbine.read_toml(Turbine.TabulatedTurbinePerformanceMap, path; group=group)
    else
        return Turbine.read_hdf5(Turbine.TabulatedTurbinePerformanceMap, path; group=group)
    end
end

function _write_map(kind::Symbol, format::Symbol, path::AbstractString, map, group::AbstractString)
    if kind == :compressor
        if format == :toml
            return Compressor.write_toml(map, path; group=group)
        else
            return Compressor.write_hdf5(map, path; group=group)
        end
    end

    if format == :toml
        return Turbine.write_toml(map, path; group=group)
    else
        return Turbine.write_hdf5(map, path; group=group)
    end
end

function _build_parser()
    settings = ArgParseSettings(
        prog="convert_performance_map_format.jl",
        description="Convert compressor/turbine tabulated performance maps between TOML and HDF5.",
    )

    @add_arg_table! settings begin
        "kind"
            help = "map kind: compressor or turbine"
            required = true
        "input_path"
            help = "input map path (.toml/.h5/.hdf5)"
            required = true
        "output_path"
            help = "output map path (.toml/.h5/.hdf5)"
            required = true
        "--input-group"
            help = "input group/table name"
            arg_type = String
        "--output-group"
            help = "output group/table name"
            arg_type = String
    end

    return settings
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _main(args::Vector{String})
    parsed = parse_args(args, _build_parser())
    isnothing(parsed["kind"]) && error("missing required argument kind")
    isnothing(parsed["input_path"]) && error("missing required argument input_path")
    isnothing(parsed["output_path"]) && error("missing required argument output_path")
    kind = _parse_kind(parsed["kind"])
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]
    default_group = _default_group(Val(kind))
    input_group = something(_parsed_opt(parsed, "input_group", "input-group"), default_group)
    output_group = something(_parsed_opt(parsed, "output_group", "output-group"), default_group)
    input_format = _infer_format(input_path)
    output_format = _infer_format(output_path)

    map = _read_map(kind, input_format, input_path, input_group)
    _write_map(kind, output_format, output_path, map, output_group)

    println("Converted $(kind) performance map:")
    println("  input:  $(input_path) [$(input_format), group=$(input_group)]")
    println("  output: $(output_path) [$(output_format), group=$(output_group)]")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
