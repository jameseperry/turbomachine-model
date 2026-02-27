#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor
const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

function _parse_object_type(raw::AbstractString)
    key = lowercase(strip(raw))
    if key in ("compressor_tabulated", "compressor")
        return :compressor_tabulated
    elseif key in ("compressor_analytic", "analytic_compressor")
        return :compressor_analytic
    elseif key in ("compressor_spec", "spec")
        return :compressor_spec
    elseif key in ("compressor_design", "design")
        return :compressor_design
    elseif key in ("turbine_tabulated", "turbine")
        return :turbine_tabulated
    end
    error(
        "unsupported object_type=$(raw) (expected compressor_tabulated|compressor_analytic|compressor_spec|compressor_design|turbine_tabulated)",
    )
end

function _default_group(object_type::Symbol)
    object_type == :compressor_tabulated && return "compressor_map"
    object_type == :compressor_analytic && return "compressor_analytic_map"
    object_type == :compressor_spec && return "compressor_spec"
    object_type == :compressor_design && return "compressor_design"
    object_type == :turbine_tabulated && return "turbine_map"
    error("unsupported object_type=$(object_type)")
end

function _build_demo_object(object_type::Symbol, interpolation::Symbol)
    if object_type == :compressor_tabulated
        return Compressor.demo_tabulated_compressor_performance_map(; interpolation=interpolation)
    elseif object_type == :compressor_analytic
        return Compressor.demo_analytic_compressor_performance_map()
    elseif object_type == :compressor_spec
        return Compressor.demo_compressor_spec()
    elseif object_type == :compressor_design
        return Compressor.demo_compressor_design()
    elseif object_type == :turbine_tabulated
        return Turbine.demo_tabulated_turbine_performance_map(; interpolation=interpolation)
    end
    error("unsupported object_type=$(object_type)")
end

function _build_parser()
    settings = ArgParseSettings(
        prog="create_demo_object.jl",
        description="Create a demo object and write it to TOML.",
    )

    @add_arg_table! settings begin
        "object_type"
            help = "demo object type: compressor_tabulated, compressor_analytic, compressor_spec, compressor_design, or turbine_tabulated"
            required = true
        "output_path"
            help = "output TOML path"
            required = true
        "--interpolation"
            help = "interpolation type: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--group"
            help = "TOML table/group name (default depends on object type)"
            arg_type = String
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    object_type = _parse_object_type(parsed["object_type"])
    output_path = parsed["output_path"]

    interpolation_raw = parsed["interpolation"]
    interpolation = Symbol(lowercase(interpolation_raw))
    interpolation in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(interpolation_raw) (expected bilinear|bicubic)")

    group = isnothing(parsed["group"]) ? _default_group(object_type) : parsed["group"]
    obj = _build_demo_object(object_type, interpolation)

    if object_type in (:compressor_tabulated, :compressor_analytic, :compressor_spec, :compressor_design)
        Compressor.write_toml(obj, output_path; group=group)
    elseif object_type == :turbine_tabulated
        Turbine.write_toml(obj, output_path; group=group)
    else
        error("unsupported object_type=$(object_type)")
    end

    println(
        "Wrote demo object: type=$(object_type), path=$(output_path), group=$(group), interpolation=$(interpolation)",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
