#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Turbomachine = TurboMachineModel.Physics.Turbomachine
const AxialMachine = TurboMachineModel.Physics.Turbomachine.AxialMachine

function _parse_object_type(raw::AbstractString)
    key = lowercase(strip(raw))
    if key in ("compressor_tabulated", "compressor")
        return :compressor_tabulated
    elseif key in ("axial_compressor_model", "axial_compressor", "compressor_meanline_model", "compressor_meanline", "meanline")
        return :axial_compressor_model
    elseif key in ("axial_turbine_model", "axial_turbine", "compressor_meanline_turbine_model", "compressor_turbine_meanline", "meanline_turbine")
        return :axial_turbine_model
    elseif key in ("turbine_tabulated", "turbine")
        return :turbine_tabulated
    end
    error(
        "unsupported object_type=$(raw) (expected compressor_tabulated|axial_compressor_model|axial_turbine_model|turbine_tabulated)",
    )
end

function _default_group(object_type::Symbol)
    object_type == :compressor_tabulated && return "performance_map"
    object_type == :axial_compressor_model && return "axial_model"
    object_type == :axial_turbine_model && return "axial_model"
    object_type == :turbine_tabulated && return "performance_map"
    error("unsupported object_type=$(object_type)")
end

function _build_demo_object(object_type::Symbol, interpolation::Symbol)
    if object_type == :compressor_tabulated
        return Turbomachine.demo_tabulated_performance_map_compressor(; interpolation=interpolation)
    elseif object_type == :axial_compressor_model
        return Turbomachine.demo_axial_compressor_model()
    elseif object_type == :axial_turbine_model
        return Turbomachine.demo_axial_turbine_model()
    elseif object_type == :turbine_tabulated
        return Turbomachine.demo_tabulated_performance_map_turbine(; interpolation=interpolation)
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
            help = "demo object type: compressor_tabulated, axial_compressor_model, axial_turbine_model, or turbine_tabulated"
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

    if object_type == :compressor_tabulated || object_type == :turbine_tabulated
        Turbomachine.write_toml(obj, output_path; group=group)
    elseif object_type in (:axial_compressor_model, :axial_turbine_model)
        AxialMachine.write_toml(obj, output_path; group=group)
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
