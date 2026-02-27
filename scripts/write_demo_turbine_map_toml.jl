#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Turbine = TurboMachineModel.Physics.Turbomachine.Turbine

function _build_parser()
    settings = ArgParseSettings(
        prog="write_demo_turbine_map_toml.jl",
        description="Write the demo turbine performance map to TOML.",
    )

    @add_arg_table! settings begin
        "output_path"
            help = "output TOML path"
            required = true
        "--interpolation"
            help = "interpolation type: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--group"
            help = "TOML table/group name"
            arg_type = String
            default = "turbine_map"
    end

    return settings
end

function _main(args::Vector{String})
    parsed = parse_args(args, _build_parser())
    output_path = parsed["output_path"]
    interpolation_raw = parsed["interpolation"]
    interpolation = Symbol(lowercase(interpolation_raw))
    group = parsed["group"]

    interpolation in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(interpolation_raw) (expected bilinear|bicubic)")

    map = Turbine.demo_turbine_performance_map(; interpolation=interpolation)
    Turbine.write_toml(map, output_path; group=group)
    println("Wrote demo turbine map TOML: $(output_path) [group=$(group), interpolation=$(interpolation)]")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
