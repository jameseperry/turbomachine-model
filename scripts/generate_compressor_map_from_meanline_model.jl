#!/usr/bin/env julia

using ArgParse
using TurboMachineModel

const Compressor = TurboMachineModel.Physics.Turbomachine.Compressor

function _parse_output_format(raw::AbstractString)
    key = lowercase(strip(raw))
    if key in ("nondimensional", "nd")
        return :nondimensional
    elseif key in ("dimensional", "dim")
        return :dimensional
    end
    error("unsupported output format=$(raw) (expected nondimensional|dimensional)")
end

function _parse_interpolation(raw::AbstractString)
    key = Symbol(lowercase(strip(raw)))
    key in (:bilinear, :bicubic) ||
        error("unsupported interpolation=$(raw) (expected bilinear|bicubic)")
    return key
end

function _parsed_opt(parsed::Dict{String,Any}, primary::String, fallback::String)
    if haskey(parsed, primary)
        return parsed[primary]
    end
    return get(parsed, fallback, nothing)
end

function _build_parser()
    settings = ArgParseSettings(
        prog="generate_compressor_map_from_meanline_model.jl",
        description="Generate a tabulated compressor performance map from a compressor meanline model TOML.",
    )

    @add_arg_table! settings begin
        "input_path"
            help = "input meanline model TOML path"
            required = true
        "output_path"
            help = "output compressor map TOML path"
            required = true
        "--input-group"
            help = "input TOML group/table name"
            arg_type = String
            default = "compressor_meanline_model"
        "--output-group"
            help = "output TOML group/table name"
            arg_type = String
            default = "compressor_map"
        "--output-format"
            help = "output map format: nondimensional or dimensional"
            arg_type = String
            default = "nondimensional"
        "--interpolation"
            help = "tabulated interpolation type: bilinear or bicubic"
            arg_type = String
            default = "bilinear"
        "--n-speed"
            help = "number of speed grid points for tabulation"
            arg_type = Int
            default = 31
        "--n-flow"
            help = "number of flow grid points for tabulation"
            arg_type = Int
            default = 41
        "--boundary-resolution"
            help = "phi probe count per speed line for surge/choke boundary detection"
            arg_type = Int
            default = 401
        "--Tt-in-ref"
            help = "reference inlet total temperature [K] for tabulation and conversion"
            arg_type = Float64
            default = 288.15
        "--Pt-in-ref"
            help = "reference inlet total pressure [Pa] for tabulation and conversion"
            arg_type = Float64
            default = 101_325.0
        "--Tt-ref"
            help = "output corrected-flow reference Tt [K] (dimensional output only)"
            arg_type = Float64
        "--Pt-ref"
            help = "output corrected-flow reference Pt [Pa] (dimensional output only)"
            arg_type = Float64
    end

    return settings
end

function _main(args::Vector{String}=ARGS)
    parsed = parse_args(args, _build_parser())
    input_path = parsed["input_path"]
    output_path = parsed["output_path"]
    input_group = something(_parsed_opt(parsed, "input_group", "input-group"), "compressor_meanline_model")
    output_group = something(_parsed_opt(parsed, "output_group", "output-group"), "compressor_map")
    output_format = _parse_output_format(parsed["output-format"])
    interpolation = _parse_interpolation(parsed["interpolation"])
    n_speed = parsed["n-speed"]
    n_flow = parsed["n-flow"]
    boundary_resolution = parsed["boundary-resolution"]
    Tt_in_ref = parsed["Tt-in-ref"]
    Pt_in_ref = parsed["Pt-in-ref"]
    Tt_ref_opt = _parsed_opt(parsed, "Tt_ref", "Tt-ref")
    Pt_ref_opt = _parsed_opt(parsed, "Pt_ref", "Pt-ref")

    meanline_model = Compressor.read_toml(
        Compressor.CompressorMeanlineModel,
        input_path;
        group=input_group,
    )
    nd_map = Compressor.tabulate_compressor_meanline_model(
        meanline_model;
        Tt_in_ref=Tt_in_ref,
        Pt_in_ref=Pt_in_ref,
        n_speed=n_speed,
        n_flow=n_flow,
        interpolation=interpolation,
        boundary_resolution=boundary_resolution,
    )

    if output_format == :nondimensional
        Compressor.write_toml(nd_map, output_path; group=output_group)
        println(
            "Generated nondimensional compressor map from meanline model: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation), n_speed=$(n_speed), n_flow=$(n_flow)",
        )
        return
    end

    Tt_ref = isnothing(Tt_ref_opt) ? Tt_in_ref : Tt_ref_opt
    Pt_ref = isnothing(Pt_ref_opt) ? Pt_in_ref : Pt_ref_opt
    dim_map = Compressor.to_tabulated_compressor_map(
        nd_map;
        Tt_in_ref=Tt_in_ref,
        Pt_in_ref=Pt_in_ref,
        Tt_ref=Tt_ref,
        Pt_ref=Pt_ref,
        interpolation=interpolation,
    )
    Compressor.write_toml(dim_map, output_path; group=output_group)
    println(
        "Generated dimensional compressor map from meanline model: input=$(input_path) group=$(input_group), output=$(output_path) group=$(output_group), interpolation=$(interpolation), n_speed=$(n_speed), n_flow=$(n_flow), Tt_ref=$(Tt_ref), Pt_ref=$(Pt_ref)",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
